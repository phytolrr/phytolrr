import os
import copy
import logging
import numpy as np
import rpy2.robjects as ro
from rpy2.robjects import pandas2ri
import rpy2.robjects.numpy2ri
import rpy2.robjects.lib.ggplot2 as ggplot2
import dao
import scripts.paper_calc_data as paper_tools
from rpy2.robjects.packages import importr
import io
import weblogolib
from Bio import SeqIO
from corebio.seq_io import array_io
from corebio import seq as corebio_seq
grdevices = importr('grDevices')

from tools import iqtree as iqtree_tools

OUTPUT_PATH = '../test_files'
MOTIF_VERSION = 3
SUBGROUPS = {'3', '7-1', '7-2', '10', '11', '12', '15'}


def get_color_scheme():
    #cs = copy.deepcopy(weblogolib.colorscheme.hydrophobicity)
    #cs.alphabet = corebio_seq.reduced_protein_alphabet
    weblogolib.colorscheme.hydrophobicity.alphabet = corebio_seq.reduced_protein_alphabet
    return weblogolib.colorscheme.hydrophobicity


def generate_weblogo(file_path, seq_strs, unit_name=None):
    seqs = weblogolib.read_seq_data(io.StringIO('\n'.join(seq_strs)), input_parser=array_io.read, alphabet=corebio_seq.reduced_protein_alphabet)
    data = weblogolib.LogoData.from_seqs(seqs)
    options = weblogolib.LogoOptions(color_scheme=get_color_scheme())
    if unit_name is not None:
        options.unit_name = unit_name
    format = weblogolib.LogoFormat(data, options)
    eps = weblogolib.eps_formatter(data, format)
    with open(file_path, 'wb') as f:
        f.write(eps)


def generate_step1_data(seqs):
    species = set()
    subgroups = set()
    for seq in seqs:
        if seq.subgroup is None:
            continue
        species.add(seq.species)
        subgroups.add(seq.subgroup)
    species_to_index = dict([(sp, i) for i, sp in enumerate(species)])
    subgroups_to_index = dict([(sg, i) for i, sg in enumerate(subgroups)])

    count_matrix = np.zeros((len(subgroups), len(species)), np.int)
    count_species = np.zeros(len(species), np.int)
    for seq in seqs:
        if seq.subgroup is None:
            continue
        species_index = species_to_index[seq.species]
        subgroup_index = subgroups_to_index[seq.subgroup]
        count_matrix[subgroup_index][species_index] += 1
        count_species[species_index] += 1

    probability_matrix = np.zeros((len(subgroups), len(species)), np.float)
    for subgroup_index in range(len(subgroups)):
        for species_index in range(len(species)):
            probability_matrix[subgroup_index][species_index] = float(count_matrix[subgroup_index][species_index]) \
                                                                / float(count_species[species_index])

    return species_to_index, subgroups_to_index, count_matrix, probability_matrix


def generate_step1_file(species_to_index, subgroups_to_index, count_matrix, probability_matrix):
    ro.numpy2ri.activate()
    R = ro.r

    R.library('pheatmap')
    R.png(file=os.path.join(OUTPUT_PATH, 'step1_count_heatmap.png'))
    R.pheatmap(count_matrix)
    R("dev.off()")
    R.png(file=os.path.join(OUTPUT_PATH, 'step1_probability_heatmap.png'))
    R.pheatmap(probability_matrix)
    R("dev.off()")


class Step3P5(object):
    def __init__(self, types):
        self._types = types
        self.subgroups_to_types_to_acc20s = {}
        self.subgroups_to_types_to_seqs_around_nsites = {}
        self.subgroups_to_types_to_n_count = {}
        self.subgroups_to_lrrs_acc20mean = {}
        self.subgroups_to_lrrs_acc20 = {}

    def add_acc20(self, subgroup, t, acc20s):
        self.subgroups_to_types_to_acc20s.setdefault(subgroup, dict([(t, []) for t in self._types]))[t].append(acc20s)

    def set_lrr_acc20s(self, subgroup, lrr_acc20s):
        self.subgroups_to_lrrs_acc20[subgroup] = lrr_acc20s
        self.set_lrr_acc20mean(subgroup, (np.mean(np.array(lrr_acc20s), axis=0).tolist(), len(lrr_acc20s)))

    def set_lrr_acc20mean(self, subgroup, mean):
        self.subgroups_to_lrrs_acc20mean[subgroup] = mean

    def add_seq_str_around_n(self, subgroup, t, seq_str):
        self.subgroups_to_types_to_seqs_around_nsites.setdefault(subgroup, dict([(t, []) for t in self._types]))[t].append(seq_str)

    def add_n_count(self, subgroup, t, count):
        self.subgroups_to_types_to_n_count.setdefault(subgroup, dict([(t, 0) for t in self._types]))[t] += count


def generate_step3_subgroup_type_data(seqs, seq_ids_to_nsites, seq_ids_to_lrrs, all_types, get_n_type, min_lrr_count=20):
    n_data = Step3P5(all_types)
    none_n_data = Step3P5(all_types)
    none_n_without9_data = Step3P5(all_types)
    subgroups_to_seqs = {}
    for seq in seqs:
        if seq.subgroup is None:
            continue
        if seq.subgroup not in subgroups_to_seqs:
            subgroups_to_seqs[seq.subgroup] = []
        subgroups_to_seqs[seq.subgroup].append(seq)

    subgroups_to_filtered_seqs = {}
    for subgroup, seqs in subgroups_to_seqs.items():
        acc20s = []
        subgroups_to_filtered_seqs[subgroup] = []
        less_than_20_lrr_count = 0
        for seq in seqs:
            lrrs = seq_ids_to_lrrs.get(seq.seq_id, [])
            acc20_list = [int(aa_acc20) for aa_acc20 in seq.acc20.split()]

            # 跳过lrr数量小于20的序列
            if len(lrrs) <= min_lrr_count:
                less_than_20_lrr_count += 1
                logging.debug(str.format("Skip sequence {} subgroup {} because the LRR count is {}",
                                         seq.seq_id, subgroup, len(seq_ids_to_lrrs.get(seq.seq_id, []))))
                continue
            subgroups_to_filtered_seqs[subgroup].append(seq.seq_id)

            # 计算取出所有LRR所在位点的acc20
            lrr_9_offsets = set()
            for lrr in lrrs:
                acc20 = acc20_list[lrr.offset:lrr.offset+24]
                if len(acc20) < 24:
                    logging.warning(str.format("Ignore LRR acc20 in seq id {} at offset {}, total length {}",
                                    seq.seq_id, lrr.offset, len(acc20_list)))
                    continue
                acc20s.append(acc20)
                lrr_9_offsets.add(lrr.offset + 8)
            logging.debug(str.format("LRR offsets for seq {}, {}", seq.seq_id, lrr_9_offsets))

            # 计算nsites位点前后各5个，附近的acc20、序列
            nsites = seq_ids_to_nsites.get(seq.seq_id)
            if len(nsites) == 0:
                continue
            nsite_poses = set([nsite.start_pos for nsite in nsites])

            for index, amino in enumerate(seq.seq):
                if seq.seq[index] != 'N':
                    continue
                acc20_around_nsites = paper_tools.get_by_start_pos(acc20_list, index)
                seq_around_nsites = paper_tools.get_by_start_pos(seq.seq, index)
                if acc20_around_nsites is None or seq_around_nsites is None:
                    assert(acc20_around_nsites is None)
                    assert(seq_around_nsites is None)
                    logging.warning(str.format("Drop acc20 and seq at start pos {} for seq {}, can not get 11 amino around the pos",
                                               index, seq.seq_id))
                    continue
                type_ = get_n_type(seq, index)
                if type_ is None:
                    logging.warning(str.format("The type for index {} seq {} is None, skip it", index, seq.seq_id))
                    continue
                if index in nsite_poses:
                    n_data.add_acc20(subgroup, type_, acc20_around_nsites)
                    n_data.add_seq_str_around_n(subgroup, type_, seq_around_nsites)
                    n_data.add_n_count(subgroup, type_, 1)
                elif index not in lrr_9_offsets:
                    none_n_data.add_acc20(subgroup, type_, acc20_around_nsites)
                    none_n_data.add_seq_str_around_n(subgroup, type_, seq_around_nsites)
                    none_n_data.add_n_count(subgroup, type_, 1)
                    none_n_without9_data.add_acc20(subgroup, type_, acc20_around_nsites)
                    none_n_without9_data.add_seq_str_around_n(subgroup, type_, seq_around_nsites)
                    none_n_without9_data.add_n_count(subgroup, type_, 1)
                else:
                    none_n_data.add_acc20(subgroup, type_, acc20_around_nsites)
                    none_n_data.add_seq_str_around_n(subgroup, type_, seq_around_nsites)
                    none_n_data.add_n_count(subgroup, type_, 1)

        if len(acc20s) == 0:
            logging.warning("The subgroup {} does not have valid acc20s, skipped/total seq count {}/{}".format(
                subgroup, less_than_20_lrr_count, len(seqs)))
            continue
        n_data.set_lrr_acc20s(subgroup, acc20s)

    for subgroup, seq_ids in subgroups_to_filtered_seqs.items():
        logging.info(str.format("Subgroup {} seq num {}, {}", subgroup, len(seq_ids), seq_ids))

    return n_data, none_n_data, none_n_without9_data


def generate_ste3_5_file_weblogo(subgroups_to_sses_to_seqs_around_nsites, file_prefix):
    total_sses_to_seqs = {}
    for subgroup, sses_to_seqs_around_nsites in subgroups_to_sses_to_seqs_around_nsites.items():
        for ss, seqs_around_nsites in sses_to_seqs_around_nsites.items():
            if ss not in total_sses_to_seqs:
                total_sses_to_seqs[ss] = []
            total_sses_to_seqs[ss] += seqs_around_nsites
            file_path = os.path.join(OUTPUT_PATH, str.format("{}_{}_{}.eps", file_prefix, subgroup, ss))
            if len(seqs_around_nsites) < 2:
                logging.warning("The size seqs_around_nsites for subgroup {} type ss type {} is {}, can not generate weblogo for it",
                                subgroup, ss, len(seqs_around_nsites))
                open(file_path, 'wb').close()
            else:
                generate_weblogo(file_path, seqs_around_nsites)

    for ss, seqs_around_nsites in total_sses_to_seqs.items():
        file_path = os.path.join(OUTPUT_PATH, str.format("{}_{}_{}.eps", file_prefix, "all", ss))
        generate_weblogo(file_path, seqs_around_nsites)
        file_path = os.path.join(OUTPUT_PATH, str.format("{}_{}_{}.eps", file_prefix, "all_probability", ss))
        generate_weblogo(file_path, seqs_around_nsites, unit_name="probability")


def _generate_step3_5_ss_acc20_line_chart(ts_to_acc20s, tname, line_chart_file_path):
    logging.debug(str.format("Begin to generate {}, data {}", line_chart_file_path, ts_to_acc20s))
    ts_to_acc20mean = calc_acc20mean_by_types(ts_to_acc20s)
    columns_to_data = {tname: [], 'site': [], 'acc20': []}
    for ss, acc20means in ts_to_acc20mean.items():
        for index, acc20mean in enumerate(acc20means):
            columns_to_data[tname].append(ss)
            columns_to_data['site'].append(index - 5)
            columns_to_data['acc20'].append(acc20mean)

    # Generate the line chart file
    r_columns_to_data = {
        tname: ro.StrVector(columns_to_data[tname]),
        'site': ro.IntVector(columns_to_data['site']),
        'acc20': ro.FloatVector(columns_to_data['acc20'])
    }
    df = ro.DataFrame(r_columns_to_data)

    logging.debug(str.format("The Data Frame for file {}: \n{}", line_chart_file_path, df))
    grdevices.png(file=line_chart_file_path, width=1024, height=512)
    gp = ggplot2.ggplot(df)
    pp = gp + \
         ggplot2.theme_bw() + \
         ggplot2.theme_classic() + \
         ggplot2.theme(**{'axis.text.x': ggplot2.element_text(size=35)}) + \
         ggplot2.theme(**{'axis.text.y': ggplot2.element_text(size=35)}) + \
         ggplot2.aes_string(x='site', y='acc20', group=tname, colour=tname) + \
         ggplot2.geom_point(size=4, shape=20) + \
         ggplot2.geom_line(size=3) + \
         ggplot2.theme(**{'legend.title': ggplot2.element_blank()}) + \
         ggplot2.theme(**{'legend.text': ggplot2.element_text(size=20)}) + \
         ggplot2.scale_x_continuous(breaks=ro.IntVector(list(range(-5, 6))),
                                    labels=ro.StrVector(['-5', '-4', '-3', '-2', '-1', 'N', '1', '2', '3', '4', '5']))
    pp.plot()
    logging.info(str.format("Output step3 file {}", line_chart_file_path))
    grdevices.dev_off()


def calc_acc20mean_by_types(types_to_acc20s):
    types_to_acc20mean = {}
    for s, acc20s in types_to_acc20s.items():
        if len(acc20s) != 0:
            acc20s_mean = np.mean(np.array(acc20s), axis=0).tolist()
            types_to_acc20mean[s] = acc20s_mean
        else:
            types_to_acc20mean[s] = []
    return types_to_acc20mean


def merge_subgroup(subgroups_to_types_to_value):
    total_types_to_value = {}
    for _, types_to_value in subgroups_to_types_to_value.items():
        for t, value in types_to_value.items():
            if t not in total_types_to_value:
                total_types_to_value[t] = copy.copy(value)
            else:
                total_types_to_value[t] += value
    return total_types_to_value


def generate_step3_5_ss_acc20_line_chart(subgroups_to_types_to_acc20s, tname, file_prefix):
    for subgroup, ts_to_acc20s in subgroups_to_types_to_acc20s.items():
        line_chart_file_path = os.path.join(OUTPUT_PATH, str.format("{}_{}.png", file_prefix, subgroup))
        _generate_step3_5_ss_acc20_line_chart(ts_to_acc20s, tname, line_chart_file_path)

    total_ts_to_acc20s = merge_subgroup(subgroups_to_types_to_acc20s)
    line_chart_file_path = os.path.join(OUTPUT_PATH, str.format("{}_{}.png", file_prefix, "all"))
    _generate_step3_5_ss_acc20_line_chart(total_ts_to_acc20s, tname, line_chart_file_path)


def generate_step3_5_lrr_acc20_line_chart(subgroups_to_lrrs_acc20mean, prefix=''):
    pandas2ri.activate()
    subgroups_to_lrr_count = {}
    columns_to_data = {'subgroup': [], 'pos': [], 'acc20': []}
    for subgroup, (acc20means, acc20_count) in subgroups_to_lrrs_acc20mean.items():
        subgroups_to_lrr_count[subgroup] = acc20_count
        for index, acc20mean in enumerate(acc20means):
            columns_to_data['subgroup'].append(subgroup)
            columns_to_data['pos'].append(index + 1)
            columns_to_data['acc20'].append(acc20mean)

    # Write the count of LRRs for each subgroup to file
    with open(os.path.join(OUTPUT_PATH, prefix + "step3_5_lrr_count.txt"), 'w') as f:
        for subgroup, lrr_count in subgroups_to_lrr_count.items():
            f.write(str.format("{}: {}\n", subgroup, lrr_count))

    # Generate the line chart file
    r_columns_to_data = {
        'subgroup': ro.StrVector(columns_to_data['subgroup']),
        'pos': ro.IntVector(columns_to_data['pos']),
        'acc20': ro.FloatVector(columns_to_data['acc20'])
    }
    df = ro.DataFrame(r_columns_to_data)

    line_chart_file_path = os.path.join(OUTPUT_PATH, prefix + "step3_5_lrr_acc20_line.png")
    logging.debug(str.format("The Data Frame for file {}: \n{}",line_chart_file_path, df))
    grdevices.png(file=line_chart_file_path, width=1024, height=512)
    gp = ggplot2.ggplot(df)
    pp = gp + \
         ggplot2.theme_bw() + \
         ggplot2.theme_classic() + \
         ggplot2.theme(**{'axis.text.x': ggplot2.element_text(size=35)}) + \
         ggplot2.theme(**{'axis.text.y': ggplot2.element_text(size=35)}) + \
         ggplot2.aes_string(x='pos', y='acc20', group='subgroup', colour='subgroup') + \
         ggplot2.geom_point(size=4, shape=20) + \
         ggplot2.geom_line(size=3) + \
         ggplot2.theme(**{'legend.title': ggplot2.element_blank()}) + \
         ggplot2.theme(**{'legend.text': ggplot2.element_text(size=20)}) + \
         ggplot2.scale_x_continuous(breaks=ro.IntVector(range(1, 25)), labels=ro.StrVector(list('LxxLxLxxNxLsGxIPxxLxxLxx')))
    pp.plot()
    logging.info(str.format("Output step3 file {}", line_chart_file_path))
    grdevices.dev_off()


def generate_sses_to_ncount_histogram(sses_to_n_count, file_name):
    columns_to_data = {
        'ss': [],
        'count': []
    }
    max_count = 0
    for ss, n_count in sses_to_n_count.items():
        columns_to_data['ss'].append(ss)
        columns_to_data['count'].append(n_count)
        if n_count > max_count:
            max_count = n_count
    r_columns_to_data = {
        'ss': ro.StrVector(columns_to_data['ss']),
        'count': ro.IntVector(columns_to_data['count'])
    }
    df = ro.DataFrame(r_columns_to_data)

    max_count = int(max_count / 1000 * 1000 * 1.5)
    histogram_file_path = os.path.join(OUTPUT_PATH, file_name)
    logging.debug(str.format("The Data Frame for file {}: \n{}",histogram_file_path, df))
    grdevices.png(file=histogram_file_path, width=1024, height=700)
    gp = ggplot2.ggplot(df)
    pp = gp + \
         ggplot2.aes_string(x='ss', y='count') + \
         ggplot2.geom_bar(position="dodge",width=0.8, stat="identity") + \
         ggplot2.theme_bw() + \
         ggplot2.theme_classic() + \
         ggplot2.theme(**{'legend.title': ggplot2.element_blank()}) + \
         ggplot2.theme(**{'legend.text': ggplot2.element_text(size=20)}) + \
         ggplot2.theme(**{'axis.text.x': ggplot2.element_text(size=40)}) + \
         ggplot2.theme(**{'axis.text.y': ggplot2.element_text(size=40)}) + \
         ggplot2.scale_y_continuous(expand=ro.IntVector([0, 0]),
                                    limits=ro.IntVector([0, max_count])) + \
         ggplot2.geom_text(ggplot2.aes_string(label='count'), size=20, angle=45, hjust=-0.2,
                           position=ggplot2.position_dodge(width=0.8),
                           vjust=-0.2)
    pp.plot()
    logging.info(str.format("Output step3_5 file {}", histogram_file_path))
    grdevices.dev_off()


def _sort_subgroup(subgroups):
    return sorted(subgroups, key=lambda subgroup: [int(num) for num in subgroup.split('-')])


def generate_histogram(subgroups_to_sses_to_n_count, tname, file_name):
    columns_to_data = {
        'subgroup': [],
        tname: [],
        'count': []
    }
    max_count = 0
    for subgroup, sses_to_n_count in subgroups_to_sses_to_n_count.items():
        for ss, n_count in sses_to_n_count.items():
            columns_to_data['subgroup'].append(subgroup)
            columns_to_data[tname].append(ss)
            columns_to_data['count'].append(n_count)
            if n_count > max_count:
                max_count = n_count
    r_columns_to_data = {
        'subgroup': ro.FactorVector(columns_to_data['subgroup'],
                                    levels=ro.StrVector(_sort_subgroup(set(columns_to_data['subgroup'])))),
        tname: ro.StrVector(columns_to_data[tname]),
        'count': ro.IntVector(columns_to_data['count'])
    }
    df = ro.DataFrame(r_columns_to_data)

    max_count = int(max_count / 1000 * 1000 + 1000)
    histogram_file_path = os.path.join(OUTPUT_PATH, file_name)
    logging.debug(str.format("The Data Frame for file {}: \n{}",histogram_file_path, df))

    grdevices.png(file=histogram_file_path, width=1200, height=800)
    gp = ggplot2.ggplot(df)
    pp = gp + \
         ggplot2.aes_string(x='subgroup', y='count', fill=tname) + \
         ggplot2.geom_bar(position="dodge",width=0.8, stat="identity") + \
         ggplot2.theme_bw() + \
         ggplot2.theme_classic() + \
         ggplot2.theme(**{'legend.title': ggplot2.element_blank()}) + \
         ggplot2.theme(**{'legend.text': ggplot2.element_text(size=40)}) + \
         ggplot2.theme(**{'axis.text.x': ggplot2.element_text(size=40,angle=45)}) + \
         ggplot2.theme(**{'axis.text.y': ggplot2.element_text(size=40)}) + \
         ggplot2.scale_y_continuous(expand=ro.IntVector([0, 0]),
                                    limits=ro.IntVector([0, max_count])) + \
         ggplot2.geom_text(ggplot2.aes_string(label='count'), size=6, angle=35, hjust=-0.1,
                           position=ggplot2.position_dodge(width=0.8),
                           vjust=-0.2)

    pp.plot()
    logging.info(str.format("Output step3 file {}", histogram_file_path))
    grdevices.dev_off()


def generate_step3_9_n_count_histogram(place_type_pos_type_to_count, file_name):
    columns_to_data = {
        'place': [],
        'pos': [],
        'count': []
    }
    max_count = 0
    for place_pos_type, n_count in place_type_pos_type_to_count.items():
        place_type, pos_type = place_pos_type.split('_')
        columns_to_data['place'].append(place_type)
        columns_to_data['pos'].append(pos_type)
        columns_to_data['count'].append(n_count)
        if n_count > max_count:
            max_count = n_count
    r_columns_to_data = {
        'place': ro.StrVector(columns_to_data['place']),
        'pos': ro.StrVector(columns_to_data['pos']),
        'count': ro.IntVector(columns_to_data['count'])
    }
    df = ro.DataFrame(r_columns_to_data)

    if max_count > 1000:
        max_count = int(max_count / 1000 * 1000 + 1000)
    else:
        max_count = int(max_count / 100 * 100 + 100)
    histogram_file_path = os.path.join(OUTPUT_PATH, file_name)
    logging.debug(str.format("The Data Frame for file {}: \n{}",histogram_file_path, df))
    grdevices.png(file=histogram_file_path, width=1024, height=512)
    gp = ggplot2.ggplot(df)
    pp = gp + \
         ggplot2.aes_string(x='pos', y='count', fill='place') + \
         ggplot2.geom_bar(position="dodge", stat="identity") + \
         ggplot2.theme_bw() + \
         ggplot2.theme_classic() + \
         ggplot2.theme(**{'axis.text.x': ggplot2.element_text(size=35)}) + \
         ggplot2.theme(**{'axis.text.y': ggplot2.element_text(size=35)}) + \
         ggplot2.scale_y_continuous(expand=ro.IntVector([0, 0]),
                                    limits=ro.IntVector([0, max_count])) + \
         ggplot2.geom_text(ggplot2.aes_string(label='count'),
                           position=ggplot2.position_dodge(width=0.8), size=10, angle=35, hjust=-0.2,
                           vjust=-0.5)
    pp.plot()
    logging.info(str.format("Output step3 file {}", histogram_file_path))
    grdevices.dev_off()


def get_seq_ids_to_indexes_to_place_type(seqs, seq_ids_to_lrrs):
    seq_ids_to_seq = dict([(seq.seq_id, seq) for seq in seqs])
    seq_ids_to_indexes_to_type = {}
    for seq_id, lrrs in seq_ids_to_lrrs.items():
        seq_len = len(seq_ids_to_seq[seq_id].seq)
        seq_ids_to_indexes_to_type[seq_id] = {}
        for lrr in lrrs:
            if lrr.offset + 24 > seq_len:
                logging.warning(str.format("The LRR at offset {} seq {} will exceeds the seq end({}) after add 24",
                                           lrr.offset, seq_id, seq_len))
                continue
            for index in range(lrr.offset, lrr.offset + 9):
                seq_ids_to_indexes_to_type[seq_id][index] = 'inside'
            for index in range(lrr.offset + 9, lrr.offset + 13):
                seq_ids_to_indexes_to_type[seq_id][index] = 'beside'
            for index in range(lrr.offset + 13, lrr.offset + 24):
                seq_ids_to_indexes_to_type[seq_id][index] = 'back'
    return seq_ids_to_indexes_to_type


def get_seq_ids_to_indexes_to_percent_pos_type(seqs, seq_ids_to_lrrs):
    seq_ids_to_indexes_to_type = {}
    for seq in seqs:
        seq_id = seq.seq_id
        lrrs = seq_ids_to_lrrs.get(seq_id, [])
        lrrs.sort(key=lambda lrr: lrr.offset)
        seq_ids_to_indexes_to_type[seq_id] = {}
        total_count = len(lrrs)
        percent_30 = int(total_count / 3)
        for i in range(0, percent_30):
            seq_ids_to_indexes_to_type[seq_id].update(
                dict([(index, 'top30') for index in range(lrrs[i].offset, lrrs[i].offset + 24)]))

        middle_end = total_count - percent_30
        for i in range(percent_30, middle_end):
            seq_ids_to_indexes_to_type[seq_id].update(
                dict([(index, 'middle') for index in range(lrrs[i].offset, lrrs[i].offset + 24)]))
        for i in range(middle_end, total_count):
            seq_ids_to_indexes_to_type[seq_id].update(
                dict([(index, 'bottom30') for index in range(lrrs[i].offset, lrrs[i].offset + 24)]))

    return seq_ids_to_indexes_to_type


class Step3P10(object):
    def __init__(self):
        self.standard_seq_id = None
        self.seq_ids_to_seq = {}
        self.seq_ids_to_ali_indexes_to_real_index = {}
        self.standard_real_indexes_to_ali_index = None

    def get_ali_indexes_to_real_index(self, seq_id):
        ali_indexes_to_real_index = self.seq_ids_to_ali_indexes_to_real_index.get(seq_id, None)
        if ali_indexes_to_real_index is None:
            ali_indexes_to_real_index = {}
            self.seq_ids_to_ali_indexes_to_real_index[seq_id] = ali_indexes_to_real_index
            real_index = 0
            for ali_index, amino in enumerate(self.seq_ids_to_seq[seq_id]):
                if amino != 'X':
                    ali_indexes_to_real_index[ali_index] = real_index
                    real_index += 1
        return ali_indexes_to_real_index

    def get_standard_real_indexes_to_ali_index(self):
        if self.standard_real_indexes_to_ali_index is None:
            self.standard_real_indexes_to_ali_index = {}
            real_index = 0
            for i, amino in enumerate(self.seq_ids_to_seq[self.standard_seq_id]):
                if amino != 'X':
                    self.standard_real_indexes_to_ali_index[real_index] = i
                    real_index += 1
        return self.standard_real_indexes_to_ali_index


def read_in_ali_seqs(seq_ids_to_file_path):
    standard_seq_ids_to_data = {}
    for seq_id, (seq_file_path, tree_file_path) in seq_ids_to_file_path.items():
        # parse relevant aka colored seq ids from the iqtree file
        tree = iqtree_tools.parse_tree_by_file(tree_file_path)
        colors_to_seq_ids = iqtree_tools.group_by_colors(tree)
        relevant_seq_ids = set([seq_id for seq_ids in colors_to_seq_ids.values() for seq_id in seq_ids])

        seq_ids_to_seq = dict([(seq.id.replace('-', ''), str(seq.seq).replace('-', 'X'))
                               for seq in SeqIO.parse(seq_file_path, "fasta")
                               if seq.id in relevant_seq_ids])
        if seq_id not in seq_ids_to_seq:
            logging.error(str.format("The standard seq {} does not in file {}", seq_id, seq_file_path))
            continue

        seq_data = Step3P10()
        seq_data.standard_seq_id = seq_id
        seq_data.seq_ids_to_seq = seq_ids_to_seq
        standard_seq_ids_to_data[seq_id] = seq_data

    return standard_seq_ids_to_data


def generate_step3_10_weblogo_per_offset(seq_ids_to_data, seq_ids_to_lrrs):
    for seq_id, data in seq_ids_to_data.items():
        real_indexes = []
        for lrr in seq_ids_to_lrrs[seq_id]:
            real_indexes.append(lrr.offset)
        for i in range(0, 24):
            all_picked_seq_str = []
            real_indexes_to_align_index = data.get_standard_real_indexes_to_ali_index()
            indexes = [real_indexes_to_align_index[real_index + i] for real_index in real_indexes]
            for tmp_seq_id, seq_str in data.seq_ids_to_seq.items():
                picked_seq_str = "".join([seq_str[index] for index in indexes])
                all_picked_seq_str.append(picked_seq_str)
                if tmp_seq_id == seq_id:
                    logging.info(str.format("Standard seq {}, picked str {}", seq_id, picked_seq_str))
            logging.info(str.format("LRRs for the {} offset, seq ids {}: {}",
                                    i, data.seq_ids_to_seq.keys(), all_picked_seq_str))
            if seq_id == 'AT4G39400.1' and i == 20:
                print(all_picked_seq_str)
            generate_weblogo(os.path.join(OUTPUT_PATH, str.format("step3_10_weblogo_{}_{}.eps", seq_id, i)),
                             all_picked_seq_str)


def generate_step3_10_weblog_lrrs(seq_ids_to_data, seq_ids_to_seq, seq_ids_to_lrrs):
    for seq_id, data in seq_ids_to_data.items():
        all_lrrs = []

        standard_real_indexes_to_align_index = data.get_standard_real_indexes_to_ali_index()
        for lrr in seq_ids_to_lrrs[seq_id]:
            real_index = lrr.offset
            align_index = standard_real_indexes_to_align_index[real_index]
            for relevant_seq_id in data.seq_ids_to_seq.keys():
                ali_indexes_to_real_index = data.get_ali_indexes_to_real_index(relevant_seq_id)
                real_index = None
                try_times = 0
                while real_index is None and try_times < 6:
                    try_times += 1
                    real_index = ali_indexes_to_real_index.get(align_index, None)
                    if real_index is None:
                        logging.warning(str.format("The align index {}(+{}) is gap in seq {}, standard seq {}",
                                                   align_index, try_times, relevant_seq_id, seq_id))
                if real_index is None:
                    logging.warning(str.format("Drop the align index {} for seq {}, standard seq {} for no real indexes",
                                             align_index, relevant_seq_id, seq_id))
                    continue

                lrr_str = seq_ids_to_seq[relevant_seq_id].seq[real_index:real_index + 24]
                if len(lrr_str) != 24:
                    logging.warning(str.format("Drop the align index {} real index {} for seq {}, standard seq {}, for lrr({}) length is less than 24",
                                             align_index, real_index, relevant_seq_id, seq_id, lrr_str))
                else:
                    all_lrrs.append(lrr_str)

        generate_weblogo(os.path.join(OUTPUT_PATH, str.format("step3_10_weblogo_lrr_ali_{}.eps", seq_id)), all_lrrs)


def generate_step3_10_weblogo_by_types(seq_ids_to_data, seq_ids_to_seq, seq_ids_to_lrrs):
    for seq_id, data in seq_ids_to_data.items():
        all_lrrs = []
        relevant_seq_ids = set(data.seq_ids_to_seq.keys())
        for relevant_seq_id in relevant_seq_ids:
            for lrr in seq_ids_to_lrrs[relevant_seq_id]:
                lrr_seq = seq_ids_to_seq[relevant_seq_id].seq[lrr.offset:lrr.offset+24]
                if len(lrr_seq) != 24:
                    logging.error(str.format("The LRR offset {} in seq {} type {}, is less than 24 {}",
                                             lrr.offset, relevant_seq_id, seq_id, lrr_seq))
                else:
                    all_lrrs.append(lrr_seq)
        logging.info(str.format("LRRs in all seq by type {}: {}", seq_id, all_lrrs))
        generate_weblogo(os.path.join(OUTPUT_PATH, str.format("step3_10_weblogo_lrr_{}.eps", seq_id)), all_lrrs)


def step1(seqs):
    species_to_index, subgroups_to_index, count_matrix, probability_matrix = generate_step1_data(seqs)
    logging.info(str.format("Step1 data generated end, row(species) index: {}, col(subgroups) index: {}",
                            species_to_index, subgroups_to_index))
    logging.info(str.format("Step1 count-matrix: \n{}\n\nprobability-matrix:\n {}",
                            count_matrix, probability_matrix))
    generate_step1_file(species_to_index, subgroups_to_index, count_matrix, probability_matrix)


def step3_5(seqs, seq_ids_to_nsites, seq_ids_to_lrrs):
    data, none_n_data, none_n_without9_data = generate_step3_subgroup_type_data(seqs, seq_ids_to_nsites, seq_ids_to_lrrs,
                                                          list('CHE'), lambda seq, index: seq.ss[index])
    logging.debug(str.format("Generate step3.5 web logo data end {}", data.subgroups_to_types_to_seqs_around_nsites))
    logging.debug(str.format("Generate step3.5 lrr acc20 data end {}", data.subgroups_to_lrrs_acc20mean))
    generate_step3_5_lrr_acc20_line_chart(data.subgroups_to_lrrs_acc20mean)
    generate_step3_5_ss_acc20_line_chart(data.subgroups_to_types_to_acc20s, 'ss', 'step3_5_nsites_line')
    generate_ste3_5_file_weblogo(data.subgroups_to_types_to_seqs_around_nsites, 'step3_5_weblogo')
    generate_histogram(data.subgroups_to_types_to_n_count, 'ss', 'step3_5_n_count_histogram.png')
    generate_sses_to_ncount_histogram(merge_subgroup(data.subgroups_to_types_to_n_count),
                                              'step3_5_n_count_byss_histogram.png')

    generate_step3_5_ss_acc20_line_chart(none_n_data.subgroups_to_types_to_acc20s, 'ss', 'step3_5_nn_nsites_line')
    generate_ste3_5_file_weblogo(none_n_data.subgroups_to_types_to_seqs_around_nsites, 'step3_5_nn_weblogo')
    generate_histogram(none_n_data.subgroups_to_types_to_n_count, 'ss', 'step3_5_nn_n_count_histogram.png')
    generate_sses_to_ncount_histogram(merge_subgroup(none_n_data.subgroups_to_types_to_n_count),
                                              'step3_5_nn_n_count_byss_histogram.png')

    generate_step3_5_ss_acc20_line_chart(none_n_without9_data.subgroups_to_types_to_acc20s, 'ss', 'step3_5_nnwo9_nsites_line')
    generate_ste3_5_file_weblogo(none_n_without9_data.subgroups_to_types_to_seqs_around_nsites, 'step3_5_nnwo9_weblogo')
    generate_histogram(none_n_without9_data.subgroups_to_types_to_n_count, 'ss', 'step3_5_nnwo9_n_count_histogram.png')
    generate_sses_to_ncount_histogram(merge_subgroup(none_n_without9_data.subgroups_to_types_to_n_count),
                                              'step3_5_nnwo9_n_count_byss_histogram.png')

    logging.info("Step3 generated end")


def step3_7(seqs, seq_ids_to_nsites, seq_ids_to_lrrs):
    seq_ids_to_indexes_to_type = get_seq_ids_to_indexes_to_place_type(seqs, seq_ids_to_lrrs)

    data, none_n_data, none_n_without9_data = generate_step3_subgroup_type_data(seqs, seq_ids_to_nsites, seq_ids_to_lrrs,
                                                          ['inside', 'beside', 'back'],
                                                          lambda seq, index: seq_ids_to_indexes_to_type[seq.seq_id].get(index, None))

    generate_ste3_5_file_weblogo(data.subgroups_to_types_to_seqs_around_nsites, 'step3_7_weblogo')
    generate_ste3_5_file_weblogo(none_n_data.subgroups_to_types_to_seqs_around_nsites, 'step3_7_nn_weblogo')
    generate_ste3_5_file_weblogo(none_n_without9_data.subgroups_to_types_to_seqs_around_nsites, 'step3_7_nnwo9_weblogo')

    generate_step3_5_ss_acc20_line_chart(data.subgroups_to_types_to_acc20s, 'place', 'step3_7_nsites_line')
    generate_step3_5_ss_acc20_line_chart(none_n_data.subgroups_to_types_to_acc20s, 'place', 'step3_7_nn_nsites_line')
    generate_step3_5_ss_acc20_line_chart(none_n_without9_data.subgroups_to_types_to_acc20s, 'place', 'step3_7_nnwo9_nsites_line')

    generate_histogram(data.subgroups_to_types_to_n_count, 'place', 'step3_7_n_count_histogram.png')
    generate_histogram(none_n_data.subgroups_to_types_to_n_count, 'place', 'step3_7_nn_n_count_histogram.png')
    generate_histogram(none_n_without9_data.subgroups_to_types_to_n_count, 'place', 'step3_7_nnwo9_n_count_histogram.png')


    generate_sses_to_ncount_histogram(merge_subgroup(data.subgroups_to_types_to_n_count),
                                              'step3_7_n_count_byss_histogram.png')
    generate_sses_to_ncount_histogram(merge_subgroup(none_n_data.subgroups_to_types_to_n_count),
                                              'step3_7_nn_n_count_byss_histogram.png')
    generate_sses_to_ncount_histogram(merge_subgroup(none_n_without9_data.subgroups_to_types_to_n_count),
                                              'step3_7_nnwo9_n_count_byss_histogram.png')

def get_step3_9_type(seq, index, seq_ids_to_indexes_to_place_type, seq_ids_to_indexes_to_percent_pos_type):
    place_type = seq_ids_to_indexes_to_place_type[seq.seq_id].get(index, None)
    pos_type = seq_ids_to_indexes_to_percent_pos_type[seq.seq_id].get(index, None)
    if place_type is None or pos_type is None:
        return None
    return place_type + '_' + pos_type


def step3_9(seqs, seq_ids_to_nsites, seq_ids_to_lrrs):
    seq_ids_to_indexes_to_place_type = get_seq_ids_to_indexes_to_place_type(seqs, seq_ids_to_lrrs)
    seq_ids_to_indexes_to_percent_pos_type = get_seq_ids_to_indexes_to_percent_pos_type(seqs, seq_ids_to_lrrs)

    all_types = []
    for place_type in ['inside', 'beside', 'back']:
        for pos_type in ['top30', 'middle', 'bottom30']:
            all_types.append(place_type + '_' + pos_type)

    data, none_n_data, none_n_without9_data = generate_step3_subgroup_type_data(seqs, seq_ids_to_nsites, seq_ids_to_lrrs, all_types,
                            lambda seq, index: get_step3_9_type(seq, index,
                                                                seq_ids_to_indexes_to_place_type,
                                                                seq_ids_to_indexes_to_percent_pos_type))

    for subgroup, types_to_n_count in data.subgroups_to_types_to_n_count.items():
        generate_step3_9_n_count_histogram(types_to_n_count, 'step3_9_place_pos_ncount_' + subgroup + '_histogram.png')

    total_types_to_n_count = merge_subgroup(data.subgroups_to_types_to_n_count)
    generate_step3_9_n_count_histogram(total_types_to_n_count, 'step3_9_place_pos_ncount_all_histogram.png')

    for subgroup, types_to_n_count in none_n_data.subgroups_to_types_to_n_count.items():
        generate_step3_9_n_count_histogram(types_to_n_count, 'step3_9_nn_place_pos_ncount_' + subgroup + '_histogram.png')

    total_types_to_n_count = merge_subgroup(none_n_data.subgroups_to_types_to_n_count)
    generate_step3_9_n_count_histogram(total_types_to_n_count, 'step3_9_nn_place_pos_ncount_all_histogram.png')

    for subgroup, types_to_n_count in none_n_without9_data.subgroups_to_types_to_n_count.items():
        generate_step3_9_n_count_histogram(types_to_n_count, 'step3_9_nnwo9_place_pos_ncount_' + subgroup + '_histogram.png')

    total_types_to_n_count = merge_subgroup(none_n_without9_data.subgroups_to_types_to_n_count)
    generate_step3_9_n_count_histogram(total_types_to_n_count, 'step3_9_nnwo9_place_pos_ncount_all_histogram.png')


def statitics_match_count(seq_ids_to_data, seq_ids_to_lrrs):
    for standard_seq_id, data in seq_ids_to_data.items():
        standard_seq_real_indexes_to_ali_index = data.get_standard_real_indexes_to_ali_index()

        standard_lrrs = seq_ids_to_lrrs[standard_seq_id]
        standard_lrrs_offset_indexes = set([lrr.offset for lrr in standard_lrrs])
        ali_lrrs_offset_indexes = set([standard_seq_real_indexes_to_ali_index[offset] for offset in standard_lrrs_offset_indexes])

        for seq_id in data.seq_ids_to_seq:
            lrrs = seq_ids_to_lrrs[seq_id]
            real_offset_indexes = set([lrr.offset for lrr in lrrs])

            match_count = 0
            ali_indexes_to_real_index = data.get_ali_indexes_to_real_index(seq_id)
            for standard_lrr_ali_index in ali_lrrs_offset_indexes:
                real_index = ali_indexes_to_real_index.get(standard_lrr_ali_index, None)
                if real_index is None:
                    continue
                if real_index in real_offset_indexes:
                    match_count += 1
            logging.info(str.format("Seq {}, standard seq {}: lrr count {}, ali count/match count {}/{}, match percent {}%",
                                    seq_id, data.standard_seq_id, len(lrrs),
                                    len(ali_lrrs_offset_indexes), match_count, match_count / len(lrrs) * 100))


def step3_10(seq_ids_to_lrrs, seqs):
    seq_ids_to_file_path = {
        'AT1G55610.2': ("mafft/brl1_homo_10_ali.fasta", "iqtree/BRL1_figtree.tree"),
        'AT4G39400.1': ("mafft/bri1_homo_100_align.fasta", "iqtree/bri1_figtree.tree"),
        'AT5G46330.1': ("mafft/fls2_homo_10_ali.fasta", "iqtree/FLS2_figtree.tree"),
        'AT4G28490.1': ("mafft/hae_top10_id_ali.fasta", "iqtree/HAE_figtree.tree"),
        'AT1G73080.1': ("mafft/pepr1_top10_ali.fasta", "iqtree/pepr1_figtree.tree"),
        'AT2G02220.1': ("mafft/pskr1_top10_ali.fasta", "iqtree/PSKR1_figtree.tree"),
        'AT5G61480.1': ("mafft/PXY_top10_ali.fasta", "iqtree/PXY_figtree.tree"),
        "AT4G26540.1": ("mafft/RGFR1_top10_ali.fasta", "iqtree/RGFR1_figtree.tree")
    }

    seq_ids_to_data = read_in_ali_seqs(seq_ids_to_file_path)
    relevant_seq_ids = set([seq_id for data in seq_ids_to_data.values()
                                   for seq_id in data.seq_ids_to_seq.keys()])

    lost_seq_ids = relevant_seq_ids - set(seq_ids_to_lrrs.keys())
    lost_seq_ids_to_lrrs = dao.find_motifs_by_seq_ids(lost_seq_ids, MOTIF_VERSION, with_wrong=False)
    seq_ids_to_lrrs.update(lost_seq_ids_to_lrrs)

    seq_ids_to_seq = dict([(seq.seq_id, seq) for seq in seqs if seq.seq_id in relevant_seq_ids])

    generate_step3_10_weblogo_per_offset(seq_ids_to_data, seq_ids_to_lrrs)
    generate_step3_10_weblog_lrrs(seq_ids_to_data, seq_ids_to_seq, seq_ids_to_lrrs)
    generate_step3_10_weblogo_by_types(seq_ids_to_data, seq_ids_to_seq, seq_ids_to_lrrs)
    statitics_match_count(seq_ids_to_data, seq_ids_to_lrrs)


def _get_subgroup(seq, seq_ids_to_lrrs):
    subgroup = seq.subgroup
    if subgroup == '10':
        if len(seq_ids_to_lrrs[seq.seq_id]) <= 20:
            return '10-1'
        else:
            return '10-2'
    return subgroup


def step3_11(seqs, seq_ids_to_lrrs):
    seq_ids = set([seq.seq_id for seq in seqs])
    seq_ids_to_nsites = dao.find_nsites_by_seq_ids(seq_ids)
    subgroups_to_ntypes_to_count = {}
    for seq in seqs:
        subgroup = _get_subgroup(seq, seq_ids_to_lrrs)
        if subgroup is None:
            continue
        if subgroup not in subgroups_to_ntypes_to_count:
            subgroups_to_ntypes_to_count[subgroup] = {'S': 0, 'T': 0}
        nsites = seq_ids_to_nsites.get(seq.seq_id, [])
        for nsite in nsites:
            if nsite.ntype == dao.nsite.S:
                subgroups_to_ntypes_to_count[subgroup]['S'] += 1
            else:
                subgroups_to_ntypes_to_count[subgroup]['T'] += 1

    generate_histogram(subgroups_to_ntypes_to_count, 'ntype', 'step3_11_n_count_by_subgroup_type.png')


def supplement_for_review(seqs, seq_ids_to_nsites, seq_ids_to_lrrs):
    data, none_n_data, none_n_without9_data = generate_step3_subgroup_type_data(seqs, seq_ids_to_nsites, seq_ids_to_lrrs,
                                                          list('CHE'), lambda seq, index: seq.ss[index], 0)
    logging.debug(str.format("Generate step3.5 web logo data end {}", data.subgroups_to_types_to_seqs_around_nsites))
    logging.debug(str.format("Generate step3.5 lrr acc20 data end {}", data.subgroups_to_lrrs_acc20mean))
    generate_step3_5_lrr_acc20_line_chart(data.subgroups_to_lrrs_acc20mean, 'supplement_')


def main():
    with dao.query_session() as session:
        seqs = dao.sequence.find_all_seqs(session)

    step1(seqs)
    logging.info("Step 1 end")

    # 3.4开始，只关注特定的5个亚家族，因此对seq做精简
    relevant_seqs = [seq for seq in seqs if seq.subgroup in SUBGROUPS]
    relevant_seq_ids = [seq.seq_id for seq in relevant_seqs]
    seq_ids_to_nsites = dao.find_nsites_by_seq_ids(relevant_seq_ids)
    seq_ids_to_lrrs = dao.find_motifs_by_seq_ids(relevant_seq_ids, MOTIF_VERSION, with_wrong=False)

    step3_5(relevant_seqs, seq_ids_to_nsites, seq_ids_to_lrrs)
    logging.info("Step 3.5 end")

    step3_7(relevant_seqs, seq_ids_to_nsites, seq_ids_to_lrrs)
    logging.info("Step 3.7 end")

    step3_9(relevant_seqs, seq_ids_to_nsites, seq_ids_to_lrrs)
    logging.info("Step 3.9 end")

    step3_10(seq_ids_to_lrrs, seqs)
    logging.info("Step 3.10 end")

    step3_11(seqs, seq_ids_to_lrrs)
    logging.info("Step 3.11 end")

    full_seq_ids = [seq.seq_id for seq in seqs]
    full_seq_ids_to_nsites = dao.find_nsites_by_seq_ids(full_seq_ids)
    full_seq_ids_to_lrrs = dao.find_motifs_by_seq_ids(full_seq_ids, MOTIF_VERSION, with_wrong=False)
    supplement_for_review(seqs, full_seq_ids_to_nsites, full_seq_ids_to_lrrs)
    logging.info("Step supplement for review end")


if __name__ == "__main__":
    logging.basicConfig(level=logging.DEBUG)
    import os
    if not os.path.exists(OUTPUT_PATH):
        os.mkdir(OUTPUT_PATH)

    root = logging.getLogger()
    fh = logging.FileHandler(os.path.join(OUTPUT_PATH, "log.log"))
    fh.setLevel(logging.DEBUG)
    fh.setFormatter(logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s'))
    root.addHandler(fh)
    main()
