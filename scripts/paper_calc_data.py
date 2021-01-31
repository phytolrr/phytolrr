import logging
import argparse
import csv
from collections import OrderedDict
import numpy as np
from tools import files
import settings
import dao


MOTIF_VERSION = 3


def read_in_config():
    parser = argparse.ArgumentParser()
    parser.add_argument('seqid', nargs="*", help='The sequence ids')
    parser.add_argument('--seqid-file', '-f', help='The path of the file which contains sequence ids')
    parser.add_argument('--output-file', '-o', help='The path of the output file which in csv format')
    args = parser.parse_args()

    seq_ids = []
    if args.seqid_file is not None:
        seq_ids += files.get_names_from_file(args.seqid_file)
    seq_ids += args.seqid
    return seq_ids, args.output_file


def query_density(seq_id):
    with dao.query_session() as session:
        motifs = dao.motif.find_motifs_by_seq_id(session, seq_id, MOTIF_VERSION, with_wrong=False)
        nsites = dao.nsite.find_nsites_by_seq_id(session, seq_id)
        seq = dao.find_seq_by_id(seq_id)
    return motifs, nsites, seq


def statistic_nsites_pos(ss, nsites):
    count = {
        'C': 0,
        'H': 0,
        'E': 0
    }
    for nsite in nsites:
        count[ss[nsite.start_pos]] += 1
    return count


def get_by_start_pos(content, start):
    begin = start - 5
    end = start + 6
    if begin < 0:
        begin = 0
    result = content[begin:end]
    if len(result) < 11:
        return None
    return result


def statistic_aminos_around_nsites(seq, nsites, ss):
    ss_to_sub_seqs = {
        'C': [],
        'H': [],
        'E': []
    }
    for nsite in nsites:
        sub_seq = get_by_start_pos(seq, nsite.start_pos)
        if sub_seq is None:
            logging.warning(str.format("Drop start pos {}", nsite.start_pos))
            continue
        ss_to_sub_seqs.get(ss[nsite.start_pos]).append(sub_seq)
    return ss_to_sub_seqs


def statistic_acc20_arround_nsites(acc20, nsites, ss):
    acc20_list = [int(aa_acc20) for aa_acc20 in acc20.split()]
    sses_to_acc20s = {
        'C': [],
        'H': [],
        'E': []
    }
    if len(nsites) == 0:
        return sses_to_acc20s

    for nsite in nsites:
        acc20_around_nsites = get_by_start_pos(acc20_list, nsite.start_pos)
        if acc20_around_nsites is None:
            logging.warning(str.format("Drop start pos {}", nsite.start_pos))
            continue

        acc20s_around_nsites = sses_to_acc20s.get(ss[nsite.start_pos])
        acc20s_around_nsites.append(acc20_around_nsites)

    sses_to_acc20mean = {}
    for s, acc20s_around_nsites in sses_to_acc20s.items():

        if len(acc20s_around_nsites) != 0:
            acc20s_matrix = np.array(acc20s_around_nsites)
            acc20s_mean = np.mean(acc20s_matrix, axis=0).tolist()
            sses_to_acc20mean[s] = acc20s_mean
        else:
            sses_to_acc20mean[s] = []
    return sses_to_acc20mean


def statistic_nsites_located_in_lrr(motifs, nsites, lrr_length=9):
    pos_in_motifs = set()
    for motif in motifs:
        poses = range(motif.offset, motif.offset + lrr_length)
        pos_in_motifs.update(poses)
    count_in_motif = 0
    count_out_motif = 0
    for nsite in nsites:
        if nsite.start_pos in pos_in_motifs:
            count_in_motif += 1
        else:
            count_out_motif += 1
    return count_in_motif, count_out_motif


def write_to_file(file_path, rows):
    if len(rows) == 0:
        logging.warning("No rows to write to file {}", file_path)
        return
    with open(file_path, 'w') as f:
        cf = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        cf.writeheader()
        cf.writerows(rows)


def main():
    seq_ids, output_file = read_in_config()

    rows = []
    for seq_id in seq_ids:
        logging.info(str.format("Begin to statistic seq {}", seq_id))
        motifs, nsites, seq = query_density(seq_id)
        ss_counts = statistic_nsites_pos(seq.ss, nsites)
        ss_to_sub_seqs = statistic_aminos_around_nsites(seq.seq, nsites, seq.ss)
        sses_to_acc20mean = statistic_acc20_arround_nsites(seq.acc20, nsites, seq.ss)
        nsites_count_in_lrr, nsites_count_outer_lrr = statistic_nsites_located_in_lrr(motifs, nsites)
        nsites_count_in_lrr16, nsites_count_outer_lrr16 = statistic_nsites_located_in_lrr(motifs, nsites, 16)
        nsites_count_in_lrr24, nsites_count_outer_lrr24 = statistic_nsites_located_in_lrr(motifs, nsites, 24)
        row = OrderedDict(
            seq_id=seq_id,
            seq_len=len(seq.seq),
            motif_count=len(motifs),
            nsite_count=len(nsites),
            c_count=ss_counts['C'],
            h_count=ss_counts['H'],
            e_count=ss_counts['E'],
            c_nsites_around=';'.join(ss_to_sub_seqs.get('C')),
            h_nsites_around=';'.join(ss_to_sub_seqs.get('H')),
            e_nsites_around=';'.join(ss_to_sub_seqs.get('E')),
            c_avr_acc20_around_nsites=' '.join([str.format('{:.2f}', avr) for avr in sses_to_acc20mean['C']]),
            h_avr_acc20_around_nsites=' '.join([str.format('{:.2f}', avr) for avr in sses_to_acc20mean['H']]),
            e_avr_acc20_around_nsites=' '.join([str.format('{:.2f}', avr) for avr in sses_to_acc20mean['E']]),
            nsites_in_lrr_count=nsites_count_in_lrr,
            nsites_outer_lrr_count=nsites_count_outer_lrr,
            nsites_in_lrr_count16=nsites_count_in_lrr16,
            nsites_outer_lrr_count16=nsites_count_outer_lrr16,
            nsites_in_lrr_count24=nsites_count_in_lrr24,
            nsites_outer_lrr_count24=nsites_count_outer_lrr24)
        rows.append(row)
    if len(rows) == 0:
        logging.warning("No sequences get")
        return
    if output_file is None:
        logging.warning(str.format("Get {} sequences total, but no output file specifed", len(rows)))
        return
    write_to_file(output_file, rows)


if __name__ == '__main__':
    logging.basicConfig(level=logging.WARNING)
    main()
