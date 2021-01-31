# -*- coding: utf-8 -*
"""
根据id list，统计在不同二级结构的Nsites前后5个氨基酸序列的acc20
"""

import logging
import argparse
import numpy as np
from tools import files
import settings
import dao


MOTIF_VERSION = 3


def read_in_config():
    parser = argparse.ArgumentParser()
    parser.add_argument('seqid', nargs="*", help='The sequence ids')
    parser.add_argument('--seqid-file', '-f', help='The path of the file which contains sequence ids')
    args = parser.parse_args()

    seq_ids = []
    if args.seqid_file is not None:
        seq_ids += files.get_names_from_file(args.seqid_file)
    seq_ids += args.seqid
    return seq_ids


def obtain_value(seq_id):
    with dao.query_session() as session:
        motifs = dao.motif.find_motifs_by_seq_id(session, seq_id, MOTIF_VERSION, with_wrong=False)
        nsites = dao.nsite.find_nsites_by_seq_id(session, seq_id)
        seq = dao.find_seq_by_id(seq_id)

    return motifs, nsites, seq


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
        acc20s_around_nsites = sses_to_acc20s.get(ss[nsite.start_pos])
        start = nsite.start_pos - 5
        end = nsite.start_pos + 6
        if start < 0:
            start = 0
        acc20_around_nsites = acc20_list[start:end]
        if len(acc20_around_nsites) < 11:
            continue
        acc20s_around_nsites.append(acc20_around_nsites)
    return sses_to_acc20s


def main():
    seq_ids = read_in_config()
    sses_to_acc20 = {}
    sses_to_acc20_one_time = {}
    for seq_id in seq_ids:
        motifs, nsites, seq = obtain_value(seq_id)
        sses_to_acc20_one_time = statistic_acc20_arround_nsites(seq.acc20, nsites, seq.ss)
        for key, value in sses_to_acc20_one_time.items():
            sses_to_acc20[key] = sses_to_acc20.get(key, []) + value
    sses_to_acc20mean = {}
    for s, acc20s_around_nsites in sses_to_acc20.items():
        if len(acc20s_around_nsites) != 0:
            acc20s_matrix = np.array(acc20s_around_nsites)
            acc20s_mean = np.mean(acc20s_matrix, axis=0).tolist()
            sses_to_acc20mean[s] = acc20s_mean
        else:
            sses_to_acc20mean[s] = []
    for key in sses_to_acc20mean:
        print(key + '\n' + str(sses_to_acc20mean[key]))




if __name__ == '__main__':
    logging.basicConfig(level=logging.INFO)
    main()