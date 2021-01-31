# -*- coding: utf-8 -*
"""
根据id list，通过二级结构，分别统计Nsites前后5个氨基酸序列，与非糖基化位点的N前后5个氨基酸序列，做weblogo
去除不等于11的序列
代码相对于compare_N_and_Nmotif.py 重构优化
"""

import logging
import argparse
from tools import files
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
        nsites = dao.nsite.find_nsites_by_seq_id(session, seq_id)
        seq = dao.find_seq_by_id(seq_id)

    return nsites, seq


def n_position_in_ss(seq_id, seq, nsites, ss):
    pos_in_nistes = set()
    ss_to_sub_seqs = {
        'C': {},
        'H': {},
        'E': {}
    }
    for nsite in nsites:
        pos_in_nistes.add(nsite.start_pos)
    for index in range(len(seq)):
        if seq[index] == "N" and index not in pos_in_nistes:
            ss_to_sub_seq = ss_to_sub_seqs.get(ss[index])

            id_key = ''
            if ss[index] == "C":
                id_key = seq_id + "nCoil" + str(index)
            if ss[index] == "H":
                id_key = seq_id + "nHelix" + str(index)
            if ss[index] == "E":
                id_key = seq_id + "nEsheet" + str(index)

            start = index - 5
            end = index + 6
            if start >= 0:
                if len(seq[start:end]) == 11:
                    ss_to_sub_seq[id_key] = seq[start:end]
        else:
            continue
    return ss_to_sub_seqs


def main():
    seq_ids = read_in_config()
    n_in_ss_onetime = {}
    n_in_ss = {
        'C': {},
        'H': {},
        'E': {}
    }

    for seq_id in seq_ids:
        nsites, seq = obtain_value(seq_id)
        n_in_ss_onetime = n_position_in_ss(seq_id, seq.seq, nsites, seq.ss)
        for key in n_in_ss_onetime:
            n_in_ss.get(key).update(n_in_ss_onetime.get(key))

    for key, value in n_in_ss.items():
        if key == 'C':
            for my_id, my_seqs in value.items():
                print('>' + my_id + '\n' + my_seqs)
        if key == 'E':
            for my_id, my_seqs in value.items():
                print('>' + my_id + '\n' + my_seqs)
        if key == 'H':
            for my_id, my_seqs in value.items():
                print('>' + my_id + '\n' + my_seqs)


if __name__ == '__main__':
    logging.basicConfig(level=logging.INFO)
    main()