# -*- coding: utf-8 -*
"""
根据id list，通过二级结构，分别统计Nsites前后5个氨基酸序列，与非糖基化位点的N前后5个氨基酸序列，做weblogo
去除不等于11的序列
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


def n_position_in_coil(seq_id, seq, nsites, ss):
    pos_in_nistes = set()
    ss_to_sub_seqs = {}
    for nsite in nsites:
        pos_in_nistes.add(nsite.start_pos)
    for index in range(len(seq)):
        if seq[index] == "N" and index not in pos_in_nistes:
            if ss[index] == "C":
                id_key = seq_id + "nCoil" + str(index)
                start = index - 5
                end = index + 6
                if start >= 0:
                    if len(seq[start:end]) == 11:
                        ss_to_sub_seqs[id_key] = seq[start:end]
                else:
                    continue
        else:
            continue
    return ss_to_sub_seqs

def n_position_in_strand(seq_id, seq, nsites, ss):
    pos_in_nistes = set()
    ss_to_sub_seqs = {}
    for nsite in nsites:
        pos_in_nistes.add(nsite.start_pos)
    for index in range(len(seq)):
        if seq[index] == "N" and index not in pos_in_nistes:
            if ss[index] == "E":
                id_key = seq_id + "nEsheet" + str(index)
                start = index - 5
                end = index + 6
                if start >= 0:
                    if len(seq[start:end]) == 11:
                        ss_to_sub_seqs[id_key] = seq[start:end]
                else:
                    continue
        else:
            continue
    return ss_to_sub_seqs

def n_position_in_helix(seq_id, seq, nsites, ss):
    pos_in_nistes = set()
    ss_to_sub_seqs = {}
    for nsite in nsites:
        pos_in_nistes.add(nsite.start_pos)
    for index in range(len(seq)):
        if seq[index] == "N" and index not in pos_in_nistes:
            if ss[index] == "H":
                id_key = seq_id + "nHelix" + str(index)
                start = index - 5
                end = index + 6
                if start >= 0:
                    if len(seq[start:end]) == 11:
                        ss_to_sub_seqs[id_key] = seq[start:end]
                else:
                    continue
        else:
            continue
    return ss_to_sub_seqs


def main():
    seq_ids = read_in_config()

    for seq_id in seq_ids:
        nsites, seq = obtain_value(seq_id)
        n_in_coil = n_position_in_coil(seq_id, seq.seq, nsites, seq.ss)
        for key in n_in_coil:
            print('>' + key + '\n' + n_in_coil[key])
    for seq_id in seq_ids:
        nsites, seq = obtain_value(seq_id)
        n_in_strand = n_position_in_strand(seq_id, seq.seq, nsites, seq.ss)
        for key in n_in_strand:
            print('>' + key + '\n' + n_in_strand[key])
    for seq_id in seq_ids:
        nsites, seq = obtain_value(seq_id)
        n_in_helix = n_position_in_helix(seq_id, seq.seq, nsites, seq.ss)
        for key in n_in_helix:
            print('>' + key + '\n' + n_in_helix[key])


if __name__ == '__main__':
    logging.basicConfig(level=logging.INFO)
    main()