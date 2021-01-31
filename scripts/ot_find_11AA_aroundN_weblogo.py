# -*- coding: utf-8 -*
"""
根据id list，通过二级结构，分别统计Nsites前后5个氨基酸序列，做weblogo
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
        motifs = dao.motif.find_motifs_by_seq_id(session, seq_id, MOTIF_VERSION, with_wrong=False)
        nsites = dao.nsite.find_nsites_by_seq_id(session, seq_id)
        seq = dao.find_seq_by_id(seq_id)

    return motifs, nsites, seq


def n_position_in_coil(seq_id, seq, nsites, ss):
    ss_to_sub_seqs = {}
    for nsite in nsites:
        if ss[nsite.start_pos] == "C":
            id_key = seq_id + "Coil" + str(nsite.start_pos)
            start = nsite.start_pos - 5
            end = nsite.start_pos + 6
            if start >= 0:
                if len(seq[start:end]) == 11:
                    ss_to_sub_seqs[id_key] = seq[start:end]
            else:
                continue
    return ss_to_sub_seqs


def n_position_in_strand(seq_id, seq, nsites, ss):
    ss_to_sub_seqs = {}
    for nsite in nsites:
        if ss[nsite.start_pos] == "E":
            id_key = seq_id + "Esheet" + str(nsite.start_pos)
            start = nsite.start_pos - 5
            end = nsite.start_pos + 6
            if start >= 0:
                if len(seq[start:end]) == 11:
                    ss_to_sub_seqs[id_key] = seq[start:end]
            else:
                continue
    return ss_to_sub_seqs


def n_position_in_helix(seq_id, seq, nsites, ss):
    ss_to_sub_seqs = {}
    for nsite in nsites:
        if ss[nsite.start_pos] == "H":
            id_key = seq_id + "Helix" + str(nsite.start_pos)
            start = nsite.start_pos - 5
            end = nsite.start_pos + 6
            if start >= 0:
                if len(seq[start:end]) == 11:
                    ss_to_sub_seqs[id_key] = seq[start:end]
            else:
                continue
    return ss_to_sub_seqs


def main():
    seq_ids = read_in_config()

    for seq_id in seq_ids:
        motifs, nsites, seq = obtain_value(seq_id)
        n_in_coil = n_position_in_coil(seq_id, seq.seq, nsites, seq.ss)
        for key in n_in_coil:
            print('>' + key + '\n' + n_in_coil[key])
    for seq_id in seq_ids:
        motifs, nsites, seq = obtain_value(seq_id)
        n_in_strand = n_position_in_strand(seq_id, seq.seq, nsites, seq.ss)
        for key in n_in_strand:
            print('>' + key + '\n' + n_in_strand[key])
    for seq_id in seq_ids:
        motifs, nsites, seq = obtain_value(seq_id)
        n_in_helix = n_position_in_helix(seq_id, seq.seq, nsites, seq.ss)
        for key in n_in_helix:
            print('>' + key + '\n' + n_in_helix[key])


if __name__ == '__main__':
    logging.basicConfig(level=logging.INFO)
    main()