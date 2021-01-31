# -*- coding: utf-8 -*
"""
根据id list，统计在LRR超螺旋内部和外部的Nsites前后5个氨基酸序列，做weblogo
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


def nsites_located_inside_lrr(seq_id, motifs, nsites, seq):
    pos_in_motifs = set()
    ss_to_sub_seqs = {}
    for motif in motifs:
        poses = range(motif.offset, motif.offset + 9)
        pos_in_motifs.add(poses)
    for nsite in nsites:
        if nsite.start_pos in pos_in_motifs:
            id_key = seq_id + "inside" + str(nsite.start_pos)
            start = nsite.start_pos - 5
            end = nsite.start_pos + 6
            if start >= 0:
                if len(seq[start:end]) == 11:
                    ss_to_sub_seqs[id_key] = seq[start:end]
            else:
                continue

    return ss_to_sub_seqs


def nsites_located_beside_lrr(seq_id, motifs, nsites, seq):
    pos_in_motifs = set()
    ss_to_sub_seqs = {}
    for motif in motifs:
        poses = range(motif.offset + 9, motif.offset + 13)
        pos_in_motifs.update(poses)
    for nsite in nsites:
        if nsite.start_pos in pos_in_motifs:
            id_key = seq_id + "beside" + str(nsite.start_pos)
            start = nsite.start_pos - 5
            end = nsite.start_pos + 6
            if start >= 0:
                if len(seq[start:end]) == 11:
                    ss_to_sub_seqs[id_key] = seq[start:end]
            else:
                continue

    return ss_to_sub_seqs


def nsites_located_back_lrr(seq_id, motifs, nsites, seq):
    pos_in_motifs = set()
    ss_to_sub_seqs = {}
    for motif in motifs:
        poses = range(motif.offset + 16, motif.offset + 24)
        pos_in_motifs.update(poses)
    for nsite in nsites:
        if nsite.start_pos in pos_in_motifs:
            id_key = seq_id + "back" + str(nsite.start_pos)
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
        nsites_inside_lrr = nsites_located_inside_lrr(seq_id, motifs, nsites, seq.seq)
        for key in nsites_inside_lrr:
            print('>' + key + '\n' + nsites_inside_lrr[key])
    for seq_id in seq_ids:
        motifs, nsites, seq = obtain_value(seq_id)
        nsites_beside_lrr = nsites_located_beside_lrr(seq_id, motifs, nsites, seq.seq)
        for key in nsites_beside_lrr:
            print('>' + key + '\n' + nsites_beside_lrr[key])
    for seq_id in seq_ids:
        motifs, nsites, seq = obtain_value(seq_id)
        nsites_back_lrr = nsites_located_back_lrr(seq_id, motifs, nsites, seq.seq)
        for key in nsites_back_lrr:
            print('>' + key + '\n' + nsites_back_lrr[key])


if __name__ == '__main__':
    logging.basicConfig(level=logging.INFO)
    main()