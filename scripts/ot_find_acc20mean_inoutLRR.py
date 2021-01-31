# -*- coding: utf-8 -*
"""
根据id list，统计在LRR超螺旋内外面的Nsites前后5个氨基酸序列的acc20
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


def nsites_located_in_lrrs(acc20, motifs, nsites):
    acc20_list = [int(aa_acc20) for aa_acc20 in acc20.split()]
    pos_in_motifs_inside = set()
    pos_in_motifs_beside = set()
    pos_in_motifs_back = set()
    acc20_to_inout_lrr = {
        'inside': [],
        'beside': [],
        'back': []
    }

    if len(nsites) == 0:
        return acc20_to_inout_lrr

    for motif in motifs:
        poses_inside = range(motif.offset, motif.offset + 9)
        pos_in_motifs_inside.update(poses_inside)
        poses_beside = range(motif.offset + 9, motif.offset + 16)
        pos_in_motifs_beside.update(poses_beside)
        poses_back = range(motif.offset + 16, motif.offset + 24)
        pos_in_motifs_back.update(poses_back)
    for nsite in nsites:
        if nsite.start_pos in pos_in_motifs_inside:
            id_key = "inside"
            acc20s_around_nsites_inside = acc20_to_inout_lrr.get(id_key)
            start = nsite.start_pos - 5
            end = nsite.start_pos + 6
            if start < 0:
                start = 0
            if len(acc20_list[start:end]) < 11:
                continue
            acc20s_around_nsites_inside.append(acc20_list[start:end])

        if nsite.start_pos in pos_in_motifs_beside:
            id_key = "beside"
            acc20s_around_nsites_beside = acc20_to_inout_lrr.get(id_key)
            start = nsite.start_pos - 5
            end = nsite.start_pos + 6
            if start < 0:
                start = 0
            if len(acc20_list[start:end]) < 11:
                continue
            acc20s_around_nsites_beside.append(acc20_list[start:end])

        if nsite.start_pos in pos_in_motifs_back:
            id_key = "back"
            acc20s_around_nsites_back = acc20_to_inout_lrr.get(id_key)
            start = nsite.start_pos - 5
            end = nsite.start_pos + 6
            if start < 0:
                start = 0
            if len(acc20_list[start:end]) < 11:
                continue
            acc20s_around_nsites_back.append(acc20_list[start:end])

        else:
            continue

    return acc20_to_inout_lrr


def main():
    seq_ids = read_in_config()
    sses_to_acc20 = {}
    sses_to_acc20_one_time = {}
    for seq_id in seq_ids:
        motifs, nsites, seq = obtain_value(seq_id)
        sses_to_acc20_one_time = nsites_located_in_lrrs(seq.acc20, motifs, nsites)
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