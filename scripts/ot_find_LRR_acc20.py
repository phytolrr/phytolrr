# -*- coding: utf-8 -*
"""
根据id list，统计LRRd的acc20
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
        seq = dao.find_seq_by_id(seq_id)

    return motifs, seq


def nsites_located_in_lrrs(acc20, motifs):
    acc20_list = [int(aa_acc20) for aa_acc20 in acc20.split()]
    acc20_to_lrr = []

    if len(motifs) == 0:
        return acc20_to_lrr

    for motif in motifs:
        start = motif.offset
        end = motif.offset + 24
        if len(acc20_list[start:end]) < 24:
            continue
        acc20_to_lrr.append(acc20_list[start:end])
    return acc20_to_lrr


def main():
    seq_ids = read_in_config()
    lrr_acc20 = []
    lrr_acc20_one_time = []
    logging.debug(str.format("All seqs {}", seq_ids))
    for seq_id in seq_ids:
        motifs, seq = obtain_value(seq_id)
        logging.debug(str.format("Motifs for seq {}， count {}", seq.seq_id, len(motifs)))
        lrr_acc20_one_time = nsites_located_in_lrrs(seq.acc20, motifs)
        lrr_acc20 += lrr_acc20_one_time
    acc20s_mean = []
    logging.info(str.format("Total motifs count {}, sequences count {}", len(lrr_acc20), len(seq_ids)))
    if len(lrr_acc20) != 0:
        acc20s_matrix = np.array(lrr_acc20)
        acc20s_mean = np.mean(acc20s_matrix, axis=0).tolist()
    else:
        acc20s_mean = []
    print(str.format("Finished, length {}, content {}", len(acc20s_mean), acc20s_mean))


if __name__ == '__main__':
    logging.basicConfig(level=logging.INFO)
    main()