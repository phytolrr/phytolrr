# -*- coding: utf-8 -*
'''
通过-f输入id list，将序列的lrr找出来用于序列比对，以期将LRR分类,通过print，直接输出fasta格式，
以序列id和数据库中的id作为name
'''

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
        seq = dao.find_seq_by_id(seq_id)

    return motifs, seq


def extract_lrrs(seq_id, motifs, seq):

    lrrs_for_ids = {}
    for motif in motifs:
        id_key = seq_id + str(motif.id)
        start = motif.offset
        end = motif.offset + 24

        lrrs_for_ids[id_key] = seq[start:end]
    return lrrs_for_ids


def main():
    seq_ids = read_in_config()

    for seq_id in seq_ids:
        motifs, seq = obtain_value(seq_id)
        lrrs_from_seqs = extract_lrrs(seq_id, motifs, seq.seq)
        for key in lrrs_from_seqs:
            print ('>'+ key + '\n' + lrrs_from_seqs [key])


if __name__ == '__main__':
    logging.basicConfig(level=logging.INFO)
    main()





