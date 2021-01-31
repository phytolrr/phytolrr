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
        seq = dao.find_seq_by_id(seq_id)

    return motifs, seq


def extract_lrrs(seq_id, motifs, seq):
    lrrs_for_ids = {
        seq_id: [],
    }
    for motif in motifs:
        start = motif.offset
        end = motif.offset + 24

        lrrs_for_ids.get(seq_id).append(seq[start:end])
    return lrrs_for_ids


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
        motifs, seq = query_density(seq_id)
        lrrs_from_seqs = extract_lrrs(seq_id, motifs, seq.seq)

        row = OrderedDict(
            seq_id=seq_id,
            LRRs=';'.join(lrrs_from_seqs.get(seq_id)))
        rows.append(row)
    if len(rows) == 0:
        logging.warning("No sequences get")
        return
    if output_file is None:
        logging.warning(str.format("Get {} sequences total, but no output file specifed", len(rows)))
        return
    write_to_file(output_file, rows)


if __name__ == '__main__':
    logging.basicConfig(level=logging.INFO)
    main()
