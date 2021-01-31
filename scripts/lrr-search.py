import logging
import argparse
import os
from multiprocessing import Pool
import dao
from tools import motifs as motif_tool, pssm_matrix

MOTIF_VERSION = 3
BASELINE_MOTIF_VERSION = 1
VALID_AMINO = {'A', 'R', 'N', 'D', 'C', 'E', 'Q', 'G', 'H', 'V', 'I', 'L', 'K', 'M', 'F', 'P', 'S', 'T', 'W', 'Y'}


class Config(object):
    def __init__(self):
        parser = argparse.ArgumentParser()
        parser.add_argument('-f', '--sequences-file', required=True, help='The sequences file')
        parser.add_argument('-t', '--training-lrr-file', required=True)
        self._parser = parser

    def parse(self):
        cfg = self._parser.parse_args()
        if not os.path.isfile(cfg.sequences_file):
            print(str.format("The sequence file {} does not exists", cfg.sequences_file))
            exit(1)
        if not os.path.isfile(cfg.training_lrr_file):
            print(str.format("The file {} does not exists", cfg.training_lrr_file))
            exit(1)
        return cfg


class CalculateTask(object):
    def __init__(self, matrix, seq):
        self.matrix = matrix
        self.seq = seq


def _print_debug(motifs, seq):
    if logging.root.level > logging.DEBUG:
        return
    offsets = [m.offset for m in motifs]
    offsets.sort()

    seq_str = str(seq.seq)
    for motif in motifs:
        offset = motif.offset
        logging.debug(str.format("Offset {}, score {} LRR {}", offset, motif.score, seq_str[offset:offset+16]))


def find_lrr_and_save_db(matrix, seq):
    seq_ids_to_motifs = dao.find_motifs_by_seq_ids([seq.seq_id,], MOTIF_VERSION)
    if len(seq_ids_to_motifs.get(seq.seq_id, [])) > 0:
        logging.warning(str.format("The sequence {} is already analysied before, skip it", seq.seq_id))
        return

    seq_str = str(seq.seq)

    invalid_amino = set(seq_str) - VALID_AMINO
    if len(invalid_amino) > 0:
        logging.error(str.format("The sequence {} contains invalid amino {}", seq.seq_id, invalid_amino))
        return

    logging.info(str.format("Begin to find lrr in seq {}, {}", seq.seq_id, seq_str))
    motifs = motif_tool.lrr_search(matrix, str(seq.seq))
    motifs = motif_tool.found_no_overlapped_motifs(motifs)
    for m in motifs:
        m.seq_id = seq.seq_id
        m.correct = True

    logging.info(str.format("Found {} LRR motifs in sequence {}", len(motifs), seq.seq_id))
    _print_debug(motifs, seq)

    # Write to db
    motif_entities = [dao.motif_entity.MotifEntityBase(**m.__dict__) for m in motifs]
    logging.debug("Begin to write to db")
    with dao.session_scope() as session:
        dao.motif.replace_motifs_by_seq(session, seq.seq_id, motif_entities, MOTIF_VERSION)


def find_lrr_and_save_db_task(task):
    logging.basicConfig(level=logging.DEBUG,
                        format='%(asctime)s:Process %(process)d:%(levelname)s:%(message)s')
    dao.reconnect()
    logging.debug(str.format("Begin to run task for seq {}", task.seq.seq_id))
    find_lrr_and_save_db(task.matrix, task.seq)


def main():
    # generate the matrix
    seqs = dao.find_baseline_seqs()
    seq_ids_to_seq_str = dict([(seq.seq_id, seq.seq) for seq in seqs])
    logging.info(str.format("Baseline sequence ids count({}): {}", len(seq_ids_to_seq_str), seq_ids_to_seq_str.keys()))

    seq_ids_to_motifs = dao.find_motifs_by_seq_ids(seq_ids_to_seq_str.keys(), BASELINE_MOTIF_VERSION, with_wrong=False)
    motif_strs = [seq_ids_to_seq_str[m.seq_id][m.offset:m.offset+16] for motifs in seq_ids_to_motifs.values() for m in motifs]
    logging.info(str.format("Baseline LRR motifs( count {}): {}", len(motif_strs), motif_strs))

    matrix = pssm_matrix.calc_pssm_matrix(motif_strs)
    logging.info(str.format("Matrix: {}", matrix))
    logging.info(str.format("PSSM: {}", matrix.pssm))

    # find the lrr
    with dao.query_session() as session:
        all_seqs = dao.sequence.find_all_seqs(session)

    tasks = [CalculateTask(matrix, seq) for seq in all_seqs]

    with Pool(12) as pool:
        pool.map(find_lrr_and_save_db_task, tasks)
    logging.info("All tasks done")


if __name__ == '__main__':
    logging.basicConfig(level=logging.DEBUG,
                        format='%(asctime)s:Process %(process)d:%(levelname)s:%(message)s')
    main()
