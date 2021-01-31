import logging
import csv
from collections import OrderedDict
import dao
from sqlalchemy import func
from Bio.Seq import Seq
from Bio import motifs
from Bio.Alphabet import IUPAC
import os
import weblogolib
from corebio.seq_io import array_io
import io


SEQ_ID = 'AT%'
LRR_COUNT = 20
MOTIF_VERSION = 3

AMINO='S'
POS = 6

OUTPUT_FILE = "../test_files/s_probability.csv"
WEBLOGO_PATH = "../test_files/"


def find_seqs():
    motif_cls = dao.motif.get_entity(MOTIF_VERSION)
    with dao.query_session() as session:
        return session.query(dao.sequence.SequenceEntity)\
                            .join(motif_cls, dao.sequence.SequenceEntity.seq_id == motif_cls.seq_id)\
                            .filter(dao.sequence.SequenceEntity.seq_id.like(SEQ_ID))\
                            .filter(motif_cls.correct)\
                            .group_by(dao.sequence.SequenceEntity.seq_id)\
                            .having(func.count(motif_cls.id) > LRR_COUNT).all()


def generate_weblogo(seq_id, seq_strs):
    seqs = weblogolib.read_seq_data(io.StringIO('\n'.join(seq_strs)), input_parser=array_io.read)
    data = weblogolib.LogoData.from_seqs(seqs)
    options = weblogolib.LogoOptions()
    options.title = seq_id
    options.unit_name = "probability"
    format = weblogolib.LogoFormat(data, options)
    eps = weblogolib.eps_formatter(data, format)
    with open(os.path.join(WEBLOGO_PATH, seq_id + '.eps'), 'wb') as f :
        f.write(eps)


def main():
    seqs = find_seqs()
    seq_ids_to_seq = dict([(seq.seq_id, seq) for seq in seqs])
    logging.info("Seqs count {}, seq_ids {}", len(seqs), seq_ids_to_seq.keys())
    seq_ids_to_motifs = dao.find_motifs_by_seq_ids(seq_ids_to_seq.keys(), MOTIF_VERSION, with_wrong=False)

    seq_ids_to_matrix = {}
    for seq_id, motif_entities in seq_ids_to_motifs.items():
        motif_entities.sort(key=lambda m:m.offset)
        seq_str = seq_ids_to_seq[seq_id].seq
        motif_seqs_str = [seq_str[m.offset:m.offset+16] for m in motif_entities]
        generate_weblogo(seq_id, motif_seqs_str)

        motif_seq = [Seq(motif_seq_str, IUPAC.protein) for motif_seq_str in motif_seqs_str]
        matrix = motifs.create(motif_seq, IUPAC.protein)

        seq_ids_to_matrix[seq_id] = matrix

    seq_ids = list(seq_ids_to_matrix.keys())
    seq_ids.sort()

    rows = []
    for seq_id in seq_ids:
        matrix = seq_ids_to_matrix.get(seq_id)
        print(matrix.pwm[AMINO][POS])
        print(str.format("seq {}:\n{}",seq_id, matrix))
        row = OrderedDict(
            seq_id=seq_id,
            s_num=matrix.counts[AMINO][POS],
            probability=matrix.pwm[AMINO][POS],
            LRRs='\n'.join([seq_ids_to_seq[seq_id].seq[m.offset:m.offset+16] for m in seq_ids_to_motifs[seq_id]])
        )
        rows.append(row)

    with open(OUTPUT_FILE, 'w') as f:
        c = csv.DictWriter(f, fieldnames=['seq_id', 's_num', 'probability', 'LRRs'])
        c.writeheader()
        c.writerows(rows)


if __name__ == "__main__":
    main()
