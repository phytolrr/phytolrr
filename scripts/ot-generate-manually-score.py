import logging
logging.basicConfig(level=logging.INFO)
import dao
from tools import motifs as motif_tool, pssm_matrix

MOTIF_VERSION = 2


with dao.session_scope() as session:
    matrix = pssm_matrix.calc_pssm_matrix(dao.find_baseline_motifs())

    logging.info(matrix.pssm)

    cls = dao.motif.get_entity(MOTIF_VERSION)
    motifs = session.query(cls).filter(cls.manually_add).all()

    results = []
    for m in motifs:
        seq = dao.find_seq_by_id(m.seq_id)
        motif_seq = seq.seq[m.offset:m.offset+24]
        score = motif_tool.calc_pssm_score(motif_seq[:16], matrix)
        m.score = score
        results.append({'seq_id': m.seq_id, 'offset': m.offset, 'seq': motif_seq, 'score': score})

    logging.info(str.format("Total manually records {}", len(results)))
    results.sort(key=lambda k: k['score'])
    for r in results:
        logging.info(str.format("Sequence {seq_id}, offset {offset}, seq {seq}, score {score}", **r))


