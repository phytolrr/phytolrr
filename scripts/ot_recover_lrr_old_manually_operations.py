# -*- coding:utf-8
'''恢复老版本motifentity表的手工操作到新版本表中
涉及到几部内容的恢复：
1. 老版本手工标记的LRR，且新版本不存在的
2. 老版本呢手工标记的错误的、且新版本也找出来的LRR，需要标记为错误的LRR
3. 老版本中存在的、且认为是正确的、且新版本中不存在的，需要恢复为手工标记字段
恢复手工标记的LRR时，仅恢复offset，score由新版本的matrix重新计算得出

顺序为，
1. 标记错误的LRR
2. 恢复老版本存在的LRR
3. 恢复老版本手工标记的LRR
做2、3时，同时检查是否出现了overlap，如果出现overlap，记录错误人工检查
'''
import logging
from sqlalchemy import and_, not_
import dao
from tools import motifs as motif_tools, pssm_matrix


OLD_VERSION = 2
NEW_VERSION = 3


def get_old_wrong_new_exists():
    old_cls = dao.motif.get_entity(OLD_VERSION)
    new_cls = dao.motif.get_entity(NEW_VERSION)
    with dao.query_session() as session:
        return session.query(old_cls, new_cls)\
            .join(new_cls, and_(old_cls.seq_id == new_cls.seq_id, old_cls.offset == new_cls.offset))\
            .filter(not_(old_cls.correct)).all()


def get_old_correct_new_not_exists():
    old_cls = dao.motif.get_entity(OLD_VERSION)
    new_cls = dao.motif.get_entity(NEW_VERSION)
    with dao.query_session() as session:
        return session.query(old_cls)\
                .outerjoin(new_cls, and_(old_cls.seq_id == new_cls.seq_id, old_cls.offset == new_cls.offset))\
                .filter(new_cls.id==None).all()


def main():
    motifs = get_old_wrong_new_exists()
    print(str.format("Get wrong in old and still exists in new, Count {}, {}", len(motifs), motifs))
    for old_m, new_m in motifs:
        print(str.format("Mark mid {}, seq_id {}, offset {} as false_discovery in version {}",
                         new_m.id, new_m.seq_id, new_m.offset, NEW_VERSION))
        #dao.update_false_discovery_by_motif(new_m.id, True, NEW_VERSION)

    matrix = pssm_matrix.calc_pssm_matrix(dao.find_baseline_motifs())
    motifs = get_old_correct_new_not_exists()
    with dao.query_session() as session:
        all_seq_ids = set([m.seq_id for m in motifs])
        seqs = dao.sequence.find_seq_by_ids(session, all_seq_ids)
        seq_ids_to_seq = dict([(seq.seq_id, seq) for seq in seqs])
        assert(len(seq_ids_to_seq) == len(all_seq_ids))

    print(str.format("Get correct in old and not exists in new, Count {}, {}", len(motifs), motifs))
    for m in motifs:
        motif_seq = seq_ids_to_seq[m.seq_id].seq[m.offset:m.offset+16]
        assert(len(motif_seq) == 16)
        score = motif_tools.calc_pssm_score(motif_seq, matrix)
        probability = motif_tools.calc_probability_by_score(score)
        new_m = dao.motif.MotifEntityBase(m.offset, m.seq_id, score, probability, 0, manually_add=True)
        new_motifs = dao.find_motifs_by_seq_ids([m.seq_id, ], NEW_VERSION, with_wrong=False)[m.seq_id]
        new_motifs.append(new_m)
        new_motifs_no_overlap = motif_tools.found_no_overlapped_motifs(new_motifs, 16)
        if len(new_motifs) != len(new_motifs_no_overlap):
            logging.error(str.format("Overlap found in seq {}, {}/{}, new offset {}",
                                     m.seq_id, len(new_motifs), len(new_motifs_no_overlap), new_m.offset))
        else:
            logging.debug(str.format("Add motif offset {}, score {}, probability {} to seq {} manually",
                                    new_m.offset, new_m.score, new_m.probability, new_m.seq_id))
            dao.add_manually_motif(new_m.seq_id, new_m.offset, NEW_VERSION, new_m.score, new_m.probability)




if __name__ == "__main__":
    logging.basicConfig(level=logging.DEBUG)
    main()
