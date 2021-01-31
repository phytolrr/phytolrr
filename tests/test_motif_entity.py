import unittest
from sqlalchemy import create_engine
from sqlalchemy.orm import sessionmaker
import logging
from dao.motif_entity import *
from dao.sequence_entity import *


MOTIF_VERSION = 2


class TestMotifEntity(unittest.TestCase):
    def setUp(self):
        logging.getLogger('sqlalchemy').setLevel(logging.DEBUG)
        Base.metadata.drop_all(engine)
        Base.metadata.create_all(engine)
        self.session = Session()
        self._add_seq_motifs('SEQ1', 10)
        self._add_seq_motifs('SEQ2', 20)
        self._add_seq_motifs('SEQ3', 10)
        self._add_seq_motifs('SEQ4', 20)
        self._add_seq_motifs('SEQ5', 10)
        self.session.commit()

    def _add_seq_motifs(self, seq_id, count):
        motifs = []
        seqs = []
        for i in range(0,count):
            motif = MotifEntityBase(i*100, seq_id, 0.1*i, 0.0, 0.0)
            motif.probability = 0.000001 * i
            motif.fdr_probability = 0.0001 * i
            if i == 0:
                motif.false_discovery = True
                motif.correct = False
            else:
                motif.false_discovery = False
                motif.correct = True
            motifs.append(motif)

        add_motifs(self.session, motifs, MOTIF_VERSION)
        add_seq(self.session, seq_id, 'SEQABC')
        self.session.commit()

    def test_find_motifs_by_seq_id(self):
        motifs = find_motifs_by_seq_id(self.session, 'SEQ1', MOTIF_VERSION)
        self.assertEqual(10, len(motifs))
        motifs = find_motifs_by_seq_id(self.session, 'SEQ4', MOTIF_VERSION)
        self.assertEqual(20, len(motifs))
        motifs = find_motifs_by_seq_id(self.session, 'SEQ1', MOTIF_VERSION, with_wrong=False)
        self.assertEqual(9, len(motifs))
        motifs = find_motifs_by_seq_id(self.session, 'SEQ4', MOTIF_VERSION, with_wrong=False)
        self.assertEqual(19, len(motifs))

    def test_find_motifs_count_by_seq_id(self):
        self.assertEqual(9, find_motifs_count_by_seq_id(self.session, 'SEQ1', MOTIF_VERSION))
        self.assertEqual(19, find_motifs_count_by_seq_id(self.session, 'SEQ2', MOTIF_VERSION))
        self.assertEqual(9, find_motifs_count_by_seq_id(self.session, 'SEQ3', MOTIF_VERSION))
        self.assertEqual(9, find_motifs_count_by_seq_id(self.session, 'SEQ5', MOTIF_VERSION))

    def test_find_seq_ids_by_motif_ids_when_not_correct(self):
        self.assertEqual(9, find_motifs_count_by_seq_id(self.session, 'SEQ1', MOTIF_VERSION))
        motifs = find_motifs_by_seq_id(self.session, 'SEQ1', MOTIF_VERSION)
        motifs[0].correct = False
        motifs[1].correct = False
        self.session.commit()
        motifs = find_motifs_by_seq_id(self.session, 'SEQ1', MOTIF_VERSION)
        self.assertEqual(8, find_motifs_count_by_seq_id(self.session, 'SEQ1', MOTIF_VERSION))


    def test_motif_exists(self):
        motifs = find_motifs_by_seq_id(self.session, 'SEQ1', MOTIF_VERSION)
        for m in motifs:
            self.assertTrue(motif_exists(self.session, m.id, MOTIF_VERSION))
            self.assertFalse(motif_exists(self.session, 1024, MOTIF_VERSION))

    def test_find_seq_ids_by_motif_ids(self):
        ids = []
        motifs = find_motifs_by_seq_id(self.session, 'SEQ1', MOTIF_VERSION)
        ids.append(motifs[1].id)
        motifs = find_motifs_by_seq_id(self.session, 'SEQ2', MOTIF_VERSION)
        ids.append(motifs[0].id)
        motifs = find_motifs_by_seq_id(self.session, 'SEQ4', MOTIF_VERSION)
        ids.append(motifs[18].id)
        ids.append(motifs[0].id)
        ids.append(motifs[10].id)
        seq_ids = find_seq_ids_by_motif_ids(self.session, ids, MOTIF_VERSION)
        self.assertEqual(3, len(seq_ids))
        self.assertEqual(set(['SEQ1', 'SEQ2', 'SEQ4']), set(seq_ids))

    def test_delete_motifs_by_seq_id(self):
        delete_motifs_by_seq_ids(self.session, ['SEQ1', 'SEQ4'], MOTIF_VERSION)
        self.session.commit()
        self.assertEqual(0, find_motifs_count_by_seq_id(self.session, 'SEQ1', MOTIF_VERSION))
        self.assertEqual(19, find_motifs_count_by_seq_id(self.session, 'SEQ2', MOTIF_VERSION))
        self.assertEqual(9, find_motifs_count_by_seq_id(self.session, 'SEQ3', MOTIF_VERSION))
        self.assertEqual(0, find_motifs_count_by_seq_id(self.session, 'SEQ4', MOTIF_VERSION))
        self.assertEqual(9, find_motifs_count_by_seq_id(self.session, 'SEQ5', MOTIF_VERSION))

    def test_replace_motifs_by_seq_add_and_remove(self):
        motifs = []
        motif = MotifEntityBase(100, 'SEQ1', 1000.0, 0.0, 0.0)
        motif.probability = 0.1
        motif.fdr_probability = 0.1
        motifs.append(motif)
        motif = MotifEntityBase(150, 'SEQ1', 1000.0, 0.0, 0.0)
        motif.probability = 0.1
        motif.fdr_probability = 0.1
        motifs.append(motif)
        replace_motifs_by_seq(self.session, 'SEQ1', motifs, MOTIF_VERSION)
        self.session.commit()
        seqs = find_motifs_by_seq_id(self.session, 'SEQ1', MOTIF_VERSION)
        self.assertEqual(2, len(seqs))
        seq = seqs[0]
        self.assertEqual('SEQ1', seq.seq_id)
        self.assertEqual(100, seq.offset)
        self.assertEqual(1000.0, seq.score)
        self.assertEqual(0.1, seq.probability)
        self.assertEqual(0.1, seq.fdr_probability)
        seq = seqs[1]
        self.assertEqual('SEQ1', seq.seq_id)
        self.assertEqual(150, seq.offset)
        self.assertEqual(1000.0, seq.score)
        self.assertEqual(0.1, seq.probability)
        self.assertEqual(0.1, seq.fdr_probability)

    def test_replace_motifs_by_seq_empty(self):
        replace_motifs_by_seq(self.session, 'SEQ2', [], MOTIF_VERSION)
        self.session.commit()
        self.assertEqual(0, find_motifs_count_by_seq_id(self.session, 'SEQ2', MOTIF_VERSION))

    def test_replace_motifs_by_seq_add(self):
        motifs_in_db = find_motifs_by_seq_id(self.session, 'SEQ3', MOTIF_VERSION)
        motifs = []
        for motif in motifs_in_db:
            new_motif = MotifEntityBase(motif.offset, motif.seq_id, motif.score, 0.0, 0.0, correct=motif.correct)
            new_motif.probability = motif.probability
            new_motif.fdr_probability = motif.fdr_probability
            motifs.append(new_motif)
        motif = MotifEntityBase(2000, 'SEQ3', 20.0, 0.0, 0.0)
        motif.probability = 0.1
        motif.fdr_probability = 0.1
        motifs.append(motif)
        expect_len = len(motifs)
        replace_motifs_by_seq(self.session, 'SEQ3', motifs, MOTIF_VERSION)
        self.session.commit()
        motifs = find_motifs_by_seq_id(self.session, 'SEQ3', MOTIF_VERSION)
        self.assertEqual(expect_len, len(motifs))

    def test_find_motifs_by_seq_ids(self):
        seq_ids_to_motifs = find_motifs_by_seq_ids(self.session, {'SEQ1', 'SEQ2', 'SEQ3'}, MOTIF_VERSION)
        self.assertEqual(3, len(seq_ids_to_motifs))
        self.assertEqual(10, len(seq_ids_to_motifs['SEQ1']))
        self.assertEqual(20, len(seq_ids_to_motifs['SEQ2']))
        self.assertEqual(10, len(seq_ids_to_motifs['SEQ3']))

        seq_ids_to_motifs = find_motifs_by_seq_ids(self.session, {'SEQ1', 'SEQ2', 'SEQ3'}, MOTIF_VERSION, False)
        self.assertEqual(3, len(seq_ids_to_motifs))
        self.assertEqual(9, len(seq_ids_to_motifs['SEQ1']))
        self.assertEqual(19, len(seq_ids_to_motifs['SEQ2']))
        self.assertEqual(9, len(seq_ids_to_motifs['SEQ3']))

    def test_delete_manually_motifs_by_ids(self):
        # 为SEQ1添加3个手工添加的motif
        ms = []
        m = MotifEntityBase(1000, 'SEQ1', 10.0, 0.0, 0.0, manually_add=True)
        ms.append(m)
        m = MotifEntityBase(1020, 'SEQ1', 12.0, 0.0, 0.0, manually_add=True)
        ms.append(m)
        m = MotifEntityBase(1040, 'SEQ1', 13.0, 0.0, 0.0, manually_add=True)
        ms.append(m)
        with session_scope() as session:
            add_motifs(session, ms, MOTIF_VERSION)

        # 检查确认添加成功
        with query_session() as session:
            motifs = find_motifs_by_seq_id(session, 'SEQ1', MOTIF_VERSION)
        self.assertEqual(13, len(motifs))
        offsets_to_motif = dict([(m.offset, m) for m in motifs])
        self.assertTrue(1000 in offsets_to_motif)
        self.assertTrue(1020 in offsets_to_motif)
        m1000 = offsets_to_motif[1000]
        m1040 = offsets_to_motif[1040]
        self.assertTrue(m1000.manually_add)
        self.assertFalse(offsets_to_motif[100].manually_add)

        # 删除其中两个手工添加的motif
        with session_scope() as session:
            ret = delete_manually_motifs_by_ids(session, [m1000.id, m1040.id], MOTIF_VERSION)

        # 检查确认删除成功、检查其他的motif没有被删除
        self.assertEqual(2, ret)
        with query_session() as session:
            motifs = find_motifs_by_seq_id(session, 'SEQ1', MOTIF_VERSION)
        self.assertEqual(11, len(motifs))
        offsets_to_motif = dict([(m.offset, m) for m in motifs])
        self.assertFalse(1000 in offsets_to_motif)
        self.assertTrue(1020 in offsets_to_motif)
        self.assertTrue(100 in offsets_to_motif)
        self.assertTrue(offsets_to_motif[1020].manually_add)
        self.assertFalse(offsets_to_motif[100].manually_add)

    def test_delete_manually_motifs_by_ids_not_manually(self):
        # 为SEQ1添加两个手工添加的motif
        ms = []
        m = MotifEntityBase(1000, 'SEQ1', 10.0, 0.0, 0.0, manually_add=True)
        ms.append(m)
        m = MotifEntityBase(1020, 'SEQ1', 12.0, 0.0, 0.0, manually_add=True)
        ms.append(m)
        with session_scope() as session:
            add_motifs(session, ms, MOTIF_VERSION)

        # 检查确认添加成功
        with query_session() as session:
            motifs = find_motifs_by_seq_id(session, 'SEQ1', MOTIF_VERSION)
        self.assertEqual(12, len(motifs))
        offsets_to_motif = dict([(m.offset, m) for m in motifs])
        self.assertTrue(1000 in offsets_to_motif)
        self.assertTrue(1020 in offsets_to_motif)
        m1000 = offsets_to_motif[1000]
        m100 = offsets_to_motif[100]
        self.assertTrue(m1000.manually_add)
        self.assertFalse(m100.manually_add)

        # 删除其中一个手工添加的motif
        with session_scope() as session:
            ret = delete_manually_motifs_by_ids(session, [m100.id, m1000.id], MOTIF_VERSION)

        # 检查确认删除成功、检查其他的motif没有被删除
        self.assertEqual(1, ret)
        with query_session() as session:
            motifs = find_motifs_by_seq_id(session, 'SEQ1', MOTIF_VERSION)
        self.assertEqual(11, len(motifs))
        offsets_to_motif = dict([(m.offset, m) for m in motifs])
        self.assertFalse(1000 in offsets_to_motif)
        self.assertTrue(1020 in offsets_to_motif)
        self.assertTrue(100 in offsets_to_motif)
        self.assertTrue(offsets_to_motif[1020].manually_add)
        self.assertFalse(offsets_to_motif[100].manually_add)

    def test_update_false_discovery_by_motif(self):
        with query_session() as session:
            motifs = find_motifs_by_seq_id(session, 'SEQ2', MOTIF_VERSION)
        self.assertTrue(motifs[0].false_discovery)
        self.assertFalse(motifs[0].correct)
        self.assertFalse(motifs[1].false_discovery)
        self.assertTrue(motifs[1].correct)
        self.assertFalse(motifs[5].false_discovery)
        self.assertTrue(motifs[5].correct)

        # 标注0的false_discovery为False
        # 标注5的fals_discovery为True
        with session_scope() as session:
            update_false_discovery_by_motif(session, motifs[0].id, False, MOTIF_VERSION)
            update_false_discovery_by_motif(session, motifs[5].id, True, MOTIF_VERSION)

        with query_session() as session:
            motifs = find_motifs_by_seq_id(session, 'SEQ2', MOTIF_VERSION)
        self.assertTrue(motifs[5].false_discovery)
        self.assertFalse(motifs[5].correct)
        self.assertFalse(motifs[1].false_discovery)
        self.assertTrue(motifs[1].correct)
        self.assertFalse(motifs[0].false_discovery)
        self.assertTrue(motifs[0].correct)

    def test_update_score(self):
        with query_session() as session:
            motifs = find_motifs_by_seq_id(session, 'SEQ1', MOTIF_VERSION)
        self.assertLess(motifs[1].score, 0.11)
        self.assertGreater(motifs[1].score, 0.09)
        self.assertLess(motifs[2].score, 0.21)
        self.assertGreater(motifs[2].score, 0.19)
        with session_scope() as session:
            update_score_by_motif(session, motifs[1].id, 100.0, MOTIF_VERSION)
        with query_session() as session:
            motifs = find_motifs_by_seq_id(session, 'SEQ1', MOTIF_VERSION)
        self.assertLess(motifs[1].score, 100.01)
        self.assertGreater(motifs[1].score, 99.99)
        self.assertLess(motifs[2].score, 0.21)
        self.assertGreater(motifs[2].score, 0.19)

if __name__ == "__main__":
    unittest.main()
