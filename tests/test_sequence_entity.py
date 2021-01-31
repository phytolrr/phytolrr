import logging
from dao.datasource import *
from dao.sequence_entity import *
import unittest
from sqlalchemy.orm import sessionmaker
from sqlalchemy.orm.exc import NoResultFound


def _set_up_baseline_seq(seq_ids):
    seqs = []
    for seq_id in seq_ids:
        seq = SequenceEntity()
        seq.seq_id = seq_id
        seq.seq = 'ABCDEFG'
        seq.baseline = True
        seqs.append(seq)
    with session_scope() as session:
        session.add_all(seqs)


class TestSequenceEntity(unittest.TestCase):
    def setUp(self):
        logging.getLogger('sqlalchemy').setLevel(logging.INFO)
        Base.metadata.drop_all(engine)
        Base.metadata.create_all(engine)
        self.session = Session()

        seqs = []
        for i in range(0, 10):
            seq = SequenceEntity()
            seq.seq_id = str.format('ABCDEFG_{}', i)
            seq.seq = str.format('ABCDEFGSEQ{}', i)
            seqs.append(seq)
        self.session.add_all(seqs)
        self.session.commit()

    def test_add_seq(self):
        add_seq(self.session, 'SEQID1', 'SEQABCDEFG')
        self.session.commit()
        seq = find_seq_by_id(self.session, 'SEQID1')
        self.assertIsNotNone(seq)
        self.assertEqual('SEQID1', seq.seq_id)
        self.assertEqual('SEQABCDEFG', seq.seq)

        add_seq(self.session, 'SEQID2', 'SEQABC', 'SEQSS', 'SEQSS8', 'SEQACC', 'SEQACC20')
        self.session.commit()
        seq = find_seq_by_id(self.session, 'SEQID2')
        self.assertIsNotNone(seq)
        self.assertEqual('SEQID2', seq.seq_id)
        self.assertEqual('SEQABC', seq.seq)
        self.assertEqual('SEQSS', seq.ss)
        self.assertEqual('SEQSS8', seq.ss8)
        self.assertEqual('SEQACC', seq.acc)
        self.assertEqual('SEQACC20', seq.acc20)

    def test_seq_exists(self):
        self.assertTrue(seq_exists(self.session, 'ABCDEFG_0'))
        self.assertTrue(seq_exists(self.session, 'ABCDEFG_1'))
        self.assertTrue(seq_exists(self.session, 'ABCDEFG_5'))
        self.assertTrue(seq_exists(self.session, 'ABCDEFG_9'))
        self.assertFalse(seq_exists(self.session, 'ABCDEFG_10'))

    def test_find_seq_by_id(self):
        self.assertIsNotNone(find_seq_by_id(self.session, 'ABCDEFG_0'))
        try:
            find_seq_by_id(self.session, 'ABC')
        except NoResultFound:
            pass
        else:
            self.assertFalse(True, "The seq ABC should not exists")

    def test_update_seq_for_scratch_result(self):
        with session_scope() as session:
            update_seq_for_scratch_result(session, 'ABCDEFG_0', 'ABCSS', 'ABCSS8', 'ABCACC', 'ABCACC20')
        seq = find_seq_by_id(self.session, 'ABCDEFG_0')
        self.assertIsNotNone(seq)
        self.assertEqual('ABCSS', seq.ss)
        self.assertEqual('ABCSS8', seq.ss8)
        self.assertEqual('ABCACC', seq.acc)
        self.assertEqual('ABCACC20', seq.acc20)

    def test_find_all_seq_ids(self):
        with query_session() as session:
            seq_ids = find_all_seq_ids(session)
        self.assertEqual(10, len(seq_ids))
        expect_seq_ids = set([str.format('ABCDEFG_{}', i) for i in range(0,10)])
        self.assertEqual(expect_seq_ids, set(seq_ids))

    def test_find_seqs_by_baseline(self):
        baseline_seq_ids = ['SEQ1', 'SEQ2', 'SEQ3']
        _set_up_baseline_seq(baseline_seq_ids)
        with query_session() as session:
            seqs = find_seqs_by_baseline(session)
        self.assertEqual(3, len(seqs))
        self.assertListEqual(baseline_seq_ids, [seq.seq_id for seq in seqs])

    def test_update_subgroup_by_seq_ids(self):
        with session_scope() as session:
            update_subgroup_by_seq_ids(session, ['ABCDEFG_1', 'ABCDEFG_2', 'ABCDEFG_5'], '7-1')

        with query_session() as session:
            seqs = find_all_seqs(session)
            seq_ids_to_seq = dict([(seq.seq_id, seq) for seq in seqs])

        seq = seq_ids_to_seq['ABCDEFG_1']
        self.assertEqual('7-1', seq.subgroup)
        self.assertEqual('ABCDEFGSEQ1', seq.seq)

        seq = seq_ids_to_seq['ABCDEFG_2']
        self.assertEqual('7-1', seq.subgroup)
        self.assertEqual('ABCDEFGSEQ2', seq.seq)

        seq = seq_ids_to_seq['ABCDEFG_3']
        self.assertIsNone(seq.subgroup)
        self.assertEqual('ABCDEFGSEQ3', seq.seq)

        seq = seq_ids_to_seq['ABCDEFG_5']
        self.assertEqual('7-1', seq.subgroup)
        self.assertEqual('ABCDEFGSEQ5', seq.seq)


    def test_find_seqs_by_subgroups(self):
        with session_scope() as session:
            update_subgroup_by_seq_ids(session, ['ABCDEFG_1', 'ABCDEFG_2', 'ABCDEFG_5'], '7-1')
            update_subgroup_by_seq_ids(session, ['ABCDEFG_0', 'ABCDEFG_6', 'ABCDEFG_8'], '10')
            update_subgroup_by_seq_ids(session, ['ABCDEFG_8', 'ABCDEFG_7'], '11')

        with query_session() as session:
            seqs = find_seqs_by_subgroups(session, ['7-1', '11'])
        self.assertEqual(5, len(seqs))
        seq_ids_to_seq = dict([(seq.seq_id, seq) for seq in seqs])
        self.assertSetEqual({'ABCDEFG_1', 'ABCDEFG_2', 'ABCDEFG_5', 'ABCDEFG_8', 'ABCDEFG_7'}, set(seq_ids_to_seq.keys()))
        self.assertEqual('7-1', seq_ids_to_seq['ABCDEFG_1'].subgroup)
        self.assertEqual('7-1', seq_ids_to_seq['ABCDEFG_2'].subgroup)
        self.assertEqual('7-1', seq_ids_to_seq['ABCDEFG_5'].subgroup)
        self.assertEqual('11', seq_ids_to_seq['ABCDEFG_7'].subgroup)
        self.assertEqual('11', seq_ids_to_seq['ABCDEFG_8'].subgroup)


if __name__ == "__main__":
    unittest.main()
