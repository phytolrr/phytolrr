import dao
from tests.test_tag_ref_entity import *
from tests.test_motif_entity import *
from tests.test_nsite_entity import *
from tests.test_sequence_entity import *

MOTIF_VERSION = 1


class TestDao(unittest.TestCase):
    def setUp(self):
        logging.getLogger('sqlalchemy').setLevel(logging.DEBUG)
        dao.Base.metadata.drop_all(dao.engine)
        dao.Base.metadata.create_all(dao.engine)

        self._set_up_seq('SEQ1', 'SEQABCDEFG', 'SP1')
        self._set_up_seq('SEQ2', 'SEQABCDEFG', 'SP1')
        self._set_up_seq('SEQ3', 'SEQABCDEFG', 'SP1')
        self._set_up_seq('SEQ4', 'SEQABCDEFG', 'SP2', start=5, count=20)
        self._set_up_seq('SEQ5', 'SEQABCDEFG', 'SP2', start=5, count=20)
        self._set_up_seq('SEQ6', 'SEQABCDEFG', 'SP2', start=5, count=20)
        self._set_up_seq('SEQ7', 'SEQABCDEFG', 'SP3', start=10, count=30)
        self._set_up_seq('SEQ8', 'SEQABCDEFG', 'SP3', start=10, count=30)
        self._set_up_seq('SEQ9', 'SEQABCDEFG', 'SP3', start=10, count=30)

    def _set_up_seq(self, seq_id, seq_str, species, start=0, count=10):
        motifs = []
        for i in range(start, count + start):
            m = dao.motif.MotifEntityBase(i * 20, seq_id, 0.1 * i, 0.0, 0.0)
            motifs.append(m)
        dao.add_seq(seq_id, seq_str, MOTIF_VERSION, motifs, species)

    def test_add_seq(self):
        seq = dao.find_seq_by_id('SEQ1')
        self.assertIsNotNone(seq)
        self.assertEqual(10, len(dao.find_correct_motifs_by_seq_ids(['SEQ1', ], MOTIF_VERSION)['SEQ1']))
        self.assertEqual('SEQABCDEFG', seq.seq)

    def test_replace_tags_by_motifs(self):
        seq_ids_to_motifs = {'SEQ1': [], 'SEQ2': []}
        with query_session() as session:
            seq_ids_to_motifs['SEQ1'] = dao.motif.find_motifs_by_seq_id(session, 'SEQ1', MOTIF_VERSION)
            seq_ids_to_motifs['SEQ2'] = dao.motif.find_motifs_by_seq_id(session, 'SEQ2', MOTIF_VERSION)
        mids_to_tags = {seq_ids_to_motifs['SEQ1'][0].id: ['a','b','c']}
        cls = dao.motif.get_entity(MOTIF_VERSION)
        with query_session() as session:
            motifs = session.query(cls).filter(cls.id.in_(mids_to_tags.keys())).all()
        self.assertEqual(1, len(motifs))
        self.assertTrue(motifs[0].correct)
        dao.replace_tags_by_motifs(mids_to_tags, MOTIF_VERSION)
        seq = dao.find_seq_by_id('SEQ1')
        self.assertIsNotNone(seq)
        self.assertEqual(10, len(dao.find_correct_motifs_by_seq_ids(['SEQ1', ], MOTIF_VERSION)['SEQ1']))
        with query_session() as session:
            motifs = session.query(cls).filter(cls.id.in_(mids_to_tags.keys())).all()
        self.assertEqual(1, len(motifs))
        self.assertTrue(motifs[0].correct)

    def unused_test_replace_tags_by_motifs_false_discovery(self):
        # 打tag前，SEQ2有10个motif
        seq_ids_to_motifs = dao.find_motifs_by_seq_ids(['SEQ2',], MOTIF_VERSION)
        self.assertTrue('SEQ2' in seq_ids_to_motifs)
        motifs = seq_ids_to_motifs['SEQ2']
        self.assertEqual(10, len(motifs))

        # 为第0个motif打标记错误的tag

        dao.replace_tags_by_motifs({motifs[0].id: ['a','inner.falsediscovery']}, MOTIF_VERSION)
        # 检查点：第0个motif重查出来后，correct字段为False
        with query_session() as session:
            motif_0_in = dao.motif.find_motifs_by_ids(session, [motifs[0].id], MOTIF_VERSION)
        self.assertEqual(1, len(motif_0_in))
        self.assertFalse(motif_0_in[0].correct)
        # 检查点：SEQ2带错查询，查出10个motif
        seq_ids_to_motifs = dao.find_motifs_by_seq_ids(['SEQ2',], MOTIF_VERSION)
        self.assertTrue('SEQ2' in seq_ids_to_motifs)
        motifs_after_tag = seq_ids_to_motifs['SEQ2']
        self.assertEqual(10, len(motifs_after_tag))
        # 检查点：SEQ2本身的motif个数为9个
        self.assertEqual(9, len(dao.find_correct_motifs_by_seq_ids(['SEQ2', ], MOTIF_VERSION)['SEQ2']))

        # 为第5个motif打标记重叠的tag，SEQ2有8个motif
        dao.replace_tags_by_motifs({motifs[5].id: ['a','inner.overlap']}, MOTIF_VERSION)
        # 检查点：不带错查询，能查出8个motif
        seq_ids_to_motifs = dao.find_correct_motifs_by_seq_ids(['SEQ2',], MOTIF_VERSION)
        self.assertTrue('SEQ2' in seq_ids_to_motifs)
        motifs_after_tag = seq_ids_to_motifs['SEQ2']
        self.assertEqual(8, len(motifs_after_tag))

        # 为第5个motif replace回不相干的tag，SEQ2有9个correct motif
        dao.replace_tags_by_motifs({motifs[5].id: ['a']}, MOTIF_VERSION)
        self.assertEqual(1, len(motif_0_in))
        self.assertFalse(motif_0_in[0].correct)
        # 检查点：SEQ2带错查询，查出10个motif
        seq_ids_to_motifs = dao.find_motifs_by_seq_ids(['SEQ2',], MOTIF_VERSION)
        self.assertTrue('SEQ2' in seq_ids_to_motifs)
        motifs_after_tag = seq_ids_to_motifs['SEQ2']
        self.assertEqual(10, len(motifs_after_tag))
        # 检查点：SEQ2本身的motif个数为9个
        self.assertEqual(9, len(dao.find_correct_motifs_by_seq_ids(['SEQ2', ], MOTIF_VERSION)['SEQ2']))

    def test_add_tags_by_names_to_ids(self):
        # 打tag前，SEQ2有10个motif
        seq_ids_to_motifs = dao.find_motifs_by_seq_ids(['SEQ2',], MOTIF_VERSION)
        self.assertTrue('SEQ2' in seq_ids_to_motifs)
        motifs = seq_ids_to_motifs['SEQ2']
        self.assertEqual(10, len(motifs))

        # 为第0个motif打重复的tag
        dao.add_tags_by_names_to_ids({'inner.overlap': [motifs[0].id]}, MOTIF_VERSION)
        # 检查点：第0个motif重查出来后，correct字段为False
        with query_session() as session:
            motif_0_in = dao.motif.find_motifs_by_ids(session, [motifs[0].id], MOTIF_VERSION)
        self.assertEqual(1, len(motif_0_in))
        self.assertFalse(motif_0_in[0].correct)
        # 检查点：SEQ2带错查询，查出10个motif
        seq_ids_to_motifs = dao.find_motifs_by_seq_ids(['SEQ2',], MOTIF_VERSION)
        self.assertTrue('SEQ2' in seq_ids_to_motifs)
        motifs_after_tag = seq_ids_to_motifs['SEQ2']
        self.assertEqual(10, len(motifs_after_tag))
        # 检查点：SEQ2本身的motif个数为9个
        self.assertEqual(9, len(dao.find_correct_motifs_by_seq_ids(['SEQ2', ], MOTIF_VERSION)['SEQ2']))

        # 为第0个motif添加无所谓的tag
        dao.add_tags_by_names_to_ids({'taga': [motifs[0].id]}, MOTIF_VERSION)
        # 检查点：第0个motif重查出来后，correct字段为False
        with query_session() as session:
            motif_0_in = dao.motif.find_motifs_by_ids(session, [motifs[0].id], MOTIF_VERSION)
        self.assertEqual(1, len(motif_0_in))
        self.assertFalse(motif_0_in[0].correct)
        # 检查点：SEQ2带错查询，查出10个motif
        seq_ids_to_motifs = dao.find_motifs_by_seq_ids(['SEQ2',], MOTIF_VERSION)
        self.assertTrue('SEQ2' in seq_ids_to_motifs)
        motifs_after_tag = seq_ids_to_motifs['SEQ2']
        self.assertEqual(10, len(motifs_after_tag))
        # 检查点: SEQ2只查正确，查出9个motif
        seq_ids_to_motifs = dao.find_correct_motifs_by_seq_ids(['SEQ2',], MOTIF_VERSION)
        self.assertTrue('SEQ2' in seq_ids_to_motifs)
        motifs_after_tag = seq_ids_to_motifs['SEQ2']
        self.assertEqual(9, len(motifs_after_tag))
        # 检查点：SEQ2本身的motif个数为9个
        self.assertEqual(9, len(dao.find_correct_motifs_by_seq_ids(['SEQ2', ], MOTIF_VERSION)['SEQ2']))

    def test_add_tags_by_names_to_ids_two_mark_on_one_motif(self):
        # 打tag前，SEQ2有10个motif
        seq_ids_to_motifs = dao.find_motifs_by_seq_ids(['SEQ2',], MOTIF_VERSION)
        self.assertTrue('SEQ2' in seq_ids_to_motifs)
        motifs = seq_ids_to_motifs['SEQ2']
        self.assertEqual(10, len(motifs))

        # 为第0个motif打重复的tag
        dao.add_tags_by_names_to_ids({'inner.overlap': [motifs[0].id]}, MOTIF_VERSION)
        # 检查点：第0个motif重查出来后，correct字段为False
        with query_session() as session:
            motif_0_in = dao.motif.find_motifs_by_ids(session, [motifs[0].id], MOTIF_VERSION)
        self.assertEqual(1, len(motif_0_in))
        self.assertFalse(motif_0_in[0].correct)
        # 检查点：SEQ2带错查询，查出10个motif
        seq_ids_to_motifs = dao.find_motifs_by_seq_ids(['SEQ2',], MOTIF_VERSION)
        self.assertTrue('SEQ2' in seq_ids_to_motifs)
        motifs_after_tag = seq_ids_to_motifs['SEQ2']
        self.assertEqual(10, len(motifs_after_tag))
        # 检查点：SEQ2本身的motif个数为9个
        self.assertEqual(9, len(dao.find_correct_motifs_by_seq_ids(['SEQ2', ], MOTIF_VERSION)['SEQ2']))

        # 为第0个motif添加手工标记错误的tag
        dao.add_tags_by_names_to_ids({'inner.falsediscovery': [motifs[0].id]}, MOTIF_VERSION)
        # 检查点：第0个motif重查出来后，correct字段为False
        with query_session() as session:
            motif_0_in = dao.motif.find_motifs_by_ids(session, [motifs[0].id], MOTIF_VERSION)
        self.assertEqual(1, len(motif_0_in))
        self.assertFalse(motif_0_in[0].correct)
        # 检查点：SEQ2带错查询，查出10个motif
        seq_ids_to_motifs = dao.find_motifs_by_seq_ids(['SEQ2',], MOTIF_VERSION)
        self.assertTrue('SEQ2' in seq_ids_to_motifs)
        motifs_after_tag = seq_ids_to_motifs['SEQ2']
        self.assertEqual(10, len(motifs_after_tag))
        # 检查点: SEQ2只查正确，查出9个motif
        seq_ids_to_motifs = dao.find_correct_motifs_by_seq_ids(['SEQ2',], MOTIF_VERSION)
        self.assertTrue('SEQ2' in seq_ids_to_motifs)
        motifs_after_tag = seq_ids_to_motifs['SEQ2']
        self.assertEqual(9, len(motifs_after_tag))
        # 检查点：SEQ2本身的motif个数为9个
        self.assertEqual(9, len(dao.find_correct_motifs_by_seq_ids(['SEQ2', ], MOTIF_VERSION)['SEQ2']))

    def test_query_sequences(self):
        filters = {'page_index': 0, 'page_size': 20}
        seqs, total = dao.query_sequences(filters, MOTIF_VERSION)
        self.assertEqual(9, len(seqs))
        self.assertEqual(9, total)

    def test_query_sequences_offset(self):
        filters = {'page_index': 0, 'page_size': 20, 'offsetlt': 40}
        seqs, total = dao.query_sequences(filters, MOTIF_VERSION)
        self.assertEqual(3, len(seqs))
        self.assertEqual(3, total)

        filters = {'page_index': 0, 'page_size': 20, 'offsetlt': 101}
        seqs, total = dao.query_sequences(filters, MOTIF_VERSION)
        self.assertEqual(6, len(seqs))
        self.assertEqual(6, total)

        filters = {'page_index': 0, 'page_size': 20, 'offsetgt': 600}
        seqs, total = dao.query_sequences(filters, MOTIF_VERSION)
        self.assertEqual(3, len(seqs))
        self.assertEqual(3, total)

    def test_query_sequences_lrr_count(self):
        filters = {'page_index': 0, 'page_size': 20, 'lrr_count_lt': 15}
        seqs, total = dao.query_sequences(filters, MOTIF_VERSION)
        self.assertEqual(3, len(seqs))
        self.assertEqual(3, total)

        filters = {'page_index': 0, 'page_size': 20, 'lrr_count_gt': 29}
        seqs, total = dao.query_sequences(filters, MOTIF_VERSION)
        self.assertEqual(3, len(seqs))
        self.assertEqual(3, total)

        filters = {'page_index': 0, 'page_size': 20, 'lrr_count_lt': 21, 'lrr_count_gt': 19}
        seqs, total = dao.query_sequences(filters, MOTIF_VERSION)
        self.assertEqual(3, len(seqs))
        self.assertEqual(3, total)

    def test_query_sequences_keyword(self):
        filters = {'page_index': 0, 'page_size': 20, 'keyword': 'SEQ1'}
        seqs, total = dao.query_sequences(filters, MOTIF_VERSION)
        self.assertEqual(1, len(seqs))
        self.assertEqual(1, total)

        filters = {'page_index': 0, 'page_size': 20, 'keyword': 'SEQ'}
        seqs, total = dao.query_sequences(filters, MOTIF_VERSION)
        self.assertEqual(9, len(seqs))
        self.assertEqual(9, total)

    def test_query_sequences_species(self):
        filters = {'page_index': 0, 'page_size': 20, 'species': ['SP1', 'SP3']}
        seqs, total = dao.query_sequences(filters, MOTIF_VERSION)
        self.assertEqual(6, len(seqs))
        self.assertEqual(6, total)

        filters = {'page_index': 0, 'page_size': 20, 'species': ['SP2',]}
        seqs, total = dao.query_sequences(filters, MOTIF_VERSION)
        self.assertEqual(3, len(seqs))
        self.assertEqual(3, total)

    def test_add_motif(self):
        # 向seq中添加motif，OK
        self.assertIsNone(dao.add_manually_motif('SEQ1', 400, MOTIF_VERSION, 10.0, 0.1))
        self.assertIsNone(dao.add_manually_motif('SEQ1', 420, MOTIF_VERSION, 20.0, 0.2))
        seq_ids_to_motifs = dao.find_motifs_by_seq_ids(['SEQ1', 'SEQ2'], MOTIF_VERSION)
        self.assertEqual(2, len(seq_ids_to_motifs))
        self.assertEqual({'SEQ1', 'SEQ2'}, seq_ids_to_motifs.keys())
        seq1_motifs = seq_ids_to_motifs['SEQ1']
        self.assertEqual(12, len(seq1_motifs))
        self.assertEqual(10, len(seq_ids_to_motifs['SEQ2']))
        offsets_to_motif = dict([(m.offset, m) for m in seq1_motifs])
        m1 = offsets_to_motif.get(400, None)
        self.assertIsNotNone(m1)
        self.assertEqual(400, m1.offset)
        self.assertTrue(10.001 > m1.score)
        self.assertTrue(9.999 < m1.score)
        m2 = offsets_to_motif.get(420, None)
        self.assertIsNotNone(m2)
        self.assertEqual(420, m2.offset)
        self.assertTrue(20.001 > m2.score)
        self.assertTrue(19.999 < m2.score)

if __name__ == "__main__":
    unittest.main()
