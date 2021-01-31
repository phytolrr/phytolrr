import logging
from dao.nsite_entity import *
import unittest


class TestNSiteEntity(unittest.TestCase):
    def setUp(self):
        logging.getLogger('sqlalchemy').setLevel(logging.DEBUG)
        Base.metadata.drop_all(engine)
        Base.metadata.create_all(engine)
        with session_scope() as session:
            replace_nsites_by_seq_id(session, 'SEQ1', [1,5,10,15,20], [S,S,S,T,T])
            replace_nsites_by_seq_id(session, 'SEQ2', [1,5,10,15,20], [S,S,S,T,T])
            replace_nsites_by_seq_id(session, 'SEQ3', [1,5,10,15,20], [S,S,S,T,T])

    def test_find_nsites_by_seq_id(self):
        with session_scope() as session:
            nsites = find_nsites_by_seq_id(session, 'SEQ1')
            self.assertEqual(5, len(nsites))
            nsites = find_nsites_by_seq_id(session, 'SEQ2')
            self.assertEqual(5, len(nsites))
            nsites = find_nsites_by_seq_id(session, 'SEQ3')
            self.assertEqual(5, len(nsites))
            nsites = find_nsites_by_seq_id(session, 'SEQ4')
            self.assertEqual(0, len(nsites))

    def test_find_nsites_count_by_seq_id(self):
        with query_session() as session:
            self.assertEqual(5, find_nsites_count_by_seq_id(session, 'SEQ1'))
            self.assertEqual(5, find_nsites_count_by_seq_id(session, 'SEQ2'))
            self.assertEqual(5, find_nsites_count_by_seq_id(session, 'SEQ3'))
            self.assertEqual(0, find_nsites_count_by_seq_id(session, 'SEQ4'))

    def test_replace_nsites_by_seq_id1(self):
        with session_scope() as session:
            replace_nsites_by_seq_id(session, 'SEQ2', [1,5,15,30], [S,T,S,S])
        with session_scope() as session:
            nsites = find_nsites_by_seq_id(session, 'SEQ2')
            self.assertEqual(4, len(nsites))
            ids = set([(n.start_pos, n.ntype) for n in nsites])
            self.assertEqual({(1, S), (5, T), (15, S), (30, S)}, ids)

    def test_replace_nsites_by_seq_id2(self):
        with session_scope() as session:
            replace_nsites_by_seq_id(session, 'SEQ2', [], [])
        with session_scope() as session:
            nsites = find_nsites_by_seq_id(session, 'SEQ2')
            self.assertEqual(0, len(nsites))

    def test_replace_nsites_by_seq_id3(self):
        ids = [1, 5, 10, 15, 20, 25, 30]
        ntypes = [S,S,S,T,T,S,T]
        with session_scope() as session:
            replace_nsites_by_seq_id(session, 'SEQ3', ids, ntypes)
        with session_scope() as session:
            nsites = find_nsites_by_seq_id(session, 'SEQ3')
            self.assertEqual(7, len(nsites))
            self.assertEqual(set(zip(ids, ntypes)), set([(n.start_pos, n.ntype) for n in nsites]))


    def test_find_nsites_by_seq_ids(self):
        with query_session() as session:
            seq_ids_to_nsites = find_nsites_by_seq_ids(session, ['SEQ1', 'SEQ2'])
        self.assertEqual(2, len(seq_ids_to_nsites))
        self.assertSetEqual({(1,S),(5,S),(10,S),(15,T),(20,T)}, set([(n.start_pos, n.ntype) for n in seq_ids_to_nsites['SEQ1']]))
        self.assertSetEqual({(1,S),(5,S),(10,S),(15,T),(20,T)}, set([(n.start_pos, n.ntype) for n in seq_ids_to_nsites['SEQ2']]))


if __name__ == "__main__":
    unittest.main()
