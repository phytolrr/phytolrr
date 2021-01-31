# -*- coding: utf-8 -*
import logging
import unittest
from dao.datasource import *
from dao.tag_ref_entity import *


class TestTagRef(unittest.TestCase):
    def setUp(self):
        logging.getLogger('sqlalchemy').setLevel(logging.INFO)
        self.session = Session()
        Base.metadata.drop_all(engine)
        Base.metadata.create_all(engine)
        target_ids_to_names = {}
        target_ids_to_names[10] = {'a','b','c','d'}
        target_ids_to_names[20] = {'a','b','c','e'}
        target_ids_to_names[30] = {'a','b','c','f'}
        target_ids_to_names[40] = {'a','b','c','g'}
        with session_scope() as session:
            add_tags_by_map(session, target_ids_to_names)

    def test_find_tags_by_motif_ids(self):
        session = self.session
        mids_to_tags = find_tags_by_motif_ids(session, {10,30})
        self.assertEqual(2, len(mids_to_tags))
        self.assertEqual({'a','b','c','d'}, set([tag.tag_name for tag in mids_to_tags[10]]))
        self.assertEqual({'a','b','c','f'}, set([tag.tag_name for tag in mids_to_tags[30]]))

    def test_find_tag_names_by_motif_ids(self):
        session = self.session
        mids_to_names = find_tag_names_by_motif_ids(session, {20, 40, 50})
        self.assertEqual(2, len(mids_to_names))
        self.assertEqual({'a', 'b', 'c', 'e'}, mids_to_names[20])
        self.assertEqual({'a', 'b', 'c', 'g'}, mids_to_names[40])

    def test_replace_tags_by_motifs1(self):
        target_ids_to_names = {}
        target_ids_to_names[10] = {}
        target_ids_to_names[20] = {}

        ids_to_tags = find_tags_by_motif_ids(self.session, {10,20,50})

        replace_tags_by_motifs(self.session, target_ids_to_names)
        self.session.commit()
        ids_to_names = find_tag_names_by_motif_ids(self.session, {10,20,50})
        self.assertEqual(0, len(ids_to_names), str(ids_to_names))

    def test_replace_tags_by_motifs2(self):
        target_ids_to_names = {}
        target_ids_to_names[10] = {'x','y'}
        target_ids_to_names[20] = {}
        replace_tags_by_motifs(self.session, target_ids_to_names)
        self.session.commit()
        ids_to_names = find_tag_names_by_motif_ids(self.session, {10,20,50})
        self.assertEqual(1, len(ids_to_names))
        self.assertEqual({'x', 'y'}, ids_to_names[10])

    def test_replace_tags_by_motifs3(self):
        target_ids_to_names = {}
        target_ids_to_names[30] = {'a','b','c','f', 'g'}
        replace_tags_by_motifs(self.session, target_ids_to_names)
        self.session.commit()
        ids_to_names = find_tag_names_by_motif_ids(self.session, {30})
        self.assertEqual(1, len(ids_to_names))
        self.assertEqual({'a','b','c','f', 'g'}, ids_to_names[30])

    def test_find_tag_refs_by_tag_name_and_motif_ids(self):
        tag_refs = find_tag_refs_by_tag_name_and_motif_ids(self.session, 'a', [10,20,30,50,70])
        except_motif_ids = {10, 20, 30}
        self.assertEqual(len(except_motif_ids), len(tag_refs))
        self.assertEqual(except_motif_ids, set([t.target_id for t in tag_refs]))

        tag_refs = find_tag_refs_by_tag_name_and_motif_ids(self.session, 'a', [50,70])
        except_motif_ids = set()
        self.assertEqual(len(except_motif_ids), len(tag_refs))
        self.assertEqual(except_motif_ids, set([t.target_id for t in tag_refs]))

    def test_add_tags_by_names_to_motifs_duplicate(self):
        # 检查点：初始状态：10，20，30，40都有tag a
        with query_session() as session:
            refs = find_tag_refs_by_tag_name_and_motif_ids(session, 'a', [10,20,30,40])
        self.assertEqual(4, len(refs))

        # 对10,20,30,40再打一次tag a，不会重复添加
        names_to_ids = {'a': [10, 20, 30, 40]}
        add_tags_by_names_to_motifs(session, names_to_ids)

        with query_session() as session:
            refs = find_tag_refs_by_tag_name_and_motif_ids(session, 'a', [10,20,30,40])
        self.assertEqual(4, len(refs))

    def test_add_tags_by_names_to_motifs(self):
        # 检查点：初始状态：10，20，30，40只有10有tag d
        with query_session() as session:
            refs = find_tag_refs_by_tag_name_and_motif_ids(session, 'd', [10,20,30,40])
        self.assertEqual(1, len(refs))

        # 对10,20,30,40打tag d
        names_to_ids = {'d': [10, 20, 30, 40]}
        with session_scope() as session:
            add_tags_by_names_to_motifs(session, names_to_ids)

        with query_session() as session:
            refs = find_tag_refs_by_tag_name_and_motif_ids(session, 'd', [10,20,30,40])
        self.assertEqual(4, len(refs))




if __name__ == "__main__":
    unittest.main()
