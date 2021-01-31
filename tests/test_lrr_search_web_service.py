# -*- coding: utf-8 -*
import logging
import unittest
import json
import random
import dao
from boddle import boddle

from web_service import lrr_search_web_service
from tools.exception import ErrorCode
from dao.sequence_entity import SequenceEntity


MOTIF_VERSION = 2
AMINOS = 'TWVILNGERPSMKHFQCAYD'


def random_seq(seq_len):
    seq5_seq = ''
    for i in range(seq_len):
        seq5_seq += random.choice(AMINOS)
    return seq5_seq


def set_up_baseline_seq(seq_id, seq_str, ss=None, start=0, count=0, step=20):
    seq = SequenceEntity()
    seq.seq_id = seq_id
    seq.seq = random_seq(3000)
    seq.baseline = True
    seq.ss = ss

    motifs = []
    for i in range(start, count + start):
        m = dao.motif.MotifEntityBase(i * step, seq_id, 0.1 * i, 10.0, 0.1)
        motifs.append(m)

    with dao.session_scope() as session:
        session.add(seq)
        dao.motif.add_motifs(session, motifs, 1)


def set_up_db():
    #logging.getLogger('sqlalchemy').setLevel(logging.DEBUG)
    dao.Base.metadata.drop_all(dao.engine)
    dao.Base.metadata.create_all(dao.engine)


def set_up_seq(seq_id, seq_str, ss=None, start=0, count=10, step=20):
    motifs = []
    for i in range(start, count + start):
        m = dao.motif.MotifEntityBase(i * step, seq_id, 0.1 * i, 0.0, 0.0)
        motifs.append(m)
    dao.add_seq(seq_id, seq_str, MOTIF_VERSION, motifs_16=motifs, ss=ss)


def mark_wrong(sid, mid):
    with boddle(json={'false_discovery': True}):
        return lrr_search_web_service.change_motif_false_discovery(MOTIF_VERSION, sid, mid)


def unmark_wrong(sid, mid):
    with boddle(json={'false_discovery': False}):
        result = lrr_search_web_service.change_motif_false_discovery(MOTIF_VERSION, sid, mid)
        return json.loads(result)


def delete_motif(mid):
    with boddle():
        return lrr_search_web_service.delete_manually_motif(MOTIF_VERSION, -1, str(mid))


class TestGetSequences(unittest.TestCase):
    def setUp(self):
        set_up_db()
        set_up_seq('SEQ2', random_seq(700), ss='CHECHECHECCC', start=10, count=20)
        set_up_seq('SEQ1', 'SEQ1ABCDABCD', ss='CHECHECHEHHH')
        set_up_seq('SEQ3', 'SEQ3ABCDABCD', ss='CHECHECHEEEE', start=50, count=30)

    def test_get_sequences(self):
        with boddle(query={"page": 0, "size": 20}):
            result = lrr_search_web_service.get_sequences(MOTIF_VERSION)
            self.assertIsNotNone(result)
            result = json.loads(result)
            self.assertIsNotNone(result.get("sequences"))
            self.assertEqual(3, len(result.get("sequences")), result)
            seq_ids_to_seq = dict([(seq['sequence_id'], seq) for seq in result.get('sequences')])
            self.assertSetEqual({'SEQ1', 'SEQ2', 'SEQ3'}, set(seq_ids_to_seq.keys()))
            self.assertEqual(seq_ids_to_seq['SEQ1']['seq'], 'SEQ1ABCDABCD')
            self.assertEqual(seq_ids_to_seq['SEQ1']['ss'], 'CHECHECHEHHH')
            self.assertEqual(len(seq_ids_to_seq['SEQ1']['motifs_16']), 10)
            self.assertEqual(len(seq_ids_to_seq['SEQ2']['seq']), 700)
            self.assertEqual(seq_ids_to_seq['SEQ2']['ss'], 'CHECHECHECCC')
            self.assertEqual(len(seq_ids_to_seq['SEQ2']['motifs_16']), 20)
            self.assertEqual(seq_ids_to_seq['SEQ3']['seq'], 'SEQ3ABCDABCD')
            self.assertEqual(seq_ids_to_seq['SEQ3']['ss'], 'CHECHECHEEEE')
            self.assertEqual(len(seq_ids_to_seq['SEQ3']['motifs_16']), 30)


    def test_get_sequences_offset_lt(self):
        with boddle(query={"page": 0, "size": 20, "offsetlt": 50}):
            result = lrr_search_web_service.get_sequences(MOTIF_VERSION)
            self.assertIsNotNone(result)
            result = json.loads(result)
            self.assertIsNotNone(result.get("sequences"), result)
            self.assertEqual(1, len(result.get("sequences")), result)
            self.assertEqual('SEQ1', result.get("sequences")[0].get("sequence_id"), result)

        with boddle(query={"page": 0, "size": 20, "offsetlt": 201}):
            result = lrr_search_web_service.get_sequences(MOTIF_VERSION)
            self.assertIsNotNone(result)
            result = json.loads(result)
            self.assertIsNotNone(result.get("sequences"), result)
            self.assertEqual(2, len(result.get("sequences")), result)
            self.assertEqual({'SEQ1','SEQ2'}, set([seq['sequence_id'] for seq in result.get("sequences")]), result)

    def test_get_sequences_offset_gt(self):
        with boddle(query={"page": 0, "size": 20, "offsetgt": 499}):
            result = lrr_search_web_service.get_sequences(MOTIF_VERSION)
            self.assertIsNotNone(result)
            result = json.loads(result)
            self.assertIsNotNone(result.get("sequences"), result)
            self.assertEqual(2, len(result.get("sequences")), result)
            self.assertEqual({'SEQ3','SEQ2'}, set([seq['sequence_id'] for seq in result.get("sequences")]), result)

        with boddle(query={"page": 0, "size": 20, "offsetgt": 600}):
            result = lrr_search_web_service.get_sequences(MOTIF_VERSION)
            self.assertIsNotNone(result)
            result = json.loads(result)
            self.assertIsNotNone(result.get("sequences"), result)
            self.assertEqual(1, len(result.get("sequences")), result)
            self.assertEqual({'SEQ3'}, set([seq['sequence_id'] for seq in result.get("sequences")]), result)


class TestFalseDiscovery(unittest.TestCase):
    def setUp(self):
        set_up_db()
        set_up_seq('SEQ2', 'SEQ2ABCDABCD', start=10, count=20)
        set_up_seq('SEQ1', 'SEQ1ABCDABCD')
        set_up_seq('SEQ3', 'SEQ3ABCDABCD', start=50, count=30)
        with dao.query_session() as session:
            seqs = dao.sequence.find_all_seqs(session)
        self.seq_ids_to_seq = dict([(seq.seq_id, seq) for seq in seqs])
        self.assertEqual(3, len(self.seq_ids_to_seq))
        self.assertEqual({'SEQ1', 'SEQ2', 'SEQ3'}, set(self.seq_ids_to_seq.keys()))
        seq_ids_to_motifs = dao.find_motifs_by_seq_ids(self.seq_ids_to_seq.keys(), MOTIF_VERSION)

    def test_tag_false_discovery(self):
        with dao.query_session() as session:
            seqs = dao.sequence.find_all_seqs(session)
            seq_ids_to_seq = dict([(seq.seq_id, seq) for seq in seqs])
        # 打上inner.falsediscovery之后，motif变为false，sequence lrr数量少1
        # 1. 打tag之前，检查SEQ1,SEQ2,SEQ3的LRR数量，及第一个motif的correct状态
        self.assertEqual(10, len(dao.find_correct_motifs_by_seq_ids(['SEQ1', ], MOTIF_VERSION)['SEQ1']))
        self.assertEqual(20, len(dao.find_correct_motifs_by_seq_ids(['SEQ2', ], MOTIF_VERSION)['SEQ2']))
        self.assertEqual(30, len(dao.find_correct_motifs_by_seq_ids(['SEQ3', ], MOTIF_VERSION)['SEQ3']))
        seq_ids_to_motifs = dao.find_motifs_by_seq_ids({'SEQ1', 'SEQ2', 'SEQ3'}, MOTIF_VERSION)
        self.assertEqual({'SEQ1', 'SEQ2', 'SEQ3'}, set(seq_ids_to_motifs.keys()))
        self.assertTrue(seq_ids_to_motifs['SEQ1'][0].correct)
        self.assertTrue(seq_ids_to_motifs['SEQ2'][0].correct)
        self.assertTrue(seq_ids_to_motifs['SEQ3'][0].correct)

        # 2. 为SEQ1、2的第一个motif打上inner.falsediscovery
        result = mark_wrong(seq_ids_to_seq['SEQ1'].id, seq_ids_to_motifs['SEQ1'][0].id)
        self.assertIsNotNone(result)
        result = json.loads(result)
        self.assertTrue(result['false_discovery'])
        result = mark_wrong(seq_ids_to_seq['SEQ2'].id, seq_ids_to_motifs['SEQ2'][0].id)
        self.assertIsNotNone(result)
        result = json.loads(result)
        self.assertTrue(result['false_discovery'])

        # 3. 重新检查SEQ1,SEQ2,SEQ3的LRR数量，及第一个motif的correct状态，SEQ1/2数量减1，correct状态为False
        self.assertEqual(9, len(dao.find_correct_motifs_by_seq_ids(['SEQ1', ], MOTIF_VERSION)['SEQ1']))
        self.assertEqual(19, len(dao.find_correct_motifs_by_seq_ids(['SEQ2', ], MOTIF_VERSION)['SEQ2']))
        self.assertEqual(30, len(dao.find_correct_motifs_by_seq_ids(['SEQ3', ], MOTIF_VERSION)['SEQ3']))
        seq_ids_to_motifs = dao.find_motifs_by_seq_ids({'SEQ1', 'SEQ2', 'SEQ3'}, MOTIF_VERSION)
        self.assertEqual({'SEQ1', 'SEQ2', 'SEQ3'}, set(seq_ids_to_motifs.keys()))
        self.assertFalse(seq_ids_to_motifs['SEQ1'][0].correct)
        self.assertFalse(seq_ids_to_motifs['SEQ2'][0].correct)
        self.assertTrue(seq_ids_to_motifs['SEQ3'][0].correct)
        self.assertEqual(10, len(seq_ids_to_motifs['SEQ1']))
        self.assertEqual(20, len(seq_ids_to_motifs['SEQ2']))
        self.assertEqual(30, len(seq_ids_to_motifs['SEQ3']))

    def test_untag_false_discovery(self):
        # 首先调用tag用例，打上tag
        self.test_tag_false_discovery()

        with dao.query_session() as session:
            seqs = dao.sequence.find_all_seqs(session)
            seq_ids_to_seq = dict([(seq.seq_id, seq) for seq in seqs])

        #  删除SEQ1中motif的tag
        seq_ids_to_motifs = dao.find_motifs_by_seq_ids({'SEQ1', 'SEQ2', 'SEQ3'}, MOTIF_VERSION)
        unmark_wrong(seq_ids_to_seq['SEQ1'].id, seq_ids_to_motifs['SEQ1'][0].id)

        # SEQ1的数量恢复
        self.assertEqual(10, len(dao.find_correct_motifs_by_seq_ids(['SEQ1', ], MOTIF_VERSION)['SEQ1']))
        self.assertEqual(19, len(dao.find_correct_motifs_by_seq_ids(['SEQ2', ], MOTIF_VERSION)['SEQ2']))
        self.assertEqual(30, len(dao.find_correct_motifs_by_seq_ids(['SEQ3', ], MOTIF_VERSION)['SEQ3']))
        seq_ids_to_motifs = dao.find_motifs_by_seq_ids({'SEQ1', 'SEQ2', 'SEQ3'}, MOTIF_VERSION)
        self.assertEqual({'SEQ1', 'SEQ2', 'SEQ3'}, set(seq_ids_to_motifs.keys()))
        self.assertTrue(seq_ids_to_motifs['SEQ1'][0].correct)
        self.assertFalse(seq_ids_to_motifs['SEQ2'][0].correct)
        self.assertTrue(seq_ids_to_motifs['SEQ3'][0].correct)

    def test_tag_multiple_times(self):
        # 首先调用tag用例，打上tag
        self.test_tag_false_discovery()

        with dao.query_session() as session:
            seq_ids_to_seq = dict([(seq.seq_id, seq) for seq in dao.sequence.find_all_seqs(session)])

        #  再打一次
        seq_ids_to_motifs = dao.find_motifs_by_seq_ids({'SEQ1', 'SEQ2', 'SEQ3'}, MOTIF_VERSION)
        mark_wrong(seq_ids_to_seq['SEQ1'].id, seq_ids_to_motifs['SEQ1'][0].id)
        mark_wrong(seq_ids_to_seq['SEQ2'].id, seq_ids_to_motifs['SEQ2'][0].id)

        # 效果与原来是一样的
        self.assertEqual(9, len(dao.find_correct_motifs_by_seq_ids(['SEQ1', ], MOTIF_VERSION)['SEQ1']))
        self.assertEqual(19, len(dao.find_correct_motifs_by_seq_ids(['SEQ2', ], MOTIF_VERSION)['SEQ2']))
        self.assertEqual(30, len(dao.find_correct_motifs_by_seq_ids(['SEQ3', ], MOTIF_VERSION)['SEQ3']))
        seq_ids_to_motifs = dao.find_motifs_by_seq_ids({'SEQ1', 'SEQ2', 'SEQ3'}, MOTIF_VERSION)
        self.assertEqual({'SEQ1', 'SEQ2', 'SEQ3'}, set(seq_ids_to_motifs.keys()))
        self.assertFalse(seq_ids_to_motifs['SEQ1'][0].correct)
        self.assertFalse(seq_ids_to_motifs['SEQ2'][0].correct)
        self.assertTrue(seq_ids_to_motifs['SEQ3'][0].correct)
        self.assertEqual(10, len(seq_ids_to_motifs['SEQ1']))
        self.assertEqual(20, len(seq_ids_to_motifs['SEQ2']))
        self.assertEqual(30, len(seq_ids_to_motifs['SEQ3']))


class TestManuallyMotif(unittest.TestCase):
    def setUp(self):
        set_up_db()
        set_up_seq('SEQ2', random_seq(700), start=10, count=20)
        set_up_seq('SEQ1', random_seq(300))
        set_up_seq('SEQ3', random_seq(1700), start=50, count=30)

    def test_add_manually_motif_invalid_input(self):
        # sid不是数字
        with boddle(json={"offset": 10}):
            result = lrr_search_web_service.add_manually_motif(MOTIF_VERSION, "a")
        self.assertIsNotNone(result)
        result = json.loads(result)
        self.assertDictEqual({'message': str.format(ErrorCode.INVALID_PARA, "sid", 'a')}, result)

        # sid不存在
        with boddle(json={"offset": 10}):
            result = lrr_search_web_service.add_manually_motif(MOTIF_VERSION, 10240)
        self.assertIsNotNone(result)
        result = json.loads(result)
        self.assertDictEqual({'message': str.format(ErrorCode.OBJECT_NOT_EXISTS, 10240)}, result)

        # offset不是数字
        with dao.query_session() as session:
            seqs = dao.sequence.find_all_seqs(session)
            seq_ids_to_seq = dict([(seq.seq_id, seq) for seq in seqs])
        seq1_sid = seq_ids_to_seq['SEQ1'].id
        seq2_sid = seq_ids_to_seq['SEQ2'].id
        with boddle(json={"offset": "b"}):
            result = lrr_search_web_service.add_manually_motif(MOTIF_VERSION, seq1_sid)
        self.assertIsNotNone(result)
        result = json.loads(result)
        self.assertDictEqual({'message': str.format(ErrorCode.INVALID_PARA, "offset", 'b')}, result)

        # offset出现重叠
        with boddle(json={"offset": 195}):
            result = lrr_search_web_service.add_manually_motif(MOTIF_VERSION, seq1_sid)
        self.assertIsNotNone(result)
        result = json.loads(result)
        self.assertDictEqual({'message': str.format(ErrorCode.OFFSET_OVERLAP)}, result)

        # offset出现重叠
        with boddle(json={"offset": 185}):
            result = lrr_search_web_service.add_manually_motif(MOTIF_VERSION, seq2_sid)
        self.assertIsNotNone(result)
        result = json.loads(result)
        self.assertDictEqual({'message': str.format(ErrorCode.OFFSET_OVERLAP)}, result)

        # 不可以在已标记为错误的offset上新增，因为直接取消错误标记就可以了
        seq_ids_to_motifs = dao.find_motifs_by_seq_ids(['SEQ1', ], MOTIF_VERSION)
        seq1_second_motif = seq_ids_to_motifs['SEQ1'][1]
        mark_wrong(seq1_sid, seq1_second_motif.id)
        with boddle(json={"offset": seq1_second_motif.offset}):
            result = lrr_search_web_service.add_manually_motif(MOTIF_VERSION, seq1_sid)
        self.assertIsNotNone(result)
        result = json.loads(result)
        self.assertDictEqual({'message': ErrorCode.OFFSET_EXISTS_WRONG}, result)

    def test_add_manually_motif_ok(self):
        set_up_baseline_seq('SEQ4', '', start=10, count=10, step=24)
        # 先把seq都查出来
        with dao.query_session() as session:
            seqs = dao.sequence.find_all_seqs(session)
            seq_ids_to_seq = dict([(seq.seq_id, seq) for seq in seqs])
        seq1_sid = seq_ids_to_seq['SEQ1'].id

        # 正常添加前，motif有10个
        with boddle(query={"page": 0, "size": 20}):
            result = lrr_search_web_service.get_sequences(MOTIF_VERSION)
        self.assertIsNotNone(result)
        result = json.loads(result)
        self.assertIsNotNone(result.get("sequences"))
        seq_ids_to_seq = dict([(seq.get("sequence_id"), seq) for seq in result.get("sequences")])
        self.assertEqual(10, len(seq_ids_to_seq['SEQ1'].get("motifs_16")))

        # 正常添加OK
        with boddle(json={"offset": 196}):
            result = lrr_search_web_service.add_manually_motif(MOTIF_VERSION, seq1_sid)
        self.assertIsNotNone(result)
        result = json.loads(result)
        self.assertTrue("motifs_16" in result)
        self.assertEqual(196, result['motifs_16'][0]['offset'])

        # 查询sequance，motifs数加1
        with boddle(query={"page": 0, "size": 20}):
            result = lrr_search_web_service.get_sequences(MOTIF_VERSION)
        self.assertIsNotNone(result)
        result = json.loads(result)
        self.assertIsNotNone(result.get("sequences"))
        seq_ids_to_seq = dict([(seq.get("sequence_id"), seq) for seq in result.get("sequences")])
        self.assertEqual(11, len(seq_ids_to_seq['SEQ1'].get("motifs_16")), seq_ids_to_seq['SEQ1'])

    def test_add_manually_motif_before_wrong_area(self):
        # 先把seq都查出来
        with dao.query_session() as session:
            seqs = dao.sequence.find_all_seqs(session)
            seq_ids_to_seq = dict([(seq.seq_id, seq) for seq in seqs])
        sid = seq_ids_to_seq['SEQ2'].id
        ms = dao.find_motifs_by_seq_ids(['SEQ2',], MOTIF_VERSION)['SEQ2']
        ms.sort(key=lambda m:m.offset)
        motif = ms[0]

        # 重叠(新增的offset在前面面)时无法添加
        with boddle(json={"offset": 185}):
            result = lrr_search_web_service.add_manually_motif(MOTIF_VERSION, sid)
        self.assertIsNotNone(result)
        result = json.loads(result)
        self.assertDictEqual({'message': str.format(ErrorCode.OFFSET_OVERLAP)}, result)

        # 标记错误
        mark_wrong(sid, motif.id)

        # 再次添加OK
        with boddle(json={"offset": 185}):
            result = lrr_search_web_service.add_manually_motif(MOTIF_VERSION, sid)
        self.assertIsNotNone(result)
        result = json.loads(result)
        self.assertTrue("motifs_16" in result)
        self.assertEqual(185, result['motifs_16'][0]['offset'])

    def test_add_manually_motif_after_wrong_area(self):
        set_up_baseline_seq('SEQ4', '', start=10, count=10, step=24)
        # 先把seq都查出来
        with dao.query_session() as session:
            seqs = dao.sequence.find_all_seqs(session)
            seq_ids_to_seq = dict([(seq.seq_id, seq) for seq in seqs])
        sid = seq_ids_to_seq['SEQ2'].id
        ms = dao.find_motifs_by_seq_ids(['SEQ2',], MOTIF_VERSION)['SEQ2']
        ms.sort(key=lambda m:m.offset)
        motif = ms[-1]

        # 重叠(新增的offset在前面面)时无法添加
        with boddle(json={"offset": 595}):
            result = lrr_search_web_service.add_manually_motif(MOTIF_VERSION, sid)
        self.assertIsNotNone(result)
        result = json.loads(result)
        self.assertDictEqual({'message': str.format(ErrorCode.OFFSET_OVERLAP)}, result)

        # 标记错误
        mark_wrong(sid, motif.id)

        # 再次添加OK
        with boddle(json={"offset": 595}):
            result = lrr_search_web_service.add_manually_motif(MOTIF_VERSION, sid)
        self.assertIsNotNone(result)
        result = json.loads(result)
        self.assertEqual(1, len(result.get('motifs_16', [])))
        self.assertEqual(595, result['motifs_16'][0].get('offset', -1))

    def _get_manually_motifs(self, seq_id):
        seq_ids_to_motifs = dao.find_motifs_by_seq_ids([seq_id, ], MOTIF_VERSION)
        self.assertTrue(seq_id in seq_ids_to_motifs)
        manually_motifs = []
        for m in seq_ids_to_motifs[seq_id]:
            if m.manually_add:
                manually_motifs.append(m)
        self.assertTrue(len(manually_motifs) > 0)
        return manually_motifs

    def test_add_manually_motif_can_not_mark_wrong(self):
        # 先把seq都查出来
        with dao.query_session() as session:
            seqs = dao.sequence.find_all_seqs(session)
            seq_ids_to_seq = dict([(seq.seq_id, seq) for seq in seqs])
        seq1_sid = seq_ids_to_seq['SEQ1'].id

        # 正常添加OK
        with boddle(json={"offset": 196}):
            result = lrr_search_web_service.add_manually_motif(MOTIF_VERSION, seq1_sid)
        self.assertIsNotNone(result)
        result = json.loads(result)
        self.assertTrue('motifs_16' in result)
        self.assertEqual(196, result['motifs_16'][0]['offset'])
        manually_motif = self._get_manually_motifs('SEQ1')[0]

        # 不允许对手动添加的motif打tag
        result = mark_wrong(seq1_sid, manually_motif.id)
        result = json.loads(result)
        self.assertTrue(result.get("message").endswith("and should not be tagged wrong again. The motif can be delete directly"))

    def test_can_not_unmark_wrong_overlapped_with_manually_added_one(self):
        # 如果手工添加的motif与原有已经mark为wrong的moiti存在覆盖，那么手动添加后，原有的motif就不可以unmark了
        # 细分的话，又分两种情况，手动motif的offset在wrong的后面，与在其前面
        # 首先对SEQ2的第一个和最后一个motif标记wrong，然后手工添加motif
        self.test_add_manually_motif_after_wrong_area()
        self.test_add_manually_motif_before_wrong_area()

        # 尝试对第一个和最后一个motif解除wrong标记，应该都要失败（因为前面添加过两个手动motif了，所以最后一个应该是-3）
        seq_ids_to_motifs = dao.find_motifs_by_seq_ids(['SEQ2',], MOTIF_VERSION)
        seq_ids_to_motifs['SEQ2'].sort(key=lambda m:m.offset)
        first_motif = seq_ids_to_motifs['SEQ2'][1]
        last_motif = seq_ids_to_motifs['SEQ2'][-2]

        with dao.query_session() as session:
            seqs = dao.sequence.find_all_seqs(session)
            seq_ids_to_seq = dict([(seq.seq_id, seq) for seq in seqs])

        result = unmark_wrong(seq_ids_to_seq['SEQ2'].id, first_motif.id)
        self.assertIsNotNone(result)
        msg = result.get('message', '')
        self.assertTrue(msg.endswith(' can not be unmark from wrong because overlapping was found'), result)

        result = unmark_wrong(seq_ids_to_seq['SEQ2'].id, last_motif.id)
        self.assertIsNotNone(result)
        msg = result.get('message', '')
        self.assertTrue(msg.endswith(' can not be unmark from wrong because overlapping was found'), result)

        # 删除手工添加的motif
        for m in self._get_manually_motifs('SEQ2'):
            delete_motif(m.id)

        # 尝试对第一个和最后一个解除wrong标记，成功
        result = unmark_wrong(seq_ids_to_seq['SEQ2'].id, first_motif.id)
        self.assertFalse(result['false_discovery'])
        result = unmark_wrong(seq_ids_to_seq['SEQ2'].id, last_motif.id)
        self.assertFalse(result['false_discovery'])

    def test_unmark_wrong_with_other_wrong(self):
        # unmark一个wrong motif的时候，考虑已存在的mark为wrong的motif，不应该影响到正常的计算
        set_up_seq('SEQ4', 'SEQ4ABCDABCDABCD', start=10, count=2, step=10)
        seq4 = dao.find_seq_by_id('SEQ4')
        motifs = dao.find_correct_motifs_by_seq_ids(['SEQ4',], MOTIF_VERSION)['SEQ4']
        self.assertEqual(2, len(motifs))
        mark_wrong(seq4.id, motifs[1].id)
        mark_wrong(seq4.id, motifs[0].id)
        self.assertEqual(0, len(dao.find_correct_motifs_by_seq_ids(['SEQ4',], MOTIF_VERSION)['SEQ4']))
        result = unmark_wrong(seq4.id, motifs[0].id)
        self.assertEqual(1, len(dao.find_correct_motifs_by_seq_ids(['SEQ4', ], MOTIF_VERSION)['SEQ4']))

    def test_delete_manually_motif_ok(self):
        # 先把seq都查出来
        with dao.query_session() as session:
            seqs = dao.sequence.find_all_seqs(session)
            seq_ids_to_seq = dict([(seq.seq_id, seq) for seq in seqs])
        seq1_sid = seq_ids_to_seq['SEQ1'].id

        # 正常添加OK
        with boddle(json={"offset": 196}):
            result = lrr_search_web_service.add_manually_motif(MOTIF_VERSION, seq1_sid)
        self.assertIsNotNone(result)
        result = json.loads(result)
        self.assertTrue("motifs_16" in result)
        self.assertEqual(196, result['motifs_16'][0]['offset'])
        manually_motif = self._get_manually_motifs('SEQ1')[0]

        # 删除OK
        result = delete_motif(manually_motif.id)
        self.assertIsNotNone(result)
        self.assertDictEqual({}, json.loads(result))

        # 重复删除，就失败了，当前暂不考虑幂等
        result = delete_motif(manually_motif.id)
        self.assertIsNotNone(result)

    def test_delete_auto_motif_failed(self):
        seq_ids_to_motifs = dao.find_motifs_by_seq_ids(['SEQ1',], MOTIF_VERSION)
        result = delete_motif(seq_ids_to_motifs['SEQ1'][0].id)
        self.assertIsNotNone(result)

        seq_ids_to_motifs = dao.find_motifs_by_seq_ids(['SEQ1',], MOTIF_VERSION)
        result = delete_motif(seq_ids_to_motifs['SEQ1'][4].id)
        self.assertIsNotNone(result)


if __name__ == '__main__':
    unittest.main()
