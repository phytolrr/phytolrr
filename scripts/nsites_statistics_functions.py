# -*- coding: utf-8 -*
import dao


class GeneralFeaturesTableRow(object):
    def __init__(self):
        self.seq_ids = []
        self.lrr_count = None
        self.nsites_count = None
        self.seq_ids_to_lrr_count = {}
        self.seq_ids_to_nsites_count = {}


def get_num_by_seq_ids(seq_ids):
    '''
    :param seq_ids: 序列ID列表
    :return: GeneralFeaturesTableRow
    '''
    row = GeneralFeaturesTableRow()
    row.seq_ids = list(seq_ids)
    row.seq_ids_to_lrr_count = dict([(seq_id, 0) for seq_id in seq_ids])
    row.seq_ids_to_nsites_count = dict([(seq_id, 0) for seq_id in seq_ids])

    return row
