# -*- coding: utf-8 -*
# THIS FILE IS PART OF phytolrr.com PROJECT.
# Copyright 2019-2021 phytolrr.com. All rights reserved.

import logging

from dao.datasource import *
import dao.sequence_entity as sequence
import dao.motif_entity as motif
import dao.nsite_entity as nsite
import dao.tag_ref_entity as tag
from sqlalchemy import delete, select, update, exists, func
from tools.exception import ValidationError


def reconnect():
    datasource.reconnect()


def find_tags_by_motif_ids(motif_ids):
    with query_session() as session:
        return tag.find_tag_names_by_motif_ids(session, motif_ids)


def find_seq_ids_by_motif_ids(motif_ids, version):
    with query_session() as session:
        return motif.find_seq_ids_by_motif_ids(session, motif_ids, version)


def _update_motif_correct(session, ids, version):
    for mid in ids:
        wrong = session.query(exists().
                              where(tag.TagRefEntity.target_id == mid).
                              where(tag.TagRefEntity.tag_name.in_(tag.WRONG_MOTIF_TAGS))).scalar()
        cls = motif.get_entity(version)
        stmt = update(cls).where(cls.id == mid).values(correct=not wrong)
        session.execute(stmt)


def replace_tags_by_motifs(ids_to_tags, version):
    with session_scope() as session:
        tag.replace_tags_by_motifs(session, ids_to_tags)
        _update_motif_correct(session, ids_to_tags.keys(), version)


def add_tags_by_names_to_ids(names_to_motif_ids, version):
    ids = set([mid for ids in names_to_motif_ids.values() for mid in ids])
    with session_scope() as session:
        tag.add_tags_by_names_to_motifs(session, names_to_motif_ids)
        _update_motif_correct(session, ids, version)


def add_manually_motif(seq_id, offset, version, score, probability):
    motif_e = motif.MotifEntityBase(offset, seq_id, score, probability, 0.0, manually_add=True)
    with session_scope() as session:
        motif.add_motifs(session, [motif_e], version)
    return None


def remove_manually_motif(mid, version):
    with session_scope() as session:
        tag.replace_tags_by_motifs(session, {mid: []})
        deleted_count = motif.delete_manually_motifs_by_ids(session, [mid], version, synchronize_session='fetch')
        if deleted_count != 1:
            raise ValidationError(str.format("Failed to delete motif by id {}", mid))


def seq_exists(seq_id):
    with query_session() as session:
        return sequence.seq_exists(session, seq_id)


def seq_exists_by_sid(sid):
    with query_session() as session:
        return sequence.seq_exists_by_sid(session, sid)


def motif_exists(id, version):
    with query_session() as session:
        return motif.motif_exists(session, id, version)


def add_seq(seq_id, seq_str, version, motifs_16=[], species=None, ss=None):
    with session_scope() as session:
        motif.replace_motifs_by_seq(session, seq_id, motifs_16, version)
        sequence.add_seq(session, seq_id, seq_str, species=species, ss=ss)


def remove_tags_by_seq_id(session, seq_id, version):
    cls = motif.get_entity(version)
    stmt = delete(tag.TagRefEntity).where(tag.TagRefEntity.target_id.in_(
        select(cls.id).select_from(cls)\
        .where(cls.seq_id == seq_id)
    ))
    session.execute(stmt)


def update_seq_for_scratch_result(seq_id, ss, ss8, acc, acc20):
    with session_scope() as session:
        sequence.update_seq_for_scratch_result(session, seq_id, ss, ss8, acc, acc20)


def find_seq_by_id(seq_id):
    with query_session() as session:
        return sequence.find_seq_by_id(session, seq_id)


def find_seq_by_sid(sid):
    with query_session() as session:
        return sequence.find_seq_by_sid(session, sid)


def replace_nsites_by_seq_id(seq_id, positions, ntypes):
    with session_scope() as session:
        nsite.replace_nsites_by_seq_id(session, seq_id, positions, ntypes)


class Sequence(object):
    def __init__(self, id, seq_id, seq):
        self.id = id
        self.seq_id = seq_id
        self.seq = seq


def _query_sequences_need_join_motif(filters):
    motif_filters = {'offsetlt', 'offsetgt', 'offseteq', 'lrr_count_lt', 'lrr_count_gt', 'lrr_count_eq'}
    return len(motif_filters & set(filters.keys())) > 0


def _add_offset_filters(filters, cls, query):
    offset_lt = filters.get('offsetlt', None)
    if offset_lt is not None:
        query = query.filter(cls.offset < offset_lt)

    offset_gt = filters.get('offsetgt', None)
    if offset_gt is not None:
        query = query.filter(cls.offset > offset_gt)

    offset_eq = filters.get('offseteq', None)
    if offset_eq is not None:
        query = query.filter(cls.offset == offset_eq)
    return query


def _add_keyword_filters(filters, query):
    keyword = filters.get('keyword', None)
    if keyword is not None:
        query = query.filter(sequence.SequenceEntity.seq_id.like('%' + keyword + '%'))
    return query


def _add_species_filter(filters, query):
    species = filters.get('species', None)
    if species is not None:
        query = query.filter(sequence.SequenceEntity.species.in_(species))
    return query


def _add_lrr_count_havings(filters, cls, query):
    lrr_count_lt = filters.get('lrr_count_lt', None)
    if lrr_count_lt is not None:
        query = query.having(func.count(cls.offset) < lrr_count_lt)

    lrr_count_gt = filters.get('lrr_count_gt', None)
    if lrr_count_gt is not None:
        query = query.having(func.count(cls.offset) > lrr_count_gt)

    lrr_count_eq = filters.get('lrr_count_eq', None)
    if lrr_count_eq is not None:
        query = query.having(func.count(cls.offset) == lrr_count_eq)
    return query


def query_sequences(filters, version):
    '''
    query sequences according to filters
    :param filters: dictionary, KEY acceptable:
    page_size, page_index: required, return page size and page number (start from 0)
    offsetlt, offsetgt, offseteq: the offset less than, greater than and equal than
    lrr_count_lt, lrr_count_gr, lrr_count_eq: the number of lrr less than, greater than and equal than
    keyword: query the sequences by ID, case insensitive
    species: species
    :return: return seqences list, ordered by ID
    '''
    cls = motif.get_entity(version)
    with query_session() as session:
        query = session.query(sequence.SequenceEntity)
        count_query = session.query(sequence.SequenceEntity.id)
        if _query_sequences_need_join_motif(filters):
            query = query.outerjoin(cls, sequence.SequenceEntity.seq_id == cls.seq_id).filter(cls.correct)
            count_query = count_query.outerjoin(cls, sequence.SequenceEntity.seq_id == cls.seq_id).filter(cls.correct)

            query = _add_offset_filters(filters, cls, query)
            count_query = _add_offset_filters(filters, cls, count_query)

        query = _add_keyword_filters(filters, query)
        count_query = _add_keyword_filters(filters, count_query)

        query = _add_species_filter(filters, query)
        count_query = _add_species_filter(filters, count_query)

        if _query_sequences_need_join_motif(filters):
            query = query.group_by(sequence.SequenceEntity.seq_id)
            count_query = count_query.group_by(sequence.SequenceEntity.seq_id)

            query = _add_lrr_count_havings(filters, cls, query)
            count_query = _add_lrr_count_havings(filters, cls, count_query)

        query = query.order_by(sequence.SequenceEntity.seq_id).limit(filters['page_size']).offset(filters['page_index'])

        return query.all(), len(count_query.all())


def find_motifs_by_seq_ids(seq_ids, version, with_wrong=True):
    with query_session() as session:
        return motif.find_motifs_by_seq_ids(session, seq_ids, version, with_wrong=with_wrong)


def find_motifs_by_offsets(seq_id, offsets, version):
    with query_session() as session:
        return motif.find_motifs_by_offsets(session, seq_id, offsets, version)


def find_correct_motifs_by_seq_ids(seq_ids, version):
    with query_session() as session:
        return motif.find_motifs_by_seq_ids(session, seq_ids, version, False)


def find_baseline_seqs():
    with query_session() as session:
        return sequence.find_seqs_by_baseline(session)


def find_baseline_motifs(baseline_version=1, with_wrong=True):
    seqs = find_baseline_seqs()
    seq_ids_to_seq_str = dict([(seq.seq_id, seq.seq) for seq in seqs])
    logging.info(str.format("Baseline sequence ids count({}): {}", len(seq_ids_to_seq_str), seq_ids_to_seq_str.keys()))

    seq_ids_to_motifs = find_motifs_by_seq_ids(seq_ids_to_seq_str.keys(), baseline_version, with_wrong)
    return [seq_ids_to_seq_str[m.seq_id][m.offset:m.offset+16] for motifs in seq_ids_to_motifs.values() for m in motifs]


def update_false_discovery_by_motif(mid, false_discovery, version):
    if version < 2:
        raise ValidationError(str.format("The version below 2 is no longer support to mark wrong"))
    with session_scope() as session:
        motif.update_false_discovery_by_motif(session, mid, false_discovery, version)


def find_motif_by_mid(mid, version):
    with query_session() as session:
        ms = motif.find_motifs_by_ids(session, [mid], version)
    if ms is None or len(ms) == 0:
        return None
    return ms[0]


def find_seqs_by_subgroups(subgroups):
    with query_session() as session:
        return sequence.find_seqs_by_subgroups(session, subgroups)


def find_nsites_by_seq_id(seq_id):
    with query_session() as session:
        return nsite.find_nsites_by_seq_id(session, seq_id)


def find_nsites_by_seq_ids(seq_ids):
    with query_session() as session:
        return nsite.find_nsites_by_seq_ids(session, seq_ids)
