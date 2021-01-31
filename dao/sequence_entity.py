# THIS FILE IS PART OF phytolrr.com PROJECT.
# Copyright 2019-2021 phytolrr.com. All rights reserved.

from sqlalchemy import Column, Integer, String, Text, Boolean
from sqlalchemy import update
from dao.datasource import *


class SequenceEntity(Base):
    __tablename__ = 'sequenceentity'

    id = Column(Integer, autoincrement=True, primary_key=True)
    seq_id = Column(String(length=256), index=True, unique=True)
    seq = Column(Text)
    ss = Column(Text, nullable=True)
    ss8 = Column(Text, nullable=True)
    acc = Column(Text, nullable=True)
    acc20 = Column(Text, nullable=True)
    species = Column(String(length=8), index=True)
    subgroup = Column(String(length=16), index=True)
    baseline = Column(Boolean, nullable=True)

    def __repr__(self):
        return str.format("<{}: {{id={}, seq_id={}}}>",
                          SequenceEntity.__class__.__name__,
                          self.id, self.seq_id)


def add_seq(session, seq_id, seq, ss=None, ss8=None, acc=None, acc20=None, species=None):
    seq_entity = SequenceEntity()
    seq_entity.seq_id = seq_id
    seq_entity.seq = seq
    seq_entity.ss = ss
    seq_entity.ss8 = ss8
    seq_entity.acc = acc
    seq_entity.acc20 = acc20
    seq_entity.species = species
    session.add(seq_entity)


def seq_exists(session, seq_id):
    return session.query(SequenceEntity.id).filter(SequenceEntity.seq_id == seq_id).scalar() is not None


def seq_exists_by_sid(session, sid):
    return session.query(SequenceEntity.id).filter(SequenceEntity.id == sid).scalar() is not None


def update_seq_for_scratch_result(session, seq_id, ss, ss8, acc, acc20):
    seq = find_seq_by_id(session, seq_id)
    seq.ss = ss
    seq.ss8 = ss8
    seq.acc = acc
    seq.acc20 = acc20
    session.add(seq)


def find_seq_by_id(session, seq_id):
    return session.query(SequenceEntity).filter(SequenceEntity.seq_id == seq_id).one()


def find_seq_by_ids(session, seq_ids):
    return session.query(SequenceEntity).filter(SequenceEntity.seq_id.in_(seq_ids)).all()


def find_seq_by_sid(session, sid):
    return session.query(SequenceEntity).filter(SequenceEntity.id == sid).one()


def find_all_seqs(session):
    return session.query(SequenceEntity).all()


def find_all_seq_ids(session):
    return [ret[0] for ret in session.query(SequenceEntity.seq_id).all()]


def find_seqs_by_baseline(session):
    return session.query(SequenceEntity).filter(SequenceEntity.baseline).all()


def update_subgroup_by_seq_ids(session, seq_ids, subgroup):
    sql = update(SequenceEntity).where(SequenceEntity.seq_id.in_(seq_ids)).values(subgroup=subgroup)
    session.execute(sql)


def find_seqs_by_subgroups(session, subgroups):
    return session.query(SequenceEntity).filter(SequenceEntity.subgroup.in_(subgroups)).all()
