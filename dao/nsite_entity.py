# THIS FILE IS PART OF phytolrr.com PROJECT.
# Copyright 2019-2021 phytolrr.com. All rights reserved.

from sqlalchemy import Column, String, Integer
from sqlalchemy import delete
from dao.datasource import *

# The type of a nsite
S = 0
T = 1


class NSiteEntity(Base):
    __tablename__ = 'nsiteentity'

    id = Column(Integer, autoincrement=True, primary_key=True)
    seq_id = Column(String(length=256), index=True)
    start_pos = Column(Integer)
    ntype = Column(Integer, index=True)

    def __init__(self, seq_id, start_pos, ntype):
        self.seq_id = seq_id
        self.start_pos = start_pos
        self.ntype = ntype


def remove_nsites_by_ids(session, ids):
    stmt = delete(NSiteEntity).where(NSiteEntity.id.in_(ids))
    session.execute(stmt)


def find_nsites_by_seq_id(session, seq_id):
    return session.query(NSiteEntity).filter(NSiteEntity.seq_id == seq_id).all()


def find_nsites_by_seq_ids(session, seq_ids):
    nsites = session.query(NSiteEntity).filter(NSiteEntity.seq_id.in_(seq_ids)).all()
    seq_ids_to_nsites = dict([(seq_id, []) for seq_id in seq_ids])
    for nsite in nsites:
        seq_ids_to_nsites[nsite.seq_id].append(nsite)
    return seq_ids_to_nsites


def find_nsites_count_by_seq_id(session, seq_id):
    return session.query(NSiteEntity.id).filter(NSiteEntity.seq_id == seq_id).count()


def replace_nsites_by_seq_id(session, seq_id, positions, ntypes):
    assert len(positions) == len(ntypes)
    nsites = find_nsites_by_seq_id(session, seq_id)
    old_poses_to_entity = dict([nsite.start_pos, nsite] for nsite in nsites)

    # Add or update
    nsites_to_add = []
    nsites_to_update = []
    for pos, ntype in zip(positions, ntypes):
        old_entity = old_poses_to_entity.get(pos, None)
        if old_entity is None:
            nsites_to_add.append(NSiteEntity(seq_id, pos, ntype))
        else:
            if old_entity.ntype != ntype:
                old_entity.ntype = ntype
                nsites_to_update.append(old_entity)

    if len(nsites_to_add) > 0:
        session.add_all(nsites_to_add)
    if len(nsites_to_update) > 0:
        session.add_all(nsites_to_update)

    # Delete
    poses_to_del = set(old_poses_to_entity.keys()) - set(positions)
    ids_to_del = set([old_poses_to_entity[pos].id for pos in poses_to_del])
    if len(ids_to_del) > 0:
        remove_nsites_by_ids(session, ids_to_del)
