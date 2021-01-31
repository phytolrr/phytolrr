# THIS FILE IS PART OF phytolrr.com PROJECT.
# Copyright 2019-2021 phytolrr.com. All rights reserved.

from sqlalchemy import Column, Integer, String, Boolean, REAL
from sqlalchemy import delete
from dao.datasource import *


class MotifEntityBase(object):
    id = Column(Integer, autoincrement=True, primary_key=True)
    seq_id = Column(String(256), index=True)
    offset = Column(Integer)
    score = Column(REAL)
    probability = Column(REAL)
    fdr_probability = Column(REAL)
    false_discovery = Column(Boolean, default=False)
    manually_add = Column(Boolean, default=False)
    correct = Column(Boolean, default=True, index=True)

    def __init__(self, offset, seq_id, score, probability, fdr_probability, manually_add=False, correct=True):
        self.offset = offset
        self.seq_id = seq_id
        self.score = score
        self.correct = correct
        self.manually_add = manually_add
        self.probability = probability
        self.fdr_probability = fdr_probability

    def __repr__(self):
        return str.format("<{} at {}> seq_id {}, offset {}, score {}, false_discovery {}, manually_add {} , correct {}",
                          self.__class__, id(self),
                          self.seq_id, self.offset, self.score,
                          self.false_discovery, self.manually_add, self.correct)


class MotifEntityV1(MotifEntityBase, Base):
    __tablename__ = 'motifentity'

    def __init__(self):
        pass


class MotifEntityV2(MotifEntityBase, Base):
    __tablename__ = 'motifentity_v2'

    def __init__(self):
        pass


class MotifEntityV3(MotifEntityBase, Base):
    __tablename__ = 'motifentity_v3'

    def __init__(self):
        pass


VERSIONS_TO_ENTITY = {
    1: MotifEntityV1,
    '1': MotifEntityV1,
    2: MotifEntityV2,
    '2': MotifEntityV2,
    3: MotifEntityV3,
    '3': MotifEntityV3
}


def get_entity(version):
    return VERSIONS_TO_ENTITY.get(version)


def convert_to_entity(inputs, version):
    cls = get_entity(version)
    entities = []
    for input in inputs:
        if not isinstance(input, MotifEntityBase):
            raise ValueError("Invalid input entity type " + type(input))
        entity = cls()
        for column in cls.__table__.columns:
            value = getattr(input, column.key)
            if type(value) == Column:
                continue
            setattr(entity, column.key, value)
        entities.append(entity)
    return entities


def add_motifs(session, motif_entities, version):
    session.add_all(convert_to_entity(motif_entities, version))


def delete_motifs_by_seq_ids(session, seq_ids, version):
    cls = get_entity(version)
    stmt = delete(cls).where(cls.seq_id.in_(seq_ids))
    session.execute(stmt)


def find_motifs_by_seq_ids(session, seq_ids, version, with_wrong=True):
    cls = get_entity(version)
    if with_wrong:
        motifs = session.query(cls).filter(cls.seq_id.in_(seq_ids)).all()
    else:
        motifs = session.query(cls).filter(cls.seq_id.in_(seq_ids)).filter(cls.correct).all()
    seq_ids_to_motifs = dict([(seq_id, []) for seq_id in seq_ids])
    for m in motifs:
        seq_ids_to_motifs[m.seq_id].append(m)
    return seq_ids_to_motifs


def find_seq_ids_by_motif_ids(session, motif_ids, version):
    cls = get_entity(version)
    seq_ids = session.query(cls.seq_id).distinct().filter(cls.id.in_(motif_ids)).all()
    return [seq_id for seq_id_tuple in seq_ids for seq_id in seq_id_tuple]


def find_motifs_by_seq_id(session, seq_id, version, with_wrong=True):
    cls = get_entity(version)
    if with_wrong:
        return session.query(cls)\
            .filter(cls.seq_id == seq_id)\
            .order_by(cls.offset).all()
    else:
        return session.query(cls)\
            .filter(cls.seq_id == seq_id)\
            .filter(cls.correct)\
            .order_by(cls.offset).all()


def motif_exists(session, id, version):
    cls = get_entity(version)
    return session.query(cls.id).filter(cls.id == id).scalar() is not None


def find_motifs_count_by_seq_id(session, seq_id, version):
    cls = get_entity(version)
    return session.query(cls.id).filter(cls.seq_id == seq_id).filter(cls.correct).count()


def replace_motifs_by_seq(session, seq_id, motifs, version):
    motifs = convert_to_entity(motifs, version)
    old_motifs = find_motifs_by_seq_id(session, seq_id, version)
    old_offsets_to_motif = dict([(m.offset, m) for m in old_motifs])
    new_offsets_to_motif = dict([(m.offset, m) for m in motifs])
    offsets_to_del = set(old_offsets_to_motif.keys()) - set(new_offsets_to_motif.keys())
    ids_to_del = set([old_offsets_to_motif[offset].id for offset in offsets_to_del])
    motifs_to_update = []
    for offset, new_motif in new_offsets_to_motif.items():
        motif = old_offsets_to_motif.get(offset, None)
        if motif is None:
            motifs_to_update.append(new_motif)
        else:
            motif.score = new_motif.score
            motif.probability = new_motif.probability
            motif.fdr_probability = new_motif.fdr_probability
            motif.correct = new_motif.correct
            motifs_to_update.append(motif)

    if len(ids_to_del) > 0:
        cls = get_entity(version)
        session.execute(delete(cls).where(cls.id.in_(ids_to_del)))
    if len(motifs_to_update) > 0:
        session.add_all(motifs_to_update)


def find_motifs_by_ids(session, motif_ids, version):
    cls = get_entity(version)
    return session.query(cls).filter(cls.id.in_(motif_ids)).all()


def delete_manually_motifs_by_ids(session, motif_ids, version, synchronize_session=False):
    cls = get_entity(version)
    return session.query(cls)\
        .filter(cls.id.in_(motif_ids)).filter(cls.manually_add)\
        .delete(synchronize_session=synchronize_session)


def find_motifs_by_offsets(session, seq_id, offsets, version):
    cls = get_entity(version)
    return session.query(cls).filter(cls.seq_id == seq_id).filter(cls.offset.in_(offsets)).all()


def update_false_discovery_by_motif(session, mid, false_discovery, version):
    cls = get_entity(version)
    session.query(cls).filter(cls.id == mid).update({"false_discovery": false_discovery, "correct": not false_discovery})


def update_score_by_motif(session, mid, score, version):
    cls = get_entity(version)
    session.query(cls).filter(cls.id == mid).update({"score": score})
