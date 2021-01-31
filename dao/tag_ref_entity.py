# THIS FILE IS PART OF phytolrr.com PROJECT.
# Copyright 2019-2021 phytolrr.com. All rights reserved.

from sqlalchemy import Column, String, Integer
from sqlalchemy import delete
from dao.datasource import *

OVERLAP_TAG = 'inner.overlap'


WRONG_MOTIF_TAGS = {OVERLAP_TAG}


class TagRefEntity(Base):
    __tablename__ = "tagrefentity"

    id = Column(Integer, autoincrement=True, primary_key=True)
    tag_name = Column(String(256), index=True)
    target_id = Column(Integer, index=True)

    def __init__(self, target_id, tag_name):
        self.tag_name = tag_name
        self.target_id = target_id


def find_tags_by_motif_ids(session, motif_ids):
    tag_entities = session.query(TagRefEntity).filter(TagRefEntity.target_id.in_(motif_ids)).all()
    motif_ids_to_tags = {}
    for tag_ref in tag_entities:
        if tag_ref.target_id in motif_ids_to_tags:
            motif_ids_to_tags[tag_ref.target_id].append(tag_ref)
        else:
            motif_ids_to_tags[tag_ref.target_id] = [tag_ref,]
    return motif_ids_to_tags


def tag_ref_exists(session, tag_name, motif_id):
    return session.query(TagRefEntity.id)\
        .filter(TagRefEntity.tag_name == tag_name)\
        .filter(TagRefEntity.target_id == motif_id).scalar() is not None


def find_tag_refs_by_tag_name_and_motif_ids(session, tag_name, motif_ids):
    return session.query(TagRefEntity)\
        .filter(TagRefEntity.tag_name == tag_name)\
        .filter(TagRefEntity.target_id.in_(motif_ids)).all()


def add_tags_by_map(session, target_ids_to_tags):
    tag_refs = []
    for target_id, tag_names in target_ids_to_tags.items():
        for tag_name in tag_names:
            tag_refs.append(TagRefEntity(target_id, tag_name))
    session.add_all(tag_refs)


def remove_by_tag_ids(session, ids):
    stmt = delete(TagRefEntity).where(TagRefEntity.id.in_(ids))
    session.execute(stmt)


def find_tag_names_by_motif_ids(session, motif_ids):
    mids_to_tags = find_tags_by_motif_ids(session, motif_ids)
    mids_to_tag_names = {}
    for mid, tags in mids_to_tags.items():
        mids_to_tag_names[mid] = set([tag.tag_name for tag in tags])
    return mids_to_tag_names


def replace_tags_by_motifs(session, ids_to_tags):
    ids_to_current_tags = find_tags_by_motif_ids(session, ids_to_tags.keys())
    ids_to_add_tags = {}
    tag_ids_to_remove = set()
    for mid in ids_to_tags:
        current_tags = ids_to_current_tags.get(mid, set())
        current_tag_names_to_id = dict([(tag.tag_name, tag.id) for tag in current_tags])
        current_tag_names = set(current_tag_names_to_id.keys())
        tags = set(ids_to_tags[mid])
        ids_to_add_tags[mid] = tags - current_tag_names
        tag_names_to_del = current_tag_names - tags
        tag_ids_to_remove |= set([current_tag_names_to_id[tag_name] for tag_name in tag_names_to_del])

    add_tags_by_map(session, ids_to_add_tags)
    remove_by_tag_ids(session, tag_ids_to_remove)


def add_tags_by_names_to_motifs(session, names_to_motif_ids):
    tag_refs_to_add = []
    for name, motif_ids in names_to_motif_ids.items():
        tag_refs = find_tag_refs_by_tag_name_and_motif_ids(session, name, motif_ids)
        exists_motif_ids = set([tag_ref.target_id for tag_ref in tag_refs])
        for motif_id_to_add in set(motif_ids) - exists_motif_ids:
            tag_refs_to_add.append(TagRefEntity(motif_id_to_add, name))
    session.add_all(tag_refs_to_add)
