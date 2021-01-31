from bottle import post, get
from bottle import request
import dao
from tools.exception import *
from web_service.service_utils import *
from web_service.motif_service import get_and_check_motif_id


def is_overlapped(offsets):
    poses = set()
    for offset in offsets:
        prev_len = len(poses)
        poses.update(range(offset, offset + 16))
        after_len = len(poses)
        if after_len - prev_len != 16:
            return True
    return False


def _check_motif_can_be_tag(mid, tag_names, version):
    with dao.query_session() as session:
        ms = dao.motif.find_motifs_by_ids(session, [mid], version)
    if ms is None or len(ms) == 0:
        raise ValidationError(str.format(ErrorCode.OBJECT_NOT_EXISTS, mid))


def get_and_check_tag_names():
    payload = request.json
    if payload is None:
        raise ValidationError(str.format(ErrorCode.PARA_NOT_EXISTS, 'tag_names'))
    tag_names = payload.get('tag_names', None)
    if tag_names is None:
        raise ValidationError(str.format(ErrorCode.PARA_NOT_EXISTS, 'tag_names'))
    if not isinstance(tag_names, list):
        raise ValidationError(str.format(ErrorCode.INVALID_PARA, 'tag_names', 'the type of tag_names must be a list'))
    return set(tag_names)


def _get_motif_tags_without_check(mids):
    mids_to_tag_names = dao.find_tags_by_motif_ids(mids)
    for mid, tags in mids_to_tag_names.items():
        mids_to_tag_names[mid] = list(tags)
        mids_to_tag_names[mid].sort()
    return mids_to_tag_names


@get('/version/<version>/sequences/<sid>/motifs/<mid>/tags')
def get_motif_tags(version, sid, mid):
    return response_ok(_get_motif_tags_without_check([mid]))


@post('/version/<version>/sequences/<sid>/motifs/<mid>/tags')
def add_tags_to_motif(version, sid, mid):
    try:
        tag_names = get_and_check_tag_names()
        mid = get_and_check_motif_id(mid, get_version_arg(version)).id
        _check_motif_can_be_tag(mid, tag_names, get_version_arg(version))
        dao.replace_tags_by_motifs({mid: tag_names}, get_version_arg(version))
        return response_ok(_get_motif_tags_without_check([mid]))
    except ValidationError as e:
        return response_error(e.message)