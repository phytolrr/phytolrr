
from bottle import post, put, delete
from bottle import request
from sqlalchemy.orm.exc import NoResultFound

import dao
import settings
from web_service.service_utils import *
from tools.exception import *
from tools import pssm_matrix
from tools import motifs as motif_tools


def motif_entity_output(motifs, ids_to_tag_names):
    output = []
    if motifs is None:
        return output
    for m in motifs:
        if dao.tag.OVERLAP_TAG in ids_to_tag_names.get(m.id, []):
            continue
        # todo: 数据库中刷新该列
        if m.manually_add is not None and m.manually_add:
            manually_add = True
        else:
            manually_add = False
        output.append({
            'id': m.id,
            'manually_add': manually_add,
            'offset': m.offset,
            'score': m.score,
            'probability': m.probability,
            'fdr_probability': m.fdr_probability,
            'false_discovery': m.false_discovery
        })
    output.sort(key=lambda k: k['offset'])
    return output


def get_and_check_int(num_str, message):
    try:
        return int(num_str)
    except ValueError:
        raise ValidationError(message)


def get_and_check_motif_id(mid_str, version):
    mid = get_and_check_int(mid_str, str.format(ErrorCode.INVALID_PARA, "mid", mid_str))
    if mid < 0:
        raise ValidationError(str.format(ErrorCode.INVALID_PARA, "mid", 'less than 0'))
    with dao.query_session() as session:
        motifs = dao.motif.find_motifs_by_ids(session, [mid], version)
    if len(motifs) == 0:
        raise ValidationError(str.format(ErrorCode.INVALID_PARA, 'mid', 'the motif does not exists'))
    return motifs[0]


def get_and_check_sid(sid):
    sid = get_and_check_int(sid, str.format(ErrorCode.INVALID_PARA, 'sid', sid))
    try:
        return dao.find_seq_by_sid(sid)
    except NoResultFound as e:
        raise ValidationError(str.format(ErrorCode.OBJECT_NOT_EXISTS, sid))


def get_and_check_offset(seq, version):
    payload = request.json
    if payload is None:
        raise ValidationError(str.format(ErrorCode.PARA_NOT_EXISTS, "offset"))
    offset = payload.get('offset', None)
    if offset is None:
        raise ValidationError(str.format(ErrorCode.PARA_NOT_EXISTS, "offset"))
    offset = get_and_check_int(offset, str.format(ErrorCode.INVALID_PARA, 'offset', offset))
    if offset < 0:
        raise ValidationError(str.format(ErrorCode.INVALID_PARA, 'offset', offset))
    if offset + 16 > len(seq.seq):
        raise ValidationError(str.format(ErrorCode.INVALID_PARA, 'offset', offset))

    seq_ids_to_motifs = dao.find_motifs_by_seq_ids([seq.seq_id, ], version)
    motifs = seq_ids_to_motifs.get(seq.seq_id, [])
    poses = set()
    wrong_poses = set()
    for m in motifs:
        if m.correct:
            poses.update(range(m.offset - 15, m.offset + 16))
        else:
            wrong_poses.add(m.offset)
    if offset in poses:
        raise ValidationError(ErrorCode.OFFSET_OVERLAP)
    if offset in wrong_poses:
        raise ValidationError(ErrorCode.OFFSET_EXISTS_WRONG)
    return offset


def _get_motifs_by_offset(seq_id, offsets, version):
    motifs = dao.find_motifs_by_offsets(seq_id, offsets, version)
    return motif_entity_output(motifs, {})


def _get_input_false_discovery():
    payload = request.json
    if payload is None:
        raise ValidationError(str.format(ErrorCode.PARA_NOT_EXISTS, 'false_discovery'))
    false_discovery = payload.get("false_discovery", None)
    if false_discovery is None:
        raise ValidationError(str.format(ErrorCode.PARA_NOT_EXISTS, 'false_discovery'))
    if not isinstance(false_discovery, bool):
        raise ValidationError(str.format(ErrorCode.INVALID_PARA, 'tag_names', 'the type of tag_names must be a boolean'))
    return false_discovery


def is_overlapped(offsets):
    poses = set()
    for offset in offsets:
        prev_len = len(poses)
        poses.update(range(offset, offset + 16))
        after_len = len(poses)
        if after_len - prev_len != 16:
            return True
    return False


MATRIX = None

if settings.mode == settings.MODE_DEV:
    @post('/version/<version>/sequences/<sid>/motifs')
    def add_manually_motif(version, sid):
        global MATRIX
        if MATRIX is None:
            MATRIX = pssm_matrix.calc_pssm_matrix(dao.find_baseline_motifs(1, False))
        try:
            seq = get_and_check_sid(sid)
            offset = get_and_check_offset(seq, get_version_arg(version))
            score = motif_tools.calc_pssm_score(seq.seq[offset:offset+16], MATRIX)
            probability = motif_tools.calc_probability_by_score(score)
            result = dao.add_manually_motif(seq.seq_id, offset, get_version_arg(version), score, probability)
            if result is None:
                return response_ok({'motifs_16': _get_motifs_by_offset(seq.seq_id, [offset, ], get_version_arg(version))})
            else:
                return response_error(result)
        except ValidationError as e:
            return response_error(e.message)


    @delete('/version/<version>/sequences/<sid>/motifs/<mid>')
    def delete_manually_motif(version, sid, mid):
        try:
            motif_entity = get_and_check_motif_id(mid, get_version_arg(version))
            dao.remove_manually_motif(motif_entity.id, get_version_arg(version))
            return response_ok({})
        except ValidationError as e:
            return response_error(e.message)


    @put('/version/<version>/sequences/<sid>/motifs/<mid>/false_discovery')
    def change_motif_false_discovery(version, sid, mid):
        try:
            version = get_version_arg(version)
            seq = get_and_check_sid(sid)
            motif = get_and_check_motif_id(mid, version)

            if seq.seq_id != motif.seq_id:
                raise ValidationError(str.format("The seq_id {} in sequence is different with motif {}",
                                                 seq.seq_id, motif.seq_id))
            if motif.manually_add:
                raise ValidationError(str.format("The motif {} offset {} in sequence {} is manually added,"
                                                 " and should not be tagged wrong again. The motif can be delete directly",
                                                 mid, motif.offset, motif.seq_id))

            false_discovery = _get_input_false_discovery()
            if not false_discovery:
                seq_ids_to_motifs = dao.find_motifs_by_seq_ids([motif.seq_id,], version)
                if is_overlapped([m.offset for m in seq_ids_to_motifs[motif.seq_id] if m.correct] + [motif.offset]):
                    raise ValidationError(str.format("The motif {} can not be unmark from wrong because overlapping was found",
                                                     mid))
            if version < 2:
                raise ValidationError("The version 1 is no longer support marking wrong")
            dao.update_false_discovery_by_motif(mid, false_discovery, version)

            motif_entity = dao.find_motif_by_mid(mid, version)
            if motif_entity is None:
                return response_error(str.format("Internal error, can not find the motif {}", mid))
            motif_output = motif_entity_output([motif_entity], {})
            return response_ok(motif_output[0])

        except ValidationError as e:
            return response_error(e.message)
