from bottle import post, request
from tools.exception import ValidationError, ErrorCode
from web_service import service_utils
from tools import motifs as motif_tool, pssm_matrix
import dao
import logging


MAX_SEQ_LENGTH = 8000
MATRIX = None


def _get_matrix():
    global MATRIX
    if MATRIX is None:
        motifs = dao.find_baseline_motifs(baseline_version=1, with_wrong=False)
        logging.info(str.format("Baseline LRR motifs for lrr-service( count {}): {}", len(motifs), motifs))
        MATRIX = pssm_matrix.calc_pssm_matrix(motifs)
    return MATRIX


def _get_sequence():
    payload = request.json
    if payload is None:
        raise ValidationError(ErrorCode.PARA_NOT_EXISTS, "seq")
    seq = payload.get("seq", None)
    if seq is None:
        raise ValidationError(ErrorCode.PARA_NOT_EXISTS, "seq")
    if len(seq) > MAX_SEQ_LENGTH:
        raise ValidationError(str.format("The length of sequence must be less than {}", MAX_SEQ_LENGTH))
    return seq


@post("/find-lrr")
def find_lrr():
    try:
        seq = _get_sequence()
        matrix = _get_matrix()
        motifs = motif_tool.lrr_search(matrix, seq)
        motifs = motif_tool.found_no_overlapped_motifs(motifs, 16)
        motifs.sort(key=lambda m:m.offset)
        return service_utils.response_ok({'LRRs': [m.__dict__ for m in motifs]})
    except ValidationError as e:
        return service_utils.response_error(e.message)
