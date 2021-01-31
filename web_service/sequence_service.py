
import logging
from bottle import get
from bottle import request
import dao
from web_service.service_utils import *
from web_service.motif_service import motif_entity_output
from web_service.nsite_service import nsites_entity_output
from tools.exception import *


__all__ = ['get_sequences', ]


def _get_offset_arg(filters):
    try:
        filters['offseteq'] = int(request.query.offset)
    except ValueError:
        pass
    else:
        return

    try:
        filters['offsetlt'] = int(request.query.offsetlt)
    except ValueError:
        pass

    try:
        filters['offsetgt'] = int(request.query.offsetgt)
    except ValueError:
        pass


def _get_lrr_count_arg(filters):
    try:
        filters['lrr_count_eq'] = int(request.query.lrr_count)
    except ValueError:
        pass
    else:
        return

    try:
        filters['lrr_count_lt'] = int(request.query.lrr_count_lt)
    except ValueError:
        pass

    try:
        filters['lrr_count_gt'] = int(request.query.lrr_count_gt)
    except ValueError:
        pass


def _get_page_arg():
    page = None
    size = None
    try:
        page = int(request.query.page)
        size = int(request.query.size)
    except ValueError:
        pass
    if page is None or type(page) is not int:
        page = 0
    if size is None or type(size) is not int:
        size = 40
    if size > 400:
        size = 400
    if page < 0:
        raise ValidationError(str.format(ErrorCode.INVALID_PARA, 'pageIndex', "less than 0"))
    if size < 0:
        raise ValidationError(str.format(ErrorCode.INVALID_PARA, 'pageSize', "less than 0"))
    return page, size


def _get_keyword_arg():
    keyword = request.query.keyword
    if keyword is not None:
        keyword = keyword.strip()
        if len(keyword) == 0:
            keyword = None
    return keyword


def _get_species_arg(filters):
    species = request.query.dict.get('species[]', None)
    if species is not None:
        if isinstance(species, str) and len(species) == 0:
            return
        if not isinstance(species, list):
            raise ValidationError(str.format(ErrorCode.INVALID_PARA, 'species', 'must be a list'))
        if len(species) != 0:
            filters['species'] = set(species)


def sequence_entity_to_output(seq, motifs, nsites, ids_to_tag_names):
    return {
        'id': seq.id,
        'sequence_id': seq.seq_id,
        'seq': seq.seq,
        'ss': seq.ss,
        'motifs_16': motif_entity_output(motifs, ids_to_tag_names),
        'nsites': nsites_entity_output(nsites)
    }


@get('/version/<version>/sequences')
def get_sequences(version):
    filters = {}
    filters['page_index'], filters['page_size'] = _get_page_arg()
    filters['keyword'] = _get_keyword_arg()
    _get_offset_arg(filters)
    _get_lrr_count_arg(filters)
    _get_species_arg(filters)

    filters['page_index'] = filters['page_index'] * filters['page_size']
    logging.debug(str.format("Filters: {}", filters))
    seqs, total = dao.query_sequences(filters, get_version_arg(version))

    # find motifs and nsites
    seq_ids = [seq.seq_id for seq in seqs]
    seq_ids_to_motifs = dao.find_motifs_by_seq_ids(seq_ids, get_version_arg(version))
    seq_ids_to_nsites = dao.find_nsites_by_seq_ids(seq_ids)

    # find tags, the overlap mark is tagged on version 1 motifs
    if get_version_arg(version) == 1:
        ids_to_tag_names = dao.find_tags_by_motif_ids(set([m.id for motifs in seq_ids_to_motifs.values() for m in motifs]))
    else:
        ids_to_tag_names = {}

    result = {'sequences': [], 'total': total}
    for seq in seqs:
        result['sequences'].append(sequence_entity_to_output(seq,
                                                             seq_ids_to_motifs.get(seq.seq_id, []),
                                                             seq_ids_to_nsites.get(seq.seq_id, []),
                                                             ids_to_tag_names))
    return response_ok(result, True)
