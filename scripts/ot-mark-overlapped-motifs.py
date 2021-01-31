import logging
import dao
from tools import motifs as motif_util


MOTIF_VERSION = 1


def _get_motifs_by_seq_id(seq_id):
    seq_ids_to_motifs = dao.find_motifs_by_seq_ids([seq_id, ], MOTIF_VERSION)
    return seq_ids_to_motifs.get(seq_id, [])


def _find_overlap(seq_id, motifs):
    no_overlap_motifs = motif_util.found_no_overlapped_motifs(motifs, 16)
    if len(no_overlap_motifs) == len(motifs):
        return []
    offsets_to_motif = dict([(m.offset, m) for m in motifs])
    prev_motif_offsets = list(offsets_to_motif.keys())
    removed_offsets = list(set(prev_motif_offsets) - set([m.offset for m in no_overlap_motifs]))

    prev_motif_offsets.sort()
    removed_offsets.sort()
    removed_ids = [offsets_to_motif.get(offset).id for offset in removed_offsets]
    logging.warning(str.format("Found overlapping in seq {}, prev offsets {}, removed {}(ids {})",
                               seq_id, prev_motif_offsets, removed_offsets, removed_ids))
    return removed_ids


def main():
    logging.basicConfig(level=logging.WARNING)
    with dao.query_session() as session:
        seq_ids = dao.sequence.find_all_seq_ids(session)

    overlap_ids = set()
    for seq_id in seq_ids:
        logging.info(str.format("Begin to analyse overlapping info for seq {}", seq_id))
        motifs = _get_motifs_by_seq_id(seq_id)
        ids = _find_overlap(seq_id, motifs)
        overlap_ids.update(ids)
    logging.warning(str.format("All overlapped ids {}, write to db...", overlap_ids))
    dao.add_tags_by_names_to_ids({'inner.overlap': overlap_ids}, MOTIF_VERSION)


if __name__ == '__main__':
    main()
