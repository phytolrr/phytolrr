# -*- coding: utf-8 -*
'''本脚本读取数据库中的所有序列，找出N糖位点，并将N糖位点写回数据库。'''
import dao
import re
import logging


def _get_ntype_from_match(match_str):
    if match_str[-1] == 'S':
        return dao.nsite.S
    elif match_str[-1] == 'T':
        return dao.nsite.T
    else:
        raise Exception(str.format("Unexpected match {}", match_str))


def search_seq(seq, regex_str):
    r = re.compile(regex_str)
    start_positions = []
    ntypes = []
    cur_pos = 0
    while True:
        m = r.search(seq)
        if m is None:
            break
        start_positions.append(m.start() + cur_pos)
        ntypes.append(_get_ntype_from_match(m.group()))

        end_pos = m.end()
        seq = seq[end_pos:]
        cur_pos += end_pos

    return start_positions, ntypes


def update_seq(seq_id, seq, regex_str):
    start_positions, ntypes = search_seq(seq, regex_str)
    if len(start_positions) == 0:
        logging.warning("No n-sites found in seq " + seq_id)
        return
    logging.info(str.format("Update nsites by seq_id {}, pos {}", seq_id, list(zip(start_positions, ntypes))))
    dao.replace_nsites_by_seq_id(seq_id, start_positions, ntypes)


def flush_all_seqs():
    regex_str = r"N[A-OQ-Z][ST]"
    with dao.session_scope() as session:
        seqs = dao.sequence.find_all_seqs(session)
        logging.info(str.format("Total seq count {}", len(seqs)))
        for seq in seqs:
            update_seq(seq.seq_id, seq.seq, regex_str)


def main():
    logging.basicConfig(level=logging.DEBUG)
    flush_all_seqs()


main()
