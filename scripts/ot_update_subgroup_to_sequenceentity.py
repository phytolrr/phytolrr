import os
import re
import logging
import roman
import dao

SUBGROUP_FILES_ROOT_PATH = '/Users/shengnan/code/cts/nsites/nsites/phy_trees/valid_tree_190510/iqtree_analisys_gt8/ids_SG_try2_190514'


def _get_subgroup_from_match(ma):
    sg_str = ma.group('sg')
    sg_left = sg_str.split('-')[0].split('_')[0]
    sg_right = sg_str.split('-')[-1].split('_')[-1]
    if sg_left == sg_right:
        sg_right = None
    sg_left = roman.fromRoman(sg_left.upper())
    if sg_right is None:
        return str(sg_left)
    else:
        return str(sg_left) + '-' + sg_right


def read_in_subgroup_seqs(path, file_pattern):
    fre = re.compile(file_pattern)
    seq_ids_to_sg = {}
    sgs_to_seq_ids = {}
    for file_name in os.listdir(path):
        file_path = os.path.join(path, file_name)
        if not os.path.isfile(file_path):
            continue
        ma = fre.match(file_name)
        if ma is None:
            continue
        sg = _get_subgroup_from_match(ma)
        if sg not in sgs_to_seq_ids:
            sgs_to_seq_ids[sg] = []
        logging.info(str.format("Begin to read in subgroup {} file {}", sg, file_path))
        with open(file_path, 'r') as f:
            for line in f.readlines():
                seq_id = line.strip()
                if seq_id in seq_ids_to_sg and seq_ids_to_sg.get(seq_id) != sg:
                    logging.error(str.format("The seq id {} was found in more than one subgroup {}, {}",
                                             seq_id, seq_ids_to_sg.get(seq_id), sg))
                    return None
                else:
                    seq_ids_to_sg[seq_id] = sg
                    sgs_to_seq_ids[sg].append(seq_id)
    return sgs_to_seq_ids


def main():
    sgs_to_seq_ids = read_in_subgroup_seqs(SUBGROUP_FILES_ROOT_PATH, r'^sg_id_(?P<sg>.*)$')
    logging.info(str.format("Read in seq ids and subgroup {}", sgs_to_seq_ids))
    with dao.session_scope() as session:
        for sg, seq_ids in sgs_to_seq_ids.items():
            dao.sequence.update_subgroup_by_seq_ids(session, seq_ids, sg)


if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO)
    main()