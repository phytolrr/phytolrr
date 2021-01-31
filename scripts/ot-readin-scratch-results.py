# -*- coding: utf-8 -*
'''本脚本读取scratch软件的输出结果，并将其刷新到数据库中，对应的sequenceentity表中的ss,ss8,acc,acc20列中。
本脚本没有参数，scratch软件的输入结果所在目录通过脚本中的PATHS变量指定，可以指定多个。
本脚本默认scratch的数据结果包含四个文件：
<seq>.fasta.out.ss
<seq>.fasta.out.ss8
<seq>.fasta.out.acc
<seq>.fasta.out.acc20
约束：
1. seq字段不强制与序列的ID一致，但是seq后面的字段必须为fasta.out.结果类型。
2. 每个文件均为fasta格式
3. 每个文件只能包含一条序列
'''
import os
import logging
from Bio import SeqIO
import dao

PATH_START = r'../scratch_analysis/ectodomain_seqs_for_scracth/run_'
PATHS = []

for index in range(1, 41):
    PATHS.append(str.format("{}{}", PATH_START, index))

# PATHS = [r'../scratch_analysis/ectodomain_seqs_for_scracth/run_1']

class ScratchResult(object):
    def __init__(self):
        self.acc = None
        self.acc20 = None
        self.ss = None
        self.ss8 = None

    def is_scratch_result_valid(self):
        if self.acc is None:
            return False
        if self.acc20 is None:
            return False
        if self.ss is None:
            return False
        if self.ss8 is None:
            return False
        return True

    def __str__(self):
        return self.__dict__.__str__()

    def __repr__(self):
        return self.__str__()


class Acc20Seq(object):
    def __init__(self):
        self.seq = None
        self.id = None


def read_in_scratch_results_file_by_path(path):
    seq_ids_to_files = {}
    for root_path, dir_names, file_names in os.walk(path):
        for file_name in file_names:
            eles = file_name.split('.fasta.out.')
            if len(eles) != 2:
                logging.error(str.format("Unexpected file name {}", file_name))
                continue
            seq_id = eles[0]
            if not seq_ids_to_files.has_key(seq_id):
                seq_ids_to_files[seq_id] = ScratchResult()
            seq_ids_to_files[seq_id].__dict__[eles[1]] = os.path.join(root_path, file_name)
    return seq_ids_to_files


def read_in_acc20(file_path):
    with open(file_path, 'r') as f:
        lines = f.readlines()
    if len(lines) < 2:
        logging.error(str.format("Failed to parse acc20 file {}, less than two lines", file_path))
        return None
    id_line = lines[0].strip()
    if not id_line.startswith('>'):
        logging.error(str.format("Failed to parse acc20 file {}, the first line is not starts with >", file_path))
        return None
    seq = Acc20Seq()
    seq.id = id_line[1:]
    seq.seq = ' '.join(lines[1:]).split()
    return seq


def read_in_scratch_result(result):
    ss = SeqIO.read(result.ss, 'fasta')
    ss8 = SeqIO.read(result.ss8, 'fasta')
    acc = SeqIO.read(result.acc, 'fasta')
    acc20 = read_in_acc20(result.acc20)
    if acc20 is None:
        return
    if ss.id != ss8.id or ss.id != acc.id or ss.id != acc20.id:
        logging.error(str.format("The result for seq {} is invalid, the seq id is different", ss.id))
        return
    if len(ss.seq) != len(ss8.seq) or len(ss.seq) != len(acc.seq) or len(ss.seq) != len(acc20.seq):
        logging.error(str.format("The result for seq {} is invalid, the seq length is different", ss.id))
        return

    seq_entity = None
    try:
        seq_entity = dao.find_seq_by_id(ss.id)
    except Exception as e:
        logging.error(str.format("Failed to save result for seq {}, can not find the record in the db", ss.id))
        return
    if seq_entity is None:
        logging.error(str.format("Failed to save result for seq {}, can not find the record in the db", ss.id))
        return
    if len(ss.seq) != len(seq_entity.seq):
        logging.error(str.format("The result for seq {} is invalid, the seq length is different with db", ss.id))
        return
    logging.info("Update scratch result for seq " + ss.id)
    dao.update_seq_for_scratch_result(ss.id, str(ss.seq), str(ss8.seq), str(acc.seq), ' '.join(acc20.seq))


def main():
    logging.basicConfig(level=logging.INFO)
    logging.debug(str.format("All paths {}", PATHS))
    seq_ids_to_files = {}
    for path in PATHS:
        one_path_results = read_in_scratch_results_file_by_path(path)
        seq_ids_to_files = dict(seq_ids_to_files, **one_path_results)
    for k, v in seq_ids_to_files.items():
        read_in_scratch_result(v)
        pass

main()
