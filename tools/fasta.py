
class Seq(object):
    def __init__(self):
        self.seq_id = None
        self.seq = None
        self.description = None


def _analyse_first_line(line):
    line.replace(">", "")


def read_in_all_seqs(buf):
    seq_ids_to_seq = {}
    seq = None
    for line in buf.splitlines():
        if line.startswith(">"):
            if seq is not None:
                seq_ids_to_seq[seq.seq_id] = seq
            seq = Seq()
