import logging
logging.basicConfig(level=logging.INFO)
import dao

FILE_PATH = r'C:\Users\titan\Desktop\ambicious_motifs_unconfirmed.txt'
MOTIF_VERSION = 2


with open(FILE_PATH, 'r') as f:
    lines = f.readlines()


ids = []
for line in lines:
    eles = line.split(',')
    seq_id = eles[0].split()[-1]
    offset = int(eles[1].split()[-1])
    m = dao.find_motifs_by_offsets(seq_id, [offset], MOTIF_VERSION)
    assert len(m) == 1
    m = m[0]
    logging.info(str.format("Begin to delete motif id {}, seq id {}, offset {}", m.id, m.seq_id, m.offset))
    ids.append(m.id)

with dao.session_scope() as session:
    dao.motif.delete_manually_motifs_by_ids(session, ids, MOTIF_VERSION)
