import os
import re
import uuid
import sqlite3
import pickle
import json
from glob import glob


def get_text_block(self, tag):
    found = False
    block = ''
    handle = open(os.path.join(os.curdir, "Y2Hreadme.txt"), 'r')
    for line in handle.readlines():
        if found:
            if line.strip() == ">>>END":
                break
            block += line
        else:
            if line.strip() == ">>>" + tag:
                found = True
    handle.close()
    return block


deepn_data = json.load(open("genome.json"))
exists = os.path.exists("deepn.sqlite")
con = sqlite3.connect("deepn.sqlite")
cur = con.cursor()
if (exists):
    cur.execute("DROP TABLE deepn")
    cur.execute("DROP TABLE exons")
    cur.execute("DROP TABLE refseq")

cur.execute(
    "CREATE TABLE deepn(deepn_id, genome, name, prey_library, orf_data, orf_lookup_table, blast_db, junction_presequence)"
)
cur.execute(
    "CREATE TABLE exons(exon_id, genome, chromosome, name, start_bp, end_bp, strand)"
)
cur.execute("CREATE TABLE refseq(refseq_id, exon_id, name)")

for f in glob("*.p"):
    print(f)
    gid = os.path.splitext(f)[0]
    for dda in deepn_data[f]:
        did = str(uuid.uuid4())[-12:]
        cur.execute(
            "INSERT INTO deepn(deepn_id, genome, name, prey_library, orf_data, orf_lookup_table, blast_db, junction_presequence) values (?,?,?,?,?,?,?,?)",
            (did, gid, dda["genome"], dda["prey_library"], dda["orf_data"],
             dda["orf_lookup_table"], dda["blast_db"],
             dda["junction_presequence"]))
    data = pickle.load(open(f, 'rb'))
    for chr in data.keys():
        cid = "chr{}".format(chr)

        for k in data[chr].keys():
            eid = str(uuid.uuid4())[-12:]
            if len(k) == 7:
                start, stop, _, d, name, nm_nums, t = k
            elif len(k) == 6:
                start, stop, _, d, name, t = k
            cur.execute(
                "INSERT INTO exons(exon_id, chromosome, genome, name, start_bp, end_bp, strand) values (?,?,?,?,?,?,?)",
                (eid, cid, gid, name, start, stop, d))
            if len(k) == 7:
                for nm in list(set(nm_nums)):
                    nid = str(uuid.uuid4())[-12:]
                    cur.execute(
                        "INSERT INTO refseq(refseq_id, exon_id, name) values (?,?,?)",
                        (nid, eid, nm))
        con.commit()
cur.execute(
    "CREATE INDEX exon_index on exons (genome, chromosome, start_bp, end_bp)")
cur.close()
con.close()