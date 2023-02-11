import os
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
    cur.execute("DROP TABLE genome")
    cur.execute("DROP TABLE deepn")
    cur.execute("DROP TABLE chromosomes")
    cur.execute("DROP TABLE exons")
    cur.execute("DROP TABLE refseq")

cur.execute("CREATE TABLE genome(genome_id, name)")
cur.execute(
    "CREATE TABLE deepn(deepn_id, genome_id, name, prey_library, orf_data, orf_lookup_table, blast_db, junction_presequence)"
)
cur.execute("CREATE TABLE chromosomes(chromosome_id, genome_id, name)")
cur.execute(
    "CREATE TABLE exons(exon_id, genome_id, chromosome_id, name, start_bp, end_bp, strand)"
)
cur.execute("CREATE TABLE refseq(refseq_id, exon_id, name)")

for f in glob("*.p"):
    print(f)
    gid = str(uuid.uuid4())[-12:]
    cur.execute("INSERT INTO genome(genome_id, name) values (?,?)",
                (gid, os.path.splitext(f)[0]))
    for dda in deepn_data[f]:
        did = str(uuid.uuid4())[-12:]
        cur.execute(
            "INSERT INTO deepn(deepn_id, genome_id, name, prey_library, orf_data, orf_lookup_table, blast_db, junction_presequence) values (?,?,?,?,?,?,?,?)",
            (did, gid, dda["genome"], dda["prey_library"], dda["orf_data"],
             dda["orf_lookup_table"], dda["blast_db"],
             dda["junction_presequence"]))
    data = pickle.load(open(f, 'rb'))
    for chr in data.keys():
        cid = str(uuid.uuid4())[-12:]
        cur.execute(
            "INSERT INTO chromosomes(chromosome_id, genome_id, name) values (?,?,?)",
            (cid, gid, chr))
        for k in data[chr].keys():
            eid = str(uuid.uuid4())[-12:]
            if len(k) == 7:
                start, stop, _, d, name, nm_nums, t = k
            elif len(k) == 6:
                start, stop, _, d, name, t = k
            cur.execute(
                "INSERT INTO exons(exon_id, chromosome_id, genome_id, name, start_bp, end_bp, strand) values (?,?,?,?,?,?,?)",
                (eid, cid, gid, name, start, stop, d))
            if len(k) == 7:
                for nm in list(set(nm_nums)):
                    nid = str(uuid.uuid4())[-12:]
                    cur.execute(
                        "INSERT INTO refseq(refseq_id, exon_id, name) values (?,?,?)",
                        (nid, eid, nm))
        con.commit()
cur.close()
con.close()