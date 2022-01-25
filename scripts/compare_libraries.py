import sys
import argparse
import json
import tbprofiler as tbp
def main(args):
    conf = tbp.get_conf_dict(args.db)
    lt2gene = tbp.rv2genes(conf["bed"])
    db1_raw = json.load(open(args.json1))
    db2_raw = json.load(open(args.json2))

    db1 = set()
    for gene in db1_raw:
        for var in db1_raw[gene]:
            for d in [ann["drug"] for ann in db1_raw[gene][var]["annotations"] if ann["type"]=="drug"]:
                db1.add((gene,var,d))
                
    db2 = set()
    for gene in db2_raw:
        for var in db2_raw[gene]:
            for d in [ann["drug"] for ann in db2_raw[gene][var]["annotations"] if ann["type"]=="drug"]:
                db2.add((gene,var,d))

    for gene,var,drug in db2-db1:
        print(lt2gene[gene],var,drug,sep="\t")

parser = argparse.ArgumentParser(description='tbprofiler script',formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument('--json1',type=str,help='')
parser.add_argument('--json2',type=str,help='')
parser.add_argument('--db',default="tbdb",type=str,help='')
parser.set_defaults(func=main)
args = parser.parse_args()
args.func(args)