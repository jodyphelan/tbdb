import sys
import argparse
import csv
import tbprofiler as tbp
import re
import json
drug_convert = {
    "AMI":"amikacin",
    "BDQ":"bedaquiline",
    "CAP":"capreomycin",
    "CFZ":"clofazimine",
    "DLM":"delamanid",
    "EMB":"ethambutol",
    "ETH":"ethionamide",
    "INH":"isoniazid",
    "KAN":"kanamycin",
    "LEV":"levofloxacin",
    "LZD":"linezolid",
    "MXF":"moxifloxacin",
    "PZA":"pyrazinamide",
    "RIF":"rifampicin",
    "STM":"streptomycin"
}

source = "https://www.who.int/publications/i/item/9789240028173"

def main(args):
    genes = set()
    dr_mutations = []
    other_mutations = []
    for row in csv.DictReader(open(args.csv)):
        if row["fail"]=="True": 
            continue
        drug = drug_convert[row["drug"]]
        genes.add((row["gene"],drug))
        if row["hgvs"][0]=="r":
            r = re.search("r.([0-9]+)([acgt]+)>([acgt]+)",row["hgvs"])
            row["hgvs"] = f"n.{r.group(1)}{r.group(2).upper()}>{r.group(3).upper()}"
        if row["classification"] in ["1","2"]:  
            dr_mutations.append({
                "Gene":row["gene"],
                "Mutation":row["hgvs"],
                "Drug":drug,
                "Confers":"resistance",
                "Confidence":row["classification"],
                "Interaction":"",
                "Literature":source
            })
        else:
            other_mutations.append({
                "Gene":row["gene"],
                "Mutation":row["hgvs"],
                "Type":"resistance_association_confidence",
                "Info":f"drug={drug};confidence={row['classification']};source={source}"
            })

    dr_mutations = [json.loads(y) for y in set(json.dumps(x) for x in dr_mutations)]
    other_mutations = [json.loads(y) for y in set(json.dumps(x) for x in other_mutations)]
    genes_list = []
    for item in genes:
        genes_list.append({
            "Gene":item[0],
            "Drug":item[1],
            "Literature":source
        })

    with open(args.out+".dr_mutations.csv","w") as O:
        writer = csv.DictWriter(O,fieldnames=list(dr_mutations[0]))
        writer.writeheader()
        writer.writerows(dr_mutations)
    
    with open(args.out+".other_mutations.csv","w") as O:
        writer = csv.DictWriter(O,fieldnames=list(other_mutations[0]))
        writer.writeheader()
        writer.writerows(other_mutations)

    with open(args.out+".watchlist.csv","w") as O:
        writer = csv.DictWriter(O,fieldnames=list(genes_list[0]))
        writer.writeheader()
        writer.writerows(genes_list)

parser = argparse.ArgumentParser(description='tbprofiler script',formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument('--csv',type=str,help='',required = True)
parser.add_argument('--out',type=str,help='',required = True)
parser.set_defaults(func=main)
args = parser.parse_args()
args.func(args)