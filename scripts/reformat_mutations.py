import argparse
import re
import csv
import pathogenprofiler as pp
from uuid import uuid4
import os 

class gene_class:
    def __init__(self,name,locus_tag,strand,start,end,length):
        self.name = name
        self.locus_tag = locus_tag
        self.strand = strand
        self.start = start
        self.end = end
        self.length = length


def load_gff(gff):
    genes = []
    for l in open(gff):
        if l[0]=="#": continue
        fields = l.rstrip().split()
        if fields[2] not in ["gene","rRNA_gene"]: continue
        strand = fields[6]
        p1 = int(fields[3])
        p2 = int(fields[4])
        gene_length = p2-p1+1
        re_obj = re.search("Name=(.+?);",l)
        gene_name = re_obj.group(1) if re_obj else "NA"
        re_obj = re.search("ID=gene:(.+?);",l)
        locus_tag = re_obj.group(1) if re_obj else "NA"
        start = p1 if strand=="+" else p2
        end =  p2 if strand=="+" else p1
        tmp = gene_class(gene_name,locus_tag,strand,start,end,gene_length)
        genes.append(tmp)
    return genes

def get_ann(variants):
    uuid = str(uuid4()) #"463545ef-71fc-449b-8f4e-9c907ee6fbf5"
    with open(uuid,"w") as O:
        O.write('##fileformat=VCFv4.2\n')
        O.write('##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n')
        O.write('##contig=<ID=Chromosome,length=4411532>\n')
        O.write('#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\ttest\n')
        for var in variants.values():
            O.write("Chromosome\t%(pos)s\t.\t%(ref)s\t%(alt)s\t255\t.\t.\tGT\t1\n" % var) 
    results = {}
    keys = list(variants.keys())
    vals = list(variants.values())
    i = 0
    for l in pp.cmd_out(f"snpEff ann Mycobacterium_tuberculosis_h37rv {uuid}"):
        if l[0]=="#": continue
        row = l.strip().split()
        for ann in row[7].split(","):
            a = ann.split("|")
            if vals[i]["gene"] in [a[3],a[4]]:
                results[keys[i]] = a[9] if vals[i]["type"]=="nucleotide" else a[10]
        i+=1
    os.remove(uuid)
    return results

    


def main(args):
    genes = load_gff(args.gff)
    refseq = pp.fasta(args.ref).fa_dict

    mutations  =  {}
    converted_mutations = {}
    for row in csv.DictReader(open(args.csv)):
        gene = [g for g in genes if g.name==row["Gene"] or g.locus_tag==row["Gene"]][0]
        mut = None
        r = re.search("r.([0-9]+)([acgt]+)>([acgt]+)",row["Mutation"])
        if r:
            converted_mutations[(row["Gene"],row["Mutation"])] = f"n.{r.group(1)}{r.group(2).upper()}>{r.group(3).upper()}"
        r = re.search("p\..+",row["Mutation"])
        if r:
            converted_mutations[(row["Gene"],row["Mutation"])] = row["Mutation"]
        r = re.search("c.-[0-9]+[ACGT]>[ACGT]",row["Mutation"])
        if r:
            converted_mutations[(row["Gene"],row["Mutation"])] = row["Mutation"]

        r = re.search("c.[0-9]+dup[ACGT]+",row["Mutation"])
        if r:
            converted_mutations[(row["Gene"],row["Mutation"])] = row["Mutation"]
        r = re.search("c.[0-9]+_[0-9]+dup[ACGT]+",row["Mutation"])
        if r:
            converted_mutations[(row["Gene"],row["Mutation"])] = row["Mutation"]
        



        r = re.search("c.([0-9]+)del",row["Mutation"])
        if r:
            # "ethA" "c.341del"
            del_start = int(r.group(1))
            del_end = int(r.group(1))
            if gene.strand == "+":
                # rpoB "c.1282_1290del"
                genome_start = gene.start + del_start - 2
                genome_end = gene.start + del_end 
            else:
                # "ethA" "c.1057_1059del"
                genome_start = gene.start - del_end
                genome_end = gene.start - del_start + 2
            ref = refseq["Chromosome"][genome_start-1:genome_end-1]
            alt = ref[0]
            mutations[(row["Gene"],row["Mutation"])] = {"pos":genome_start, "ref":ref, "alt":alt,"gene":row["Gene"],"type":"nucleotide"}

        r = re.search("c.([0-9]+)_([0-9]+)del",row["Mutation"])
        if r:
            del_start = int(r.group(1))
            del_end = int(r.group(2))
            if gene.strand == "+":
                # rpoB "c.1282_1290del"
                genome_start = gene.start + del_start - 2
                genome_end = gene.start + del_end 
            else:
                # "ethA" "c.1057_1059del"
                genome_start = gene.start - del_end
                genome_end = gene.start - del_start + 2
            ref = refseq["Chromosome"][genome_start-1:genome_end-1]
            alt = ref[0]
            mutations[(row["Gene"],row["Mutation"])] = {"pos":genome_start, "ref":ref, "alt":alt,"gene":row["Gene"],"type":"nucleotide"}

        r = re.search("c.-([0-9]+)del",row["Mutation"])
        if r:
            del_start = int(r.group(1))
            del_end = int(r.group(1))
            if gene.strand == "+":
               # "embA" "c.-29_-28del"
                genome_start = gene.start + del_start - 1
                genome_end = gene.start + del_end + 1
            else:
                # "alr" "c.-283_-280delCAAT"
                genome_start = gene.start - del_end - 1
                genome_end = gene.start - del_start + 1
            ref = refseq["Chromosome"][genome_start-1:genome_end-1]
            alt = ref[0]
            mutations[(row["Gene"],row["Mutation"])] = {"pos":genome_start, "ref":ref, "alt":alt,"gene":row["Gene"],"type":"nucleotide"}

        
        r = re.search("c.(-[0-9]+)_(-[0-9]+)del",row["Mutation"])
        if r:
            del_start = int(r.group(1))
            del_end = int(r.group(2))
            if gene.strand == "+":
               # "embA" "c.-29_-28del"
                genome_start = gene.start + del_start - 1
                genome_end = gene.start + del_end + 1
            else:
                # "alr" "c.-283_-280delCAAT"
                genome_start = gene.start - del_end - 1
                genome_end = gene.start - del_start + 1
            ref = refseq["Chromosome"][genome_start-1:genome_end-1]
            alt = ref[0]
            mutations[(row["Gene"],row["Mutation"])] = {"pos":genome_start, "ref":ref, "alt":alt,"gene":row["Gene"],"type":"nucleotide"}

        

        r = re.search("c.([0-9]+)_([0-9]+)ins([ACGT]+)", row["Mutation"])
        if r:
            ins_start = int(r.group(1))
            ins_end = int(r.group(2))
            ins_seq = r.group(3)
            if gene.strand == "+":
                # "rpoB" "c.1296_1297insTTC"
                genome_start = gene.start + ins_start - 1 
                genome_end = gene.start + ins_end - 1
            else:
                # "pncA" "c.521_522insT"
                ins_seq = pp.revcom(ins_seq)
                genome_start = gene.start - ins_start 
                genome_end = gene.start - ins_end + 2

            ref = refseq["Chromosome"][genome_start-1:genome_end-1]
            alt = ref + ins_seq
            mutations[(row["Gene"],row["Mutation"])] = {"pos":genome_start, "ref":ref, "alt":alt,"gene":row["Gene"],"type":"nucleotide"}
        

        if row["Mutation"] == "frameshift":
            converted_mutations[(row["Gene"],row["Mutation"])] = row["Mutation"]
        if row["Mutation"] == "large_deletion":
            converted_mutations[(row["Gene"],row["Mutation"])] = row["Mutation"]
        if row["Mutation"][:19] == "any_missense_codon_":
            converted_mutations[(row["Gene"],row["Mutation"])] = row["Mutation"]
        
    print("Converting %s mutations" % len(mutations))
    mutation_conversion = get_ann(mutations)
    for key in mutation_conversion:
        converted_mutations[key] = mutation_conversion[key]

    with open(args.out+".csv","w")  as O:
        with open(args.out+".log","w")  as L:
            writer = csv.DictWriter(O,fieldnames=list(row))
            writer.writeheader()
            for row in csv.DictReader(open(args.csv)):
                key = (row["Gene"],row["Mutation"])
                if row["Mutation"]!= converted_mutations[key]:
                        L.write(f'Recoded {row["Gene"]} {row["Mutation"]} as {converted_mutations[key]}\n')
                
                row["Mutation"] = converted_mutations[key]
                    
                writer.writerow(row)
        
    
            
        

parser = argparse.ArgumentParser(description='tbprofiler script',formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument('--csv',type=str,help='File with samples')
parser.add_argument('--out',type=str,help='File with samples')
parser.add_argument('--gff',type=str,help='File with samples')
parser.add_argument('--ref',type=str,help='File with samples')
parser.set_defaults(func=main)

args = parser.parse_args()
args.func(args)