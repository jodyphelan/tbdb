import csv
import json
import re
from collections import defaultdict
import argparse
from shutil import copyfile

def fa2dict(filename):
	fa_dict = {}
	seq_name = ""
	for l in open(filename):
		line = l.rstrip()
		if line[0] == ">":
			seq_name = line[1:].split()[0]
			fa_dict[seq_name] = []
		else:
			fa_dict[seq_name].append(line)
	result = {}
	for seq in fa_dict:
		result[seq] = "".join(fa_dict[seq])
	return result

def parse_mutation(mut,gene,fasta_dict):
	aa_long2short = {
	"Ala":"A","Arg":"R","Asn":"N","Asp":"D","Cys":"C",
	"Gln":"Q","Glu":"E","Gly":"G","His":"H","Ile":"I",
	"Leu":"L","Lys":"K","Met":"M","Phe":"F","Pro":"P",
	"Ser":"S","Thr":"T","Trp":"W","Tyr":"Y","Val":"V",
	"Stop":"*", "-":"-"
	}
	# AA change
	re_obj = re.search("p.([A-Z][a-z][a-z])([0-9]+)([A-Z][a-z][a-z])",mut)
	if re_obj:
		ref_aa = aa_long2short[re_obj.group(1)]
		alt_aa = aa_long2short[re_obj.group(3)]
		codon_num = re_obj.group(2)
		return "%s%s>%s%s" % (codon_num,ref_aa,codon_num,alt_aa)
	# Stop codon
	re_obj = re.search("p.([A-Z][a-z][a-z])([0-9]+)(\*)",mut)
	if re_obj:
		ref_aa = aa_long2short[re_obj.group(1)]
		alt_aa = re_obj.group(3)
		codon_num = re_obj.group(2)
		return "%s%s>%s%s" % (codon_num,ref_aa,codon_num,alt_aa)
	# Deletion
	re_obj = re.search("c.([0-9]+)_([0-9]+)del",mut)
	if re_obj:
		start_nt = int(re_obj.group(1))
		end_nt = int(re_obj.group(2))
		seq = fasta_dict[gene][start_nt-2:end_nt]
		return "%s%s>%s" % (start_nt-1,seq,seq[0])
	# Insertion
	re_obj = re.search("c.([0-9]+)_([0-9]+)ins([A-Z]+)",mut)
	if re_obj:
		start_nt = int(re_obj.group(1))
		seq_ins = re_obj.group(3)
		seq_start = fasta_dict[gene][start_nt-1]
		return "%s%s>%s" % (start_nt,seq_start,seq_start+seq_ins)
	## Promoter Mutation
	## c.-16G>C
	re_obj = re.search("c.(\-[0-9]+)([A-Z])>([A-Z])",mut)
	if re_obj:
		nt_pos = re_obj.group(1)
		ref_nt = re_obj.group(2)
		alt_nt = re_obj.group(2)
		return "%s%s>%s" % (nt_pos,ref_nt,alt_nt)
	## ncRNA Mutation
	## r.514a>c
	re_obj = re.search("r.([0-9]+)([a-z])>([a-z])",mut)
	if re_obj:
		nt_pos = re_obj.group(1)
		ref_nt = re_obj.group(2)
		alt_nt = re_obj.group(3)
		return "%s%s>%s" % (nt_pos,ref_nt.upper(),alt_nt.upper())
def write_bed(gene_dict,gene_info,outfile):
	O = open(outfile,"w")
	for gene in gene_dict:
		if gene not in gene_info:
			sys.stderr.write("%s not found in the 'gene_info' dictionary... Exiting!" % gene)
			quit()
		O.write("Chromosome\t%s\t%s\t%s\t%s\t%s\n" % (gene_info[gene]["start"],gene_info[gene]["end"],gene_info[gene]["locus_tag"],gene_info[gene]["gene"],",".join(gene_dict[gene])))
	O.close()

def load_gene_info(filename):
	gene_info = {}
	for l in open(filename):
		row = l.rstrip().split()
		gene_info[row[0]] = {"locus_tag":row[0],"gene":row[1],"start":row[2],"end":row[3],"strand":row[4]}
		gene_info[row[1]] = {"locus_tag":row[0],"gene":row[1],"start":row[2],"end":row[3],"strand":row[4]}
	return gene_info

def main(args):
	fasta_dict = fa2dict("genes.fasta")
	gene_info = load_gene_info("genes.txt")
	db = {}
	locus_tag_to_drug_dict = defaultdict(set)
	for row in csv.DictReader(open(args.csv)):
		locus_tag = gene_info[row["Gene"]]["locus_tag"]
		mut = parse_mutation(row["Mutation"],locus_tag,fasta_dict)
		locus_tag_to_drug_dict[locus_tag].add(row["Drug"].lower())
		if locus_tag not in db:
			db[locus_tag] = {}
		if mut not in db[locus_tag]:
			db[locus_tag][mut] = {"drugs":[]}
		db[locus_tag][mut]["drugs"].append(row["Drug"].lower())
		annotation_columns = set(row.keys()) - set(["Gene","Mutation","Drug"])
		for col in annotation_columns:
			if col not in db[locus_tag][mut]:
				db[locus_tag][mut][col] = [row[col]]
			else:
				db[locus_tag][mut][col].append(row[col])

	copyfile("genome.fasta", "%s.fasta" % args.prefix)
	copyfile("genome.gff", "%s.gff" % args.prefix)
	copyfile("ann.txt", "%s.ann.txt" % args.prefix)
	copyfile("barcode.bed", "%s.barcode.bed" % args.prefix)
	write_bed(locus_tag_to_drug_dict,gene_info,"%s.bed"%args.prefix)
	json.dump(db,open("%s.dr.json"%args.prefix,"w"))

parser = argparse.ArgumentParser(description='TBProfiler pipeline',formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument('--prefix','-p',default="tbdb",type=str,help='NGS Platform')
parser.add_argument('--csv','-c',default="tbdb.csv",type=str,help='NGS Platform')
parser.set_defaults(func=main)
args = parser.parse_args()
args.func(args)
