import csv
import json
import re
from collections import defaultdict
import argparse
from shutil import copyfile
import os.path
import sys
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

def revcom(s):
        """Return reverse complement of a sequence"""
        def complement(s):
                        basecomplement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', 'N': 'N'}
                        letters = list(s)
                        letters = [basecomplement[base] for base in letters]
                        return ''.join(letters)
        return complement(s[::-1])

def write_gene_pos(infile,genes,outfile):
	O = open(outfile,"w")
	for l in open(infile):
		row = l.strip().split()
		rv,gene,chr_start,chr_end,gene_start,gene_end = [row[0],row[1]]+[int(row[i]) for i in range(2,6)]
		if rv in genes:
			y = 0
			for i,chr_pos in enumerate(range(chr_start,chr_end+1)):
				x = 1 if gene_start<gene_end else -1
				if gene_start+(x*i)==0:
					y = 1 if gene_start<gene_end else -1
				O.write("Chromosome\t%s\t%s\t%s\n" % (chr_pos,rv,gene_start+(x*i)+y))
	O.close()

def parse_mutation(mut,gene,fasta_dict,gene_info):
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
		return ["%s%s>%s%s" % (codon_num,ref_aa,codon_num,alt_aa)]
	# Stop codon
	re_obj = re.search("p.([A-Z][a-z][a-z])([0-9]+)(\*)",mut)
	if re_obj:
		ref_aa = aa_long2short[re_obj.group(1)]
		alt_aa = re_obj.group(3)
		codon_num = re_obj.group(2)
		return ["%s%s>%s%s" % (codon_num,ref_aa,codon_num,alt_aa)]
	# Deletion single base
	re_obj = re.search("c.([\-0-9]+)del",mut)
	if re_obj:
		gene_start_nt = int(re_obj.group(1))
		strand = gene_info[gene]["strand"]
		if strand=="-":
			chr_start_nt = gene_info[gene]["end"] + gene_info[gene]["gene_end"] - gene_start_nt + (1 if gene_info[gene]["gene_end"]<0 else 0)
		else:
			chr_start_nt = gene_info[gene]["start"] - gene_info[gene]["gene_start"] + gene_start_nt - (0 if gene_start_nt<0 else 1)
		seq = fasta_dict["Chromosome"][chr_start_nt-2:chr_start_nt]
		return ["%s%s>%s" % (chr_start_nt-1,seq,seq[0])]
	# Deletion multi base
	re_obj = re.search("c.([\-0-9]+)_([\-0-9]+)del",mut)
	if re_obj:
		gene_start_nt = int(re_obj.group(1))
		gene_end_nt = int(re_obj.group(2))
		del_len = gene_end_nt-gene_start_nt+1
		strand = gene_info[gene]["strand"]
		if strand=="-":
			chr_start_nt = gene_info[gene]["end"] + gene_info[gene]["gene_end"] - gene_start_nt - (del_len-1) + (1 if gene_info[gene]["gene_end"]<0 else 0)
		else:
			chr_start_nt = gene_info[gene]["start"] - gene_info[gene]["gene_start"] + gene_start_nt - (0 if gene_start_nt<0 else 1)
		chr_end_nt = chr_start_nt+del_len-1
		seq = fasta_dict["Chromosome"][chr_start_nt-2:chr_end_nt]
		return ["%s%s>%s" % (chr_start_nt-1,seq,seq[0])]
	# Insertion
	re_obj = re.search("c.([0-9]+)_([0-9]+)ins([A-Z]+)",mut)
	if re_obj:
		gene_start_nt = int(re_obj.group(1))
		seq_ins = re_obj.group(3)
		strand = gene_info[gene]["strand"]
		if strand=="-":
			chr_start_nt = gene_info[gene]["end"] + gene_info[gene]["gene_end"] - gene_start_nt + 1
			seq_ins = revcom(seq_ins)
		else:
			chr_start_nt = gene_info[gene]["start"] - gene_info[gene]["gene_start"] + gene_start_nt - 1
		seq_start = fasta_dict["Chromosome"][chr_start_nt-1]
		return ["%s%s>%s" % (chr_start_nt,seq_start,seq_start+seq_ins)]
	## Promoter Mutation
	## c.-16G>C
	re_obj = re.search("c.(\-[0-9]+)([A-Z])>([A-Z])",mut)
	if re_obj:
		nt_pos = re_obj.group(1)
		ref_nt = re_obj.group(2)
		alt_nt = re_obj.group(3)
		return ["%s%s>%s" % (nt_pos,ref_nt,alt_nt)]
	## ncRNA Mutation
	## r.514a>c
	re_obj = re.search("r.([0-9]+)([a-z]+)>([a-z]+)",mut)
	if re_obj:
		nt_pos = re_obj.group(1)
		ref_nt = re_obj.group(2)
		alt_nt = re_obj.group(3)
		return ["%s%s>%s" % (nt_pos,ref_nt.upper(),alt_nt.upper())]
	## frameshift
	## frameshift
	re_obj = re.search("frameshift",mut)
	if re_obj:
		return ["frameshift"]
	## premature_stop
	## premature_stop
	re_obj = re.search("premature_stop",mut)
	if re_obj:
		return ["premature_stop"]
	## codon range
	## any_missense_codon_425_452
	re_obj = re.search("any_missense_codon_([0-9]+)_([0-9]+)",mut)
	if re_obj:
		start = int(re_obj.group(1))
		end = int(re_obj.group(2))
		return ["any_missense_codon_%s" % i for i in range(start,end+1)]
	## Large deletion
	## "large_deletion"
	re_obj = re.search("large_deletion",mut)
	if re_obj:
		return ["large_deletion"]
	return []
#	sys.exit("%s is not a valid formatted mutation... Exiting!" % mut)
def write_bed(gene_dict,gene_info,outfile):
	O = open(outfile,"w")
	lines = []
	for gene in gene_dict:
		if gene not in gene_info:
			sys.stderr.write("%s not found in the 'gene_info' dictionary... Exiting!" % gene)
			quit()
		lines.append(["Chromosome",int(gene_info[gene]["start"]),int(gene_info[gene]["end"]),gene_info[gene]["locus_tag"],gene_info[gene]["gene"],",".join(gene_dict[gene])])
	for line in sorted(lines,key=lambda x: x[1]):
		line[1] = str(line[1])
		line[2] = str(line[2])
		O.write("%s\n" %"\t".join(line))
	O.close()

def load_gene_info(filename):
	gene_info = {}
	for l in open(filename):
		row = l.rstrip().split()
		strand = "-" if row[0][-1]=="c" else "+"
		gene_info[row[0]] = {"locus_tag":row[0],"gene":row[1],"start":int(row[2]),"end":int(row[3]),"gene_start":int(row[4]),"gene_end":int(row[5]),"strand":strand}
		gene_info[row[1]] = {"locus_tag":row[0],"gene":row[1],"start":int(row[2]),"end":int(row[3]),"gene_start":int(row[4]),"gene_end":int(row[5]),"strand":strand}
	return gene_info

def main(args):
	fasta_dict = fa2dict("genome.fasta")
	gene_info = load_gene_info("genes.txt")
	db = {}
	locus_tag_to_drug_dict = defaultdict(set)
	for row in csv.DictReader(open(args.csv)):
		print(row)
		locus_tag = gene_info[row["Gene"]]["locus_tag"]
		muts = parse_mutation(row["Mutation"],locus_tag,fasta_dict,gene_info)
		print(muts)
		#	sys.exit()
		for mut in muts:
			locus_tag_to_drug_dict[locus_tag].add(row["Drug"].lower())
			if locus_tag not in db:
				db[locus_tag] = {}
			if mut not in db[locus_tag]:
				db[locus_tag][mut] = {"drugs":{}}
			db[locus_tag][mut]["drugs"][row["Drug"].lower()] = {}
			annotation_columns = set(row.keys()) - set(["Gene","Mutation","Drug"])
			for col in annotation_columns:
				if row[col]=="":continue
				db[locus_tag][mut]["drugs"][row["Drug"].lower()][col.lower()] = row[col]
	conf_file = "%s.config.json" % args.prefix
	genome_file = "%s.fasta" % args.prefix
	gff_file = "%s.gff" % args.prefix
	ann_file = "%s.ann.txt" % args.prefix
	barcode_file = "%s.barcode.bed" % args.prefix
	bed_file = "%s.bed" % args.prefix
	json_file = "%s.dr.json" % args.prefix
	conf = {
		"gff": os.path.abspath(gff_file), "ref": os.path.abspath(genome_file),
		"ann": os.path.abspath(ann_file), "barcode": os.path.abspath(barcode_file),
		"bed": os.path.abspath(bed_file), "json_db": os.path.abspath(json_file)
	}

	copyfile("genome.fasta", genome_file)
	copyfile("genome.gff", gff_file)
	copyfile("barcode.bed", barcode_file)
	write_gene_pos("genes.txt",list(db.keys()),ann_file)
	write_bed(locus_tag_to_drug_dict,gene_info,bed_file)
	json.dump(db,open(json_file,"w"))
	json.dump(conf,open(conf_file,"w"))

parser = argparse.ArgumentParser(description='Script to generate the files required to run TBProfiler',formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument('--prefix','-p',default="tbdb",type=str,help='The input CSV file containing the mutations')
parser.add_argument('--csv','-c',default="tbdb.csv",type=str,help='The prefix for all output files')
parser.set_defaults(func=main)
args = parser.parse_args()
args.func(args)
