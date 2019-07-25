import csv
from collections import defaultdict
import re

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

def load_gene_info(filename):
	results = {}
	for l in open(filename):
		row = l.strip().split()
		rv,gene,chr_start,chr_end,gene_start,gene_end = [row[0],row[1]]+[int(row[i]) for i in range(2,6)]
		results[gene] = {}

		y = 0
		for i,chr_pos in enumerate(range(chr_start,chr_end+1)):
			x = 1 if gene_start<gene_end else -1
			if gene_start+(x*i)==0:
				y = 1 if gene_start<gene_end else -1
			gene_pos = gene_start+(x*i)+y
			results[gene][gene_pos] = chr_pos
		results[rv] = results[gene]
	return results

def revcom(s):
	"""Return reverse complement of a sequence"""
	def complement(s):
			basecomplement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', 'N': 'N'}
			letters = list(s)
			letters = [basecomplement[base] for base in letters]
			return ''.join(letters)
	return complement(s[::-1])


aa_long2short = {"Ala":"A","Arg":"R","Asn":"N","Asp":"D","Cys":"C","Gln":"Q","Glu":"E","Gly":"G","His":"H","Ile":"I","Leu":"L","Lys":"K","Met":"M","Phe":"F","Pro":"P","Ser":"S","Thr":"T","Trp":"W","Tyr":"Y","Val":"V","*":"*"}

codon2aa = {
    'ATA':'I', 'ATC':'I', 'ATT':'I', 'ATG':'M',
    'ACA':'T', 'ACC':'T', 'ACG':'T', 'ACT':'T',
    'AAC':'N', 'AAT':'N', 'AAA':'K', 'AAG':'K',
    'AGC':'S', 'AGT':'S', 'AGA':'R', 'AGG':'R',
    'CTA':'L', 'CTC':'L', 'CTG':'L', 'CTT':'L',
    'CCA':'P', 'CCC':'P', 'CCG':'P', 'CCT':'P',
    'CAC':'H', 'CAT':'H', 'CAA':'Q', 'CAG':'Q',
    'CGA':'R', 'CGC':'R', 'CGG':'R', 'CGT':'R',
    'GTA':'V', 'GTC':'V', 'GTG':'V', 'GTT':'V',
    'GCA':'A', 'GCC':'A', 'GCG':'A', 'GCT':'A',
    'GAC':'D', 'GAT':'D', 'GAA':'E', 'GAG':'E',
    'GGA':'G', 'GGC':'G', 'GGG':'G', 'GGT':'G',
    'TCA':'S', 'TCC':'S', 'TCG':'S', 'TCT':'S',
    'TTC':'F', 'TTT':'F', 'TTA':'L', 'TTG':'L',
    'TAC':'Y', 'TAT':'Y', 'TAA':'*', 'TAG':'*',
    'TGC':'C', 'TGT':'C', 'TGA':'*', 'TGG':'W'}

aa2codon = defaultdict(list)
for codon in codon2aa:
	aa2codon[codon2aa[codon]].append(codon)

gene_info = load_gene_info("genes.txt")
gene_strand = {}
for l in open("genes.txt"):
	row = l.rstrip().split()
	rv,gene = row[:2]
	gene_strand[rv] = "-" if rv[-1]=="c" else "+"
	gene_strand[gene] = "-" if rv[-1]=="c" else "+"

ref_fa = fa2dict("genome.fasta")
for row in csv.DictReader(open("tbdb.csv")):
	mut = row["Mutation"]
	strand = gene_strand[row["Gene"]]
	# print("-"*40)
	# print(row)
	if mut[0]=="p":
		re_obj = re.search("p.([A-Za-z]+)([0-9]+)([A-Za-z\*]+)",mut)
		ref_aa = aa_long2short[re_obj.group(1)]
		codon_num = int(re_obj.group(2))
		alt_aa = aa_long2short[re_obj.group(3)]

		if strand=="-":
			gene_pos = [(codon_num*3),(codon_num*3)-1,(codon_num*3)-2]
		else:
			gene_pos = [(codon_num*3)-2,(codon_num*3)-1,codon_num*3]
		chr_pos = [gene_info[row["Gene"]][i] for i in gene_pos]
		ref_codon = "".join([ref_fa["Chromosome"][i-1] for i in chr_pos])
		if strand=="-":
			ref_codon = revcom(ref_codon)
		if ref_aa!=codon2aa[ref_codon]:
			exit("Conversion wrong")
		alt_codons = aa2codon[alt_aa]
		for alt_codon in alt_codons:
			changed_codon_positions = [i for i in [0,1,2] if ref_codon[i]!=alt_codon[i]]
			if strand=="+":
				changed_chr_positions = [chr_pos[i] for i in changed_codon_positions]
				ref_nucleotides = [ref_codon[i] for i in changed_codon_positions]
				alt_nucleotides = [alt_codon[i] for i in changed_codon_positions]
			else:
				changed_chr_positions = [chr_pos[::-1][i] for i in changed_codon_positions]
				ref_nucleotides = [ref_codon[i] for i in changed_codon_positions]
				alt_nucleotides = [alt_codon[i] for i in changed_codon_positions]
			print("%s\t%s\t%s\t%s\t%s\t%s" % (row["Gene"],row["Mutation"],row["Drug"],",".join([str(x) for x in changed_chr_positions]),",".join(ref_nucleotides),",".join(alt_nucleotides)))
	elif mut[0]=="c" and "ins" not in mut and "del" not in mut:
		re_obj = re.search("c.(-[0-9]+)([A-Z])>([A-Z])",mut)
		ref_nuc = re_obj.group(2)
		gene_pos = int(re_obj.group(1))
		alt_nuc = re_obj.group(3)
		chr_pos = gene_info[row["Gene"]][gene_pos]
		if strand=="-":
			alt_nuc = revcom(alt_nuc)
			ref_nuc = revcom(ref_nuc)
		print("%s\t%s\t%s\t%s\t%s\t%s" % (row["Gene"],row["Mutation"],row["Drug"],chr_pos,ref_nuc,alt_nuc))
	elif mut[0]=="r":
		re_obj = re.search("r.([0-9]+)([a-z])>([a-z])",mut)
		ref_nuc = re_obj.group(2).upper()
		gene_pos = int(re_obj.group(1))
		alt_nuc = re_obj.group(3).upper()
		chr_pos = gene_info[row["Gene"]][gene_pos]
		print("%s\t%s\t%s\t%s\t%s\t%s" % (row["Gene"],row["Mutation"],row["Drug"],chr_pos,ref_nuc,alt_nuc))
