#! /usr/bin/env python
from __future__ import division
import json
from collections import defaultdict
import re
import argparse
import os
from tqdm import tqdm
import sys
import csv
import numpy as np
import statsmodels.api as sm
import math
import subprocess

def get_codon_number(x):
	re_obj = re.search("p.[A-Za-z]+([0-9]+)[A-Za-z]+",x)
	return re_obj.group(1)

def download(args):
	subprocess.call("wget http://pathogenseq.lshtm.ac.uk/downloads/tbprofiler_results.tgz",shell=True)
	subprocess.call("tar -xvf tbprofiler_results.tgz",shell=True)


def main(args):

	gene2locustag = {}
	drug2genes = defaultdict(set)
	for l in open(args.bed):
		row = l.rstrip().split()
		gene2locustag[row[4]] = row[3]
		for d in row[5].split(","):
			drug2genes[d].add(row[3])

	multi_change_codons = defaultdict(list)
	if args.add_to_db:
		for row in csv.DictReader(open(args.add_to_db)):
			if "any_missense_codon" in row["Mutation"]:
				multi_change_codons[gene2locustag[row["Gene"]]].append(row["Mutation"].replace("any_missense_codon_",""))

	meta = {}
	for row in csv.DictReader(open(args.meta)):
		meta[row["id"]] = row



	if args.samples:
		samples = [x.rstrip() for x in open(args.samples).readlines()]
	else:
		samples = [x.replace(".results.json","") for x in os.listdir("%s/" % args.dir) if x[-13:]==".results.json"]

	variants = defaultdict(lambda:defaultdict(list))
	mutation_types = defaultdict(dict)
	sys.stderr.write("Loading tb-profiler results\n")
	for s in tqdm(samples):
		tmp = json.load(open("%s/%s.results.json" % (args.dir,s)))
		for var in tmp["dr_variants"]:
			variants[var["locus_tag"]][var["change"]].append(s)
			if "large_deletion" in var["type"]:
				variants[var["locus_tag"]]["large_deletion"].append(s)
			elif  "frameshift" in var["type"]:
				variants[var["locus_tag"]]["frameshift"].append(s)
			if var["locus_tag"] in multi_change_codons and var["type"]=="missense":
				if get_codon_number(var["change"]) in multi_change_codons[var["locus_tag"]]:
					variants[var["locus_tag"]]["any_missense_codon_"+get_codon_number(var["change"])].append(s)
		for var in tmp["other_variants"]:
			variants[var["locus_tag"]][var["change"]].append(s)
			if "large_deletion" in var["type"]:
				variants[var["locus_tag"]]["large_deletion"].append(s)
			elif  "frameshift" in var["type"]:
				variants[var["locus_tag"]]["frameshift"].append(s)
			if var["locus_tag"] in multi_change_codons and var["type"]=="missense":
				if get_codon_number(var["change"]) in multi_change_codons[var["locus_tag"]]:
					variants[var["locus_tag"]]["any_missense_codon_"+get_codon_number(var["change"])].append(s)
	print("Collected %s unique variants in %s genes" % (sum([len(variants[x]) for x in variants]),len(variants)))
	results = []
	for drug in sorted(drug2genes):
		if drug not in meta[samples[0]]: continue
		for gene in sorted(drug2genes[drug]):
			if gene not in variants: continue
			sys.stderr.write("Calculating metrics for %s with %s\n" % (gene,drug))
			for mutation in tqdm(variants[gene]):
				result = {"gene":gene,"drug":drug,"mutation":mutation}
				t = [
						[0.5,0.5],
						[0.5,0.5]
					 ]
				for s in samples:
					if s not in meta: continue
					if meta[s][drug]=="1" and s in variants[gene][mutation]: 		t[0][0]+=1
					if meta[s][drug]=="0" and s in variants[gene][mutation]: 		t[0][1]+=1
					if meta[s][drug]=="1" and s not in variants[gene][mutation]: 	t[1][0]+=1
					if meta[s][drug]=="0" and s not in variants[gene][mutation]: 	t[1][1]+=1
				t2 = sm.stats.Table2x2(np.asarray(t))
				result["OR"] = t2.oddsratio
				result["OR_pval"] = t2.oddsratio_pvalue()
				result["RR"] = t2.riskratio
				result["RR_pval"] = t2.riskratio_pvalue()
				result["table"] = t
				results.append(result)

	OR_corrected_pvals = sm.stats.multipletests([r["OR_pval"] for r in results],method="fdr_bh")[1]
	RR_corrected_pvals = sm.stats.multipletests([r["RR_pval"] for r in results],method="fdr_bh")[1]
	sys.stderr.write("Writing results\n")
	organised_results = defaultdict(lambda:defaultdict(dict))
	with open(args.out,"w") as O:
		columns = ["drug","gene","mutation","table","OR","OR_pval","OR_pval_corrected","RR","RR_pval","RR_pval_corrected","confidence"]
		writer = csv.DictWriter(O,fieldnames=columns)
		writer.writeheader()
		for i in tqdm(range(len(results))):
			organised_results[results[i]["gene"]][results[i]["mutation"]][results[i]["drug"]] = results[i]
			results[i]["OR_pval_corrected"] = OR_corrected_pvals[i]
			results[i]["RR_pval_corrected"] = RR_corrected_pvals[i]
			if results[i]["OR"]>10 and results[i]["OR_pval_corrected"]<args.pval_cutoff and results[i]["RR"]>1 and results[i]["RR_pval_corrected"]<args.pval_cutoff:
				results[i]["confidence"] = "high"
			elif 5<results[i]["OR"]<=10 and results[i]["OR_pval_corrected"]<args.pval_cutoff and results[i]["RR"]>1 and results[i]["RR_pval_corrected"]<args.pval_cutoff:
				results[i]["confidence"] = "moderate"
			elif 1<results[i]["OR"]<=5 and results[i]["OR_pval_corrected"]<args.pval_cutoff and results[i]["RR"]>1 and results[i]["RR_pval_corrected"]<args.pval_cutoff:
				results[i]["confidence"] = "low"
			elif (results[i]["OR"]<=1 and results[i]["OR_pval_corrected"]<args.pval_cutoff) or (results[i]["RR"]<=1 and results[i]["RR_pval_corrected"]<args.pval_cutoff):
				results[i]["confidence"] = "no_association"
			else:
				results[i]["confidence"] = "indeterminate"
			writer.writerow(results[i])

	if args.add_to_db:
		with open(args.add_to_db.replace(".csv")+".confidence.csv","w") as O:
			reader = csv.DictReader(open(args.add_to_db))
			columns = reader.fieldnames + ["confidence"]
			writer = csv.DictWriter(O,fieldnames = columns)
			writer.writeheader()
			for row in reader:
				if row["Gene"] not in gene2locustag:
					sys.stderr.write(str(row)+"\n")
					continue
				locus_tag = gene2locustag[row["Gene"]]

				if locus_tag in organised_results and row["Mutation"] in organised_results[locus_tag] and row["Drug"] in organised_results[locus_tag][row["Mutation"]]:
					row["confidence"] = organised_results[locus_tag][row["Mutation"]][row["Drug"]]["confidence"]
				else:
					row["confidence"] = "indeterminate"
				writer.writerow(row)


parser = argparse.ArgumentParser(description='TBProfiler pipeline',formatter_class=argparse.ArgumentDefaultsHelpFormatter)
subparsers = parser.add_subparsers(help="Task to perform")

parser_sub = subparsers.add_parser('calculate', help='Collate results form multiple samples together', formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser_sub.add_argument('--meta',type=str,help="Meta data",required=True)
parser_sub.add_argument('--bed',type=str,help="BED file with genes",required=True)
parser_sub.add_argument('--out',type=str,help="Output file",required=True)
parser_sub.add_argument('--samples',type=str,help='NGS Platform')
parser_sub.add_argument('--dir',default="tbprofiler_results/",type=str,help='NGS Platform')
parser_sub.add_argument('--add-to-db',default=None,type=str,help='Choose db CSV file to add to')
parser_sub.add_argument('--pval-cutoff',default=0.05,type=float,help='Pvalue cutoff to use for the corrected OR and RR p-vaule significance')
parser_sub.set_defaults(func=main)

parser_sub = subparsers.add_parser('download', help='Download phenotypic and genotypic data for >16000 strains', formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser_sub.set_defaults(func=download)

args = parser.parse_args()
if vars(args)=={}:
	parser.print_help(sys.stderr)
else:
	args.func(args)
