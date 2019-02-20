import re
import sys

genes = {
	"non-coding":["rrs","rrl",],
	"coding":[
		"ahpC","embA","embB","embC","embR","ethA","ethR","folC",
		"gid","gyrA","gyrB","inhA","kasA","katG","panD","pncA","ribD","rplC",
		"rpoB","rpoC","rpsA","rpsL","Rv0678","thyA","tlyA",],
	"promoter":[
		"ahpC_promoter","embA_promoter","fabG1_promoter","katG_promoter","pncA_promoter",
		"eis_promoter",
	]
	}
promoter2gene = {"ahpC_promoter":"ahpC","embA_promoter":"embA","fabG1_promoter":"fabG1","katG_promoter":"katG","pncA_promoter":"pncA",
"eis_promoter":"eis",}
print("Gene\tMutation\tDrug")
for l in open("drdb.txt"):
	row = l.rstrip().split()
	if l[0]=="#":sys.stderr.write(l);continue
	if "N" in row[2] or "N" in row[3]: continue
	gene = row[4]
	drug = row[0].lower()
	if gene in genes["non-coding"]:
		re_obj = re.search("([A-Z])([0-9]+)([A-Z])",row[5])
		mut = "r.%s%s>%s" % (re_obj.group(2),re_obj.group(1).lower(),re_obj.group(3).lower())
	elif gene in genes["coding"]:
		if len(row[2])==len(row[3]): # substitution
			re_obj = re.search("([A-Z][a-z][a-z])([0-9]+)([A-Z][a-z][a-z])",row[5])
			if not re_obj: sys.exit("Strange substitution")
			ref = re_obj.group(1)
			alt = re_obj.group(3)
			pos = int(re_obj.group(2))
			if alt=="Sto": alt="*"
			mut = "p.%s%s%s" % (ref,pos,alt)
		elif len(row[2])>len(row[3]): # deletion
			re_obj = re.search("([A-Z]+)([0-9]+)([A-Z])",row[5])
			if not re_obj: sys.exit("Strange deletion")
			ref = re_obj.group(1)
			alt = re_obj.group(3)
			pos = int(re_obj.group(2))
			deletion_len = len(ref)-1
			if deletion_len==1:
				mut = "c.%sdel" % (pos+1)
			else:
				mut = "c.%s_%sdel" % (pos+1,pos-1+len(ref))
		elif len(row[2])<len(row[3]): # insertion
			re_obj = re.search("([A-Z])([0-9]+)([A-Z]+)",row[5])
			if not re_obj: sys.exit("Strange insertion")
			ref = re_obj.group(1)
			alt = re_obj.group(3)
			pos = int(re_obj.group(2))
			mut = "c.%s_%sins%s" % (pos,pos+1,alt[1:])
		else:
			sys.exit("Unknown coding mutation")
	elif gene in genes["promoter"]:
		re_obj = re.search("([A-Z])(\-[0-9]+)([A-Z])",row[5])
		if not re_obj: sys.exit("Strange promoter")
		ref = re_obj.group(1)
		pos = re_obj.group(2)
		alt = re_obj.group(3)
		mut = "c.%s%s>%s" % (pos,ref,alt)
		gene = promoter2gene[gene]
	else:
		sys.exit("Unknown Mutation\n%s" %l)
	print("%s\t%s\t%s" % (gene,mut,drug))
