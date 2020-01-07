# TBDB: A repository for the TBProfiler library

This repository contains the scripts and data to generate all files required to run [TBProfiler](https://github.com/jodyphelan/TBProfiler/).

## Why is there a seperate github repository?

With analysis pipelines pretty much standardised, it is evident that accuracy of prediction is affected mostly by the underlying library of mutations. As new evidence for the inclusion or exclusion of mutations is generated, there is a constant need to  update and re-evaluate the mutation library. Moreover, it is important for the control of the library to be put in the hands of the end-users. By hosting the library on a separate repository (rather than buried deep in the profiling tool code) it makes it easier to find out exactly which mutations are present. Additionally, github has a number of useful features which can be utilised:
 - All changes to the library are tracked and can be easily be reverted to previous versions.
 - Users can raise concerns or discuss the library using the [Issues](https://github.com/jodyphelan/tbdb/issues) section of the github repository.
 - Different versions of the library can be maintained using [Forks](https://help.github.com/en/articles/fork-a-repo). Allowing users to experiment with the library without affecting the main project. These changes can then be merged into the main repo after the changes are reviewed.
 - Multiple users/developers can contribute towards the library.

**tl;dr** - Hosting it here makes it easier to update the library.

## Want to contribute?

If you think a mutation should be removed or added please raise and issue [here](https://github.com/jodyphelan/tbdb/issues).
If you want to help curate the library, leave a comment [here](https://github.com/jodyphelan/tbdb/issues/4).

#### Adding/removing mutations
Mutations can be added by submitting a pull request on a branch modified `tbdb.csv` file. If that previous sentence made no sense to you then you can suggest a change using an [issue](https://github.com/jodyphelan/tbdb/issues) and we will try help. On submitting a pull request the `tbdb_bot` will automatically calculate the confidence of the mutations in question and submit the results as a comment on the pull request (like [this](https://github.com/jodyphelan/tbdb/pull/5)). All `tbdb_bot` checks should pass, at least two reviews should be requested and upon review can be merged into the master branch

## How does it work?

The mutations are listed in [tbdb.csv](https://github.com/jodyphelan/tbdb/blob/master/tbdb.csv). These are parsed by `parse_db.py` to generate the json formatted database used by TBProfiler along with a few more files. Mutations can be removed and added from tbdb.csv and a new library can be built using `parse_db.py`.

#### tbdb.csv
This is a CSV file which must contain the following column headings:
1. **Gene** - These can be the gene names (e.g. *rpoB*) or locus tag (e.g. *Rv0667*).
2. **Mutation** - These must follow the [hgvs nomenclature](http://varnomen.hgvs.org/). More info down below.
3. **Drug** - Name of the drug
4. **Literature** - Any literature which provides evidence for the mutation. I would recommend using pubmed IDs (e.g. PMC3315572) but in theory anything can be put here. Multiple entries can be separated with ";".

The first three columns must contain a value, however literature may remain empty. Additional columns may be added and will be built into the json library, and can be output in the tb-profiler results.

###### Mutation format
Mutations must follow the HGVS nomenclature. Information on this format can be found [here](http://varnomen.hgvs.org/). The following types of mutations are currently allowed:
* Amino acid substitutions. *Example: S450L in rpoB would be p.Ser450Leu*
* Deletions in genes. *Example: Deletion of nucleotide 758 in tlyA would be c.758del*
* Insertion in genes. *Example: Insertion of GT between nucleotide 1850 and 1851 in katG would be c.1850_1851insGT*
* SNPs in non-coding RNAs. *Example: A to G at position 1401 in rrs would be r.1401a>g*
* SNPs in gene promoters. *Example: A to G 7 bases 5' of the start codon in pncA c.-7A>G*

**Important! The mutations and resulting library files are in reference to the H37Rv (NC_000962.3/AL123456.3) reference genome**

## Confidence values

Confidence values for mutations will be included upon building the database. By default these are source from the [tbdb.confidence.csv](https://github.com/jodyphelan/tbdb/blob/master/tbdb.confidence.csv) file. These confidence values have been precomputed using a database of >16,000 isolates with phenotypic and genotypic data available. The database was created as specified in the [tb-profiler paper](https://genomemedicine.biomedcentral.com/articles/10.1186/s13073-019-0650-x). The confidence values were calculated using a method similar to that as proposed by the ReSeqTB consortium and is detailed in supplementary material 5 of [this paper](https://erj.ersjournals.com/content/50/6/1701354). In short, the risk ratio is used to infer if mutations are positively or negatively associated with resistance and the odds ratio is used to characterise the effect size. Associated p-value for these metrics are calculated and if any fall below 0.05, the confidence is graded as _indeterminate_. If the p-values are significant the OR is used to grade the mutations in high (OR>10), moderate (5<OR<=10), low (1<OR<=5) or no_association (OR<=1). The script used to download the necessary data and calculate the metrics can be found [here](https://github.com/jodyphelan/tbdb/blob/master/scripts/generate_confidence.py).

### Can I use my own grading system?

Of course! If you would prefer to assign your own confidence values you can supply a CSV file with the columns drug, gene, mutation and confidence to the `--confidence` flag when using `parse_db.py`. If a mutation if found in the confidence list it will be included, if not it will be set to indeterminate.

## Generating a new library

##### Install

Just download the repository using `https://github.com/jodyphelan/tbdb.git`. The script is written in pure python with no dependencies.

##### Run
Once the tbdb.csv file is finalised the `parse_db.py` script can be run.
When run without any optional arguments like this: `python parse_db.py`, it will by defualt read the **tbdb.csv** file and produce all output with "tbdb" as the file prefix. These can be changed using the `--csv` flag to specify a different CSV file and `--prefix` to specify the output file prefixes. Use the `--help` flag to print the script help:

```
usage: parse_db.py [-h] [--prefix PREFIX] [--csv CSV]

Script to generate the files required to run TBProfiler

optional arguments:
  -h, --help            show this help message and exit
  --prefix PREFIX, -p PREFIX
                        The input CSV file containing the mutations (default:
                        tbdb)
  --csv CSV, -c CSV     The prefix for all output files (default: tbdb.csv)
```

## Output files

After running `parse_db.py` you will have 6 output files (considering you did not change the output files prefix):
1. tbdb.fasta - The reference genome used for mapping
2. tbdb.gff - The reference genome annotation
3. tbdb.barcode.bed - The lineage specific marked used to assign lineages
4. tbdb.bed - A bed file containing all candidate genes and the associated drugs
5. tbdb.ann.txt - A file listing the genomic and corresponding gene positions
6. tbdb.config.json - A json file listing the absolute path of all the above files.

The **tbdb.config.json** file can then be loaded by TBProfiler using `tb-profiler load_library drdb.config.json` and is then ready to be used.

## Using alternate reference names

The tbdb database will assume you have mapped your data to a reference with "Chromosome" as the sequence name. If your reference sequence is the same but has a differenct name e.g "NC_000962.3". You can generate a database with an alternate sequence name using the `--seqname` flag.

## Watchlist

There are some genes it may be of interest to record mutations even if we do not have any specific associated mutaitons. To allow this funcitonality we have included a "watchlist" file. To include genes just add them and the associated drug(s) to the `tbdb.watchlist.csv` file.
## Future plans

- Build in customisation of the lineage assigning SNPs
