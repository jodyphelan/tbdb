# TBDB: A repository for the TBProfiler library

This repository contains the scripts and data to generate all files required to run [TBProfiler2](https://github.com/jodyphelan/TBProfiler2/).

## Why is there a seperate github repository?

With analysis pipelines pretty much standardised, it is evident that accuracy of prediction is affected mostly by the underlying library of mutations. As new evidence for the inclusion or exclusion of mutations is generated, there is a constant need to  update and re-evaluate the mutation library. Moreover, it is important for the control of the library to be put in the hands of the end-users. By hosting the library on a separate repository (rather than buried deep in the profiling tool code) it makes it easier to find out exactly which mutations are present. Additionally, github has a number of useful features which can be utilised:
 - All changes to the library are tracked and can be easily be reverted to previous versions.
 - Users can raise concerns or discuss the library using the [Issues](https://github.com/jodyphelan/tbdb/issues) section of the github repository.
 - Different versions of the library can be maintained using [Forks](https://help.github.com/en/articles/fork-a-repo). Allowing users to experiment with the library without affecting the main project. These changes can then be merged into the main repo after the changes are reviewed.
 - Multiple users/developers can contribute towards the library.

**tl;dr** - Hosting it here makes it easier to update the library.

## How does it work?

The mutations are listed in [tbdb.csv](https://github.com/jodyphelan/tbdb/blob/master/tbdb.csv). These are parsed by `parse_db.py` to generate the json formatted database used by TBProfiler along with a few more files. Mutations can be removed and added from tbdb.csv and a new library can be built using `parse_db.py`.

#### tbdb.csv
This is a CSV file which must contain the following column headings:
1. **Gene** - These can be the gene names (e.g. *rpoB*) or locus tag (e.g. *Rv0667*).
2. **Mutation** - These must follow the [hgvs nomenclature](http://varnomen.hgvs.org/). More info down below.
3. **Drug** - Name of the drug
4. **Literature** - Any literature which provides evidence for the mutation. I would recommend using pubmed IDs (e.g. PMC3315572) but in theory anything can be put here. Multiple entries can be separated with ";".

The first three columns must contain a value, however literature may remain empty. Additional columns may be added and will be built into the json library. Currently, only the columns listed above are actually used by TBProfiler to generate the results but we are working on allowing a more flexible output in the future.

###### Mutation format
Mutations must follow the HGVS nomenclature. Information on this format can be found [here](http://varnomen.hgvs.org/). The following types of mutations are currently allowed:
* Amino acid substitutions. *Example: S450L in rpoB would be p.Ser450Leu*
* Deletions in genes. *Example: Deletion of nucleotide 758 in tlyA would be c.758del*
* Insertion in genes. *Example: Insertion of GT between nucleotide 1850 and 1851 in katG would be c.1850_1851insGT*
* SNPs in non-coding RNAs. *Example: A to G at position 1401 in rrs would be r.1401a>g*
* SNPs in gene promoters. *Example: A to G 7 bases 5' of the start codon in pncA c.-7A>G*

**Important! The mutations and resulting library files are in reference to the H37Rv reference genome**

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

The tbdb database will assume you have mapped your data to a reference with "Chromosome" as the sequenc name. If your reference sequence is the same but has a differenct name e.g "NC_000962.3". You can generate a database with an alternate sequence name using the `--seqname` flag.

## FAQ

Will populate this once we get some frequently asked questions!

## Future plans

- Build in customisation of the lineage assigning SNPs
