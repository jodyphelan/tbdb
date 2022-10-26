# TBDB: A repository for the TBProfiler library

This repository contains the scripts and data to generate all files required to run [tb-profiler](https://github.com/jodyphelan/TBProfiler/). 

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
4. **Literature** - Any literature which provides evidence for the mutation. Pubmed IDs (e.g. PMC3315572)or DOIs are recommended but in theory anything can be put here. Multiple entries can be separated with ";".

The first three columns must contain a value, however literature may remain empty. Additional columns may be added and will be built into the json library, and can be output in the tb-profiler results.

###### Mutation format
Mutations must follow the HGVS nomenclature. Information on this format can be found [here](http://varnomen.hgvs.org/). The following types of mutations are currently allowed:
* Amino acid substitutions. *Example: S450L in rpoB would be p.Ser450Leu*
* Deletions in genes. *Example: Deletion of nucleotide 758 in tlyA would be c.758del*
* Insertion in genes. *Example: Insertion of GT between nucleotide 1850 and 1851 in katG would be c.1850_1851insGT*
* SNPs in non-coding RNAs. *Example: A to G at position 1401 in rrs would be n.1401a>g*
* SNPs in gene promoters. *Example: A to G 7 bases 5' of the start codon in pncA c.-7A>G*

In addition, sequence ontology terms can also be used in place of the mutation. Supported sequence ontology terms can be found at http://pcingola.github.io/SnpEff/se_inputoutput/#effect-prediction-details. For example the following line is used to denote that any frameshift in _katG_ confers resistance.

| Gene | Mutation           | Drug      | Confers    | Interaction | Literature | WHO Confidence |
|------|--------------------|-----------|------------|-------------|------------|----------------|
| katG | frameshift_variant | isoniazid | resistance |             |            |                |

Additionally, ranges can be applied to sequence ontology terms to limit the resistance association of any particular term to a certain region within the gene or protein or non-coding RNA. To use gene coordinates add "_c.X_Y" where *X* and *Y* are the ranges between which the variant should occur in. Similarly, for protein you can use "_p.X_Y" and for RNA use "_n.X_Y". For example to use any missense variant between codon 430 and 470, the term will be *missense_variant_p.450_470*.

**Important! The mutations and resulting library files are in reference to the H37Rv (NC_000962.3/AL123456.3) reference genome**

## Watchlist

There are some genes it may be of interest to record mutations even if we do not have any specific associated mutaitons. To allow this funcitonality we have included a "watchlist" file. To include genes just add them and the associated drug(s) to the `tbdb.watchlist.csv` file.

## TB-Profiler

This repo contains all the files required to generate a library for `tb-profiler`. To find out more about how to build and load the library please visit the [tb-profiler repo](https://github.com/jodyphelan/TBProfiler)
