# MarkovDeconv tutorial
## Installation and configuration
First, make sure to have [`wgbstools`](https://github.com/nloyfer/wgbs_tools) installed and configure the `mm9` genome
```bash
git clone https://github.com/nloyfer/wgbs_tools.git
cd wgbs_tools
python setup.py
wgbstools init_genome mm9
```

It is recommended to add wgbstools to your $PATH, E.g,
```bash
export PATH=${PATH}:$PWD
```

Next, install `MarkovDeconv`
```bash
git clone https://github.com/nloyfer/MarkovDeconv.git
cd MarkovDeconv/counter/
make
cd ..
```

## All set. Let's begin
### Data and region
For this short tutorial, we will be demonstrating how MarkovDeconv works on a subset of mouse cardiomyocyte-specific methylation blocks. 

```bash
$ cd tutorial
$ head cardio_testmarkers.bed
```

These cardiomyocyte-specific methylation blocks are closely linked to the regulation of cardiomyocyte-specific gene functions. To identify similar regions, `wgbstools find_markers` command can be used to find cell-type specific methylation blocks for two or more groups of samples.

```bash
$ head cardio_testmarkers_info.txt | cut -f 1-3, 15-16 
```

These test markers were identified using the following publicly available WGBS data from healthy normal mouse tissues and cell-types. The fastq files were downloaded from The Sequence Read Archive [SRA](https://www.ncbi.nlm.nih.gov/sra), mapped to mm9 using [Bismark](https://github.com/FelixKrueger/Bismark), processed using [`wgbstools`](https://github.com/nloyfer/wgbs_tools) `bam2pat` and finally sliced to to the regions of these top 17 cardiomyocyte-specific methylation markers.  
| Data Availability  | Tissue or Cell-type |  Replicates | PMID |
|---|---|---|---|
| [GSE100262](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE100262) |  Bcells          |	3	|	[29326230](https://pubmed.ncbi.nlm.nih.gov/29326230/)
