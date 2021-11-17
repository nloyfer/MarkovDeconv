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
| Data Availability  | Tissue or Cell-type |  Samples | PMID |
|---|---|---|---|
| [GSE100262](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE100262) |  Bcells          |	3	|	[29326230](https://pubmed.ncbi.nlm.nih.gov/29326230/)
| [PRJEB14591](https://www.ebi.ac.uk/ena/browser/view/PRJEB14591)	|	Tcells	|	6	|	[28783152](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5912503/)
| [PRJNA229470](https://www.ebi.ac.uk/ena/browser/view/PRJNA229470) |  Cardiomyocyte	|	3	|	[25335909](https://pubmed.ncbi.nlm.nih.gov/25335909/)
| [PRJNA310298](https://www.ebi.ac.uk/ena/browser/view/PRJNA310298) |  Hepatocyte	|	3	|	[27380908](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4934005/)
| [GSE42836](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE42836) |  Cerebellum, Colon, Kidney, Intestine	|	4	|	[23995138](https://pubmed.ncbi.nlm.nih.gov/23995138/)
| [ENCODE](https://www.encodeproject.org/search/?type=Experiment&control_type!=*&related_series.@type=ReferenceEpigenome&replicates.library.biosample.donor.organism.scientific_name=Mus+musculus&assay_title=WGBS&limit=all) |  Cerebellum, Intestine, Kidney	|	9	|	[32728249](https://www.encodeproject.org)
| [GSE67386](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE67386) |  Mammary Epithelial	|	3	|	[25959817](https://pubmed.ncbi.nlm.nih.gov/25959817/)
| [PRJNA329552](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA329552) |  Hypothalamus	|	3	|	[28498846](https://pubmed.ncbi.nlm.nih.gov/28498846/)
| [PRJNA344551](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA344551) |  Lung endothelial	|	3	|	[29749927](https://pubmed.ncbi.nlm.nih.gov/29749927/)
| in-progress |  Bcell, CD4Tcell, CD8Tcell, Neutrophil, buffycoat      |	7	|	in-progress

