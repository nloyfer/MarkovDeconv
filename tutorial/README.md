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
#chr    start   end     startCpG        endCpG
chr2	76783821	76784035	1280271	1280278
chr3	138213689	138214595	2452558	2452576
chr3	153604387	153605925	2544184	2544201
chr4	129989114	129990365	3140758	3140779
chr5	72836518	72837756	3682357	3682371
chr5	113324207	113325325	3901749	3901768
chr6	24547545	24548897	4313192	4313212
chr9	66880594	66881081	6745430	6745449
chr9	71459209	71459240	6776082	6776085
chr12	53981538	53981815	8779720	8779737
```

These cardiomyocyte-specific methylation blocks are closely linked to the regulation of cardiomyocyte-specific gene functions. To identify similar regions, `wgbstools find_markers` command can be used to find cell-type specific methylation blocks for two or more groups of samples.

```bash
$ head cardio_testmarkers_info.txt | cut -f 1-3, 15-16 
chr	start	end	Gene	Function
chr13	12423161	12423293	Actn2	Cardiomyocyte differentiation;Cardiogenesis
chr3	138213689	138214595	Eif4e	Morphology of the heart; Cardiac hypertrophy
chr19	40571151	40571988	Sorbs1	Actomyosin structure organization
chr4	129989114	129990365	Fabp3	Morphology of the heart; Cardiac hypertrophy
chr9	66880594	66881081	Tpm1	Function of cardiac muscle; Cardiomyocyte differentiation
chr2	76783821	76784035	Ttn	Morphology of cardiomyocytes
chr3	153604387	153605925	Acadm	Cardiomyocyte differentiation;Cardiogenesis
chr15	27406708	27407121	Ank	Morphology of the heart; Cardiac hypertrophy
chr6	24547545	24548897	Lmod2	Myofibril assembly; Cardiogenesis
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

For this tutorial, the reference data sliced to to these test regions, can be found at:
```bash
$ cd Data/Train/
$ ls -lh *.pat.gz | head -5
cardiosplt1_95_cardiotestmarkers.pat.gz
cardiosplt2_95_cardiotestmarkers.pat.gz
cardiosplt3_95_cardiotestmarkers.pat.gz
hon_cerebellum_cardiotestmarkers.pat.gz
hon_colon_cardiotestmarkers.pat.gz
hon_intestine_cardiotestmarkers.pat.gz
hon_kidney_cardiotestmarkers.pat.gz
mCD19B__cardiotestmarkers.pat.gz
mCD4T__cardiotestmarkers.pat.gz
mCD8T__cardiotestmarkers.pat.gz
```

### Generate Mixin Test Data
Reference mouse cardiomyocyte WGBS data was split 0.95 train and 0.05 test. Reads from the 0.05 cardiomyocyte split were in-silico mixed into a background of reads from mouse buffy coat (or lymphocyte) WGBS datasets using [`wgbstools`](https://github.com/nloyfer/wgbs_tools) `mix_pat`. We performed three replicates for each admixture ratio assessed (0.1%, 0.5%, 1%, 2%, 5%, 10%, 15%, 50%).
```bash
$ cd Data/Test/Mixin/
$ ls -lh *.pat.gz | head -5
mcardiotest01mixWBC_1.pat.gz
mcardiotest01mixWBC_2.pat.gz
mcardiotest01mixWBC_3.pat.gz
mcardiotest05mixWBC_1.pat.gz
mcardiotest05mixWBC_2.pat.gz
```

### Visualization 
For example, the cardiomyocyte-specific methylation pattern at chr13:12423161-12423293 (Actn2) is hypomethylated in cardiomyocyte
$wgbstools vis --genome mm9 -r chr13:12423161-12423293 Data/Train/cardiosplt1_95_cardiotestmarkers.pat.gz --min_len 3 --yebl

<!--![alt text](Images/cardiomyocyte.png "Cardiomyocyte_Actn2")-->
<img src="Images/cardiomyocyte.png" width="600" height="600" />

