# FindMarkers_Human README

This command takes as input:
1. group file: a csv table\ text file with three columns "name","sname","group" where name = name of file minus file extension, sname = sample name (or same as name), group = cell-type/tissue represented.

2. lbeta files: a binary file with a fixed size of 2 * 8 * NR_SITES bytes. Can be produced from BAM files using wgbstools bam2pat and then wgbstools beta_to_blocks --lbeta. 
For each of the NR_SITES CpG sites in the genome, it holds 2 values: the #meth and #covered. Meaning, the i'th row in the matrix corresponds to the i'th CpG site:
#meth: the number of times the i'th site (CpGi) site is observed in a methylated state.
#coverage: the total number of times i'th site (CpGi) is observed. #coverage==0 is equivalent to a missing value (NaN).
CpGi's beta value is obtained by dividing #meth / #coverage.

3.blocks file: a bed file with 2 extra columns for CpG indexes. Could be the output of the wgbstools segment command, or any custom bed file once you added the [startCpG, endCpG] columns with wgbstools convert -L BED_FILE.


## Usage Example
```bash
Rscript find_markers.R -g groups.csv --bins_dir lbeta/ -o OUTPUT_DIR/ -b blocks.bed.gz
```
