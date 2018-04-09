

## MHap

#### Maximum likelihood inference of methylation haplotypes<br/><br/>
Download MHap [here]  
Download the user guide [here](https://github.com/bzmartinelli/MHap/blob/master/MHap_user_guide.pdf)<br/><br/>

MHap is used to infer methylation haplotypes and to estimate methylation haplotype frequencies as well as methylation entropy from DNA methylation data derived from short-read sequencing.

Genomic regions of interested are given by the user. MHap uses a sliding window to select windows with n number of CpG with intermediate methylation. For each genomic region, a set of windows is retrieved and the analysis is performed in each one of them.

Making use of the Expectation-Maximization (EM) algorithm, MHap identifies all methylation haplotypes consistent with the sequence reads and estimates the frequency of each haplotype within the windows. Once the algorithm has converged, the inferred haplotype frequencies are used to calculate the Shannon entropy. Also the methylation fraction of individual CpGs is calculated using the methylation states provided by the sequencing reads and using the haplotype frequencies.

#### Usage
`python PATH/MHap.py -gr <PATH/genomic_regions_file> -data <PATH/cpg_reads_file> 
-data_from <bissnp or bismark> [options] `


