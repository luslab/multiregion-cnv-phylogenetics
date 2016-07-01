# multiregion-cnv-phylogenetics
This respository details the analysis associated with the publication 'Intra-Tumor Genetic Heterogeneity in Wilms Tumor: Clonal Evolution and Clinical Implications'.

# Bioinformatics Supplementary Methodology
# Log R ratio genomic waves detection and correction
For each autosomal SNP array probe, the mean Log R ratio (log2[observed intensity/reference intensity], LRR) across all normal tissue arrays in our dataset (normal kidney and blood samples, n = 22) was calculated, generating a mean normal profile. For each autosome a moving average profile was generated from the mean normal profile using a window size of 10 probes. The moving average profile from the normal samples represents the genomic waves of each autosome in the array. Genomic wave detection was applied separately for male and female X chromosome LRRs. To centre the male X chromosome signal around zero, the profile was subtracted by the median value. 
LRR values were corrected for genomic waves in each array for all autosome and X chromosome probes. The contribution of the genomic waves to genome-wide LRR variation in each array is determined by subtracting the genomic waves (G) from corresponding array probe LRRs (S) using the coefficient (ɑ) which minimises LRR population variation across all genomic probes to generate normalised LRR (N),

arg min {var(Si - aGi)} = amin,

Si-aGi = Ni.

# Regions of copy number changes per case
The normalised LRRs from each array were segmented and copy number states were called using the ‘CGHcall’ R package[1] in Bioconductor[2] using the ‘sdundo’ option for the undo.splits parameter for segmentation (undo.SD=5, clen=13, relSDlong=8.33). 
Called copy number segment boundaries were smoothed and segments were condensed into genomic regions with copy number state per array per case using the ‘CGHregions’ R package[3] in Bioconductor[2] (averror = 0). Genomic regions that contained less than 100 probes were removed. Genomic regions were retained if the number of base pairs per array probe was within the median probe density plus/minus double interquartile range of all genomic regions of all cases analysed (103.5 to 104.4 base pairs per probe in this series, i.e. median(all regions) ±2*IQR[all regions]) to filter genomic regions with extreme tiling densities.

# Tumour-specific mirrored B-allele frequency profile
In order to generate a tumour specific mirrored B-allele frequency (mBAF) profile the estimation of aberrant cell fraction (A) generated by ASCAT[4] (min. 80, max. 95%) was used to define a threshold (T) above which probe mBAF values are considered non-informative homozygous SNPs and removed,

1/1+(1-A) = T.

The remaining mBAF values (B) that are greater than 0.56 are scaled between 0.5 and the non-informative homozygous SNP threshold (T) to give a tumour specific mBAF value (M),

((0.5*B) - 0.5) / (T - 0.5) + 0.5 = M.

Mean mBAF was calculated for each copy number loss (copy number = 1) and gain (copy number = 3) state genomic region and for each region, calls were reversed to normal if the mean mBAF was not greater than 0.58 and 0.64 for gain and loss states respectively.

# Copy number neutral loss of heterozygosity detection
As ‘CGHcall’ uses LRR data to generate copy number segments, copy number neutral loss of heterozygosity events, that do not alter LRR, are not detected. Therefore we developed an algorithm to identify LOH in genomic regions of a total copy number of 2. 
Here we counted the number of probes with a processed mBAF <= 0.66 in windows of 100 SNP probes, sliding by 10 probes per chromosome. Windows with 1% or less probes with mBAF <= 0.66 were considered to not be in the ‘AB state’ of normal allelic balance. Consecutive windows of allelic imbalance were merged for each array and boundaries across different arrays from the same tumour were smoothed to make these regions comparable.
These newly identified regions were then incorporated into the genomic regions generated from the LRR data, if (1) an allelic imbalance region was called within a genome region with a copy number state of 2 (as described earlier) and (2) the allelic imbalance region contained equal to or more than 1000 SNP probes (more than 90 consecutive windows called as allelic imbalance).

# Calculation of major and minor copy numbers
Once the genomic regions of copy number states for a case were calculated, we calculated major (A) and minor (B) allele copy numbers for each region using the total copy number (T - as called by ‘CGHcall’) and the mean mBAF (M) for the probes within each genomic region. We used the following equation to determine major and minor allele copy numbers;

round(TM) = A,
T-A = B.

The major and minor allele copy numbers per chromosome per array were used as input for MEDICC[5].

1	van de Wiel MA, Kim KI, Vosse SJ, van Wieringen WN, Wilting SM, Ylstra B. CGHcall: calling aberrations for array CGH tumor profiles. Bioinformatics 2007; 23: 892–4.

2	Gentleman RC, Carey VJ, Bates DM, et al. Bioconductor: open software development for computational biology and bioinformatics. Genome Biol 2004; 5: R80.

3	van de Wiel MA, van Wieringen WN. CGHregions: dimension reduction for array CGH data with minimal information loss. Cancer Inform 2007. http://www.ncbi.nlm.nih.gov/pmc/articles/PMC2675846/.

4	Van Loo P, Nordgard SH, Lingjærde OC, et al. Allele-specific copy number analysis of tumors. Proc Natl Acad Sci U S A 2010; 107: 16910–5.

5	Schwarz RF, Trinh A, Sipos B, Brenton JD, Goldman N, Markowetz F. Phylogenetic quantification of intra-tumour heterogeneity. PLoS Comput Biol 2014; 10: e1003535.
