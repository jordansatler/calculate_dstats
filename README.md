# calculate_dstats

This script will calculate D-statistics (ABBA/BABA) for sets of four taxa. Data file is assumed to be in phylip format and SNPs are treated as if they are unlinked. An input phylogeny is required as is the designation of an outgroup taxon. Script requires the ETE Toolkit (Huerta-Cepas et al. 2016) to be installed . Script will generate all possible combinations of 3 ingroup taxa, then using the phylogeny as rooted by the outgroup taxon, will place taxa in phylogenetic order for D-stat calculations. To assess significance, script uses 1,000 bootstrap replicates to test for significant deviation from D = 0.


User needs to provide a snps data file, a tree file, and the name of the outgroup taxon.

usage:  
```python
    python dstats.py --data snps.file --tree tree.file --outgroup outgroup
```

***
References:

Huerta-Cepas, J., Serra, F., & Bork, P. 2016. ETE 3: reconstruction, analysis, and visualization of phylogenomic data. Molecular Biology and Evolution, 33:1635-1638.
