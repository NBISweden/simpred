# simpred

SNV impact predictor

## Synopsis
Simpred takes a snpEff-generated annotation file and for every SNV except the ones annotated as "modifier" it computes the potential substitution effect based on the following sets of descriptors:

1. The aaI7 set of descriptors-based effect. [[1]](#1) 
2. The Exchangeability of Amino Acids in Proteins. [[2]](#2) 
3. The Sneath dissimilarity index. [[3]](#3)

## Usage

```simpred mySnpEff_file.vcf > result.csv```

## Output

A csv file containing the result:

Scaffold | Coord | Ref | Var | Type Effect | Transcript | Ref_aa | Coord_aa | Var_aa | Ref_aa_abbrev | Var_aa_abbrev | aaI7 | exchgb_ref_var | exchgb_var_ref | sneath_dissim
-------- | ----- | --- | --- | ----------- | ---------- | ------ | -------- | ------ | ------------- | ------------- | ---- | -------------- | -------------- | -------------
Sc9M7eS_1_HRSCAF_2 | 362173 | A | C | missense_variant | MODERATE | mRNA20769 | Ile | 46 | Leu | I | L | 0.05 | 0.52 | 0.34 | 0.11
Sc9M7eS_1_HRSCAF_2 | 412292 | T | A | missense_variant | MODERATE | mRNA20769 | Asp | 163 | Glu | D | E | 0.16 | 0.16 | 0.46 | 0.16
Sc9M7eS_1_HRSCAF_2 | 456817 | T | C | missense_variant | MODERATE | mRNA20770 | Lys | 3194 | Arg | K | R | 0.25 | 0.52 | 0.3 | 0.31
Sc9M7eS_1_HRSCAF_2 | 459865 | G | A | missense_variant | MODERATE | mRNA20770 | Ser | 2178 | Leu | S | L | 0.54 | 0.67 | 0.53 | 0.51
Sc9M7eS_1_HRSCAF_2 | 461297 | A | C | missense_variant | MODERATE | mRNA20770 | Ser | 1701 | Ala | S | A | 0.15 | 0.39 | 0.23 | 0.36
Sc9M7eS_1_HRSCAF_2 | 463630 | T | C | missense_variant | MODERATE | mRNA20770 | Asn | 923 | Ser | N | S | 0.28 | 0.38 | 0.44 | 0.33
Sc9M7eS_1_HRSCAF_2 | 464138 | G | A | missense_variant | MODERATE | mRNA20770 | Arg | 754 | Cys | R | C | 1.0 | 0.52 | 0.64 | 0.8
Sc9M7eS_1_HRSCAF_2 | 478428 | C | T | missense_variant | MODERATE | mRNA20770 | Val | 389 | Ile | V | I | 0.15 | 0.22 | 0.15 | 0.16
Sc9M7eS_1_HRSCAF_2 | 797450 | C | T | missense_variant | MODERATE | mRNA20771 | Pro | 323 | Ser | P | S | 0.24 | 0.61 | 0.38 | 0.53  

## References
<a id="1">[1]</a> [Rudnicki WR, Komorowski J (2010). Feature Synthesis and Extraction for the Construction of Generalized Properties of Amino Acids. Rough Sets and Current Trends in Computing Vol.3066, ed Tsumoto S., Słowiński R., Komorowski J. G-BJ. (Springer, Berlin, Heidelberg), pp 786-791.](https://www.researchgate.net/publication/220801316_Feature_Synthesis_and_Extraction_for_the_Construction_of_Generalized_Properties_of_Amino_Acids)  
<a id="2">[2]</a> [Lev Y. Yampolsky and Arlin Stoltzfus. GENETICS August 1, 2005 vol. 170 no. 4 1459-1472;](https://doi.org/10.1534/genetics.104.039107)  
<a id="3">[3]</a> [Sneath, P. H. (1966-11-01). Relations between chemical structure and biological activity in peptides. Journal of Theoretical Biology. 12 (2): 157-195.](https://doi.org/10.1016/0022-5193(66)90112-3)  