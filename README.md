# simpred

SNV impact predictor

## Synopsis
Simpred takes a snpEff-generated annotation file and for every SNV annotated as "MODERATE" it computes the potential substitution effect based on the following sets of descriptors:

1. The aaI7 set of descriptors-based effect. [[1]](#1) 
2. The Exchangeability of Amino Acids in Proteins. [[2]](#2) 
3. Sneath index. [[3]](#3)



## References
<a id="1">[1]</a> [Rudnicki WR, Komorowski J (2010). Feature Synthesis and Extraction for the Construction of Generalized Properties of Amino Acids. Rough Sets and Current Trends in Computing Vol.3066, ed Tsumoto S., Słowiński R., Komorowski J. G-BJ. (Springer, Berlin, Heidelberg), pp 786-791.](https://www.researchgate.net/publication/220801316_Feature_Synthesis_and_Extraction_for_the_Construction_of_Generalized_Properties_of_Amino_Acids)
<a id="2">[2]</a> [Lev Y. Yampolsky and Arlin Stoltzfus. GENETICS August 1, 2005 vol. 170 no. 4 1459-1472;](https://doi.org/10.1534/genetics.104.039107)
<a id="3">[3]</a> [Sneath, P. H. (1966-11-01). Relations between chemical structure and biological activity in peptides. Journal of Theoretical Biology. 12 (2): 157-195.](https://doi.org/10.1016/0022-5193(66)90112-3)