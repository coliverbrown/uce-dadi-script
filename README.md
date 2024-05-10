# uce-dadi-script
dadi model scripts used for Beringia UCE population genetics analysis  

We used diffusion analysis for demographic inference (δaδi v2.3.0; Gutenkunst et al., 2009) to identify best-fit demographic models. For demographic analyses, we used the thinned VCF file with Z-linked loci removed. We tested eight models of divergence: A) neutral (populations never diverged); B) split with no migration (divergence without gene flow); C) split with migration (divergence with gene flow that is bidirectionally symmetric); D) split with bidirectional migration (divergence with gene flow that is bidirectionally asymmetric; E) secondary contact with migration; F) secondary contact with bidirectional migration; G) split with exponential population growth and no migration; H) split with exponential population growth and migration (Figure 2). 

We then ran the best-fit model with 100 bootstrap replicate data sets using a bootstrapping script. 
