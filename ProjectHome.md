The input is a multiple sequence alignment (MSA) of protein sequences and the output is a division of the MSA into subsets of MSA rows corresponding to predicted subfamilies. Each row in the MSA is used to derive an amino acid profile, followed by agglomerative clustering of the profiles into a hierarchical tree using relative entropy between profile distributions. A minimum encoding cost measure is used to derive the point during the agglomeration to determine a division of the tree into subtrees, defining subfamilies. The algorithm is described in [Brown DP, Krishnamurthy N, Sjölander K, "Automated Protein Subfamily Identification and Classification," PLoS Computational Biology 2007, 3(8): e160 doi:10.1371/journal.pcbi.0030160](http://www.ploscompbiol.org/article/info:doi/10.1371/journal.pcbi.0030160)


---


This software has been produced jointly by the [Berkeley Phylogenomics Group (BPG)](http://phylofacts.berkeley.edu/) headed by Dr Kimmen Sjölander and DuPont Pharmaceuticals Bioinformatics Group, headed by Dr Jean-Francois Tomb.