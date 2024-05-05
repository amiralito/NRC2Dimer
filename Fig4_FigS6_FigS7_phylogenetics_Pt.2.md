1. To extract the NRC clade and NRC helpers, we first made a [phylogenetic tree of all NLRs from 124 genomes based on the NBARC domains](phylogenetics/NBARC_ref_filtered.newick). Then a well-supported branch comprising both NRC helpers and sensors was selected and extracted *Dendroscope* (Options > Advanced Options > Extract Subnetwork...) from the phylogenetic tree:

![NBARC_all](extras/NBARC_all.png)

2. The extracted sequences were aligned again and using the alignment a new phylogenetic tree of the NRC network was made. Finally, the well-supported branch containing all reference NRC helpers was selected and extracted in the same way as before for the NRC helper clade:

![NRC_all](extras/NRC_all.png)