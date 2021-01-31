# About

## Plant LRR-RLKs structure

In plants, the leucine-rich repeat receptor-like protein kinases (LRR-RLKs) are composed of an extracellular LRR domain (ECD) responsible for ligand binding, a single membrane-spanning helix (TM), and a cytoplasmic kinase domain (KD). The extracellular LRR motifs share canonical plant specific LRR consensus sequences (CS): LxxLxLxxNxL(s/t)GxLPxxLGxLxx (x represents any amino acid), where the LxxLxLxxNxL(s/t)GxLP segments are the plant-specific LRR highly conserved segment (HCS).

## The Phyto-LRR prediction program

Since the program is based on the position specific scoring matrix (PSSM) algorithm and its training dataset are 4000 LRR motifs HCSs from LRR-RLKs of 17 land plant species (Tab. 1), this [plant LRR prediction program](/findlrr) could identify plant LRRs effectively, especially for LRR motifs in plant LRR-RLKs. The program also has an offline version in the form of a python package, which has been published on PyPI: https://pypi.org/project/predict-phytolrr/.

## The database of the plant LRR repeats

In this [database](/), a total of 55457 LRR motifs were identified from 3987 LRR-RLK genes in 17 represented land plants (Tab. 1) using the "Phyto-LRR-predict" program, which is freely accessible at the [Predict LRRs page](/findlrr). The sequence second structure and the soluble accessibility were obtained from the SSpro and ACCpro20 programs in the SCRATCH suite ([Magnan CN & Baldi P., 2014](https://academic.oup.com/bioinformatics/article/30/18/2592/2475618)). Moreover, since the ectodomains are heavily asparagine-linked (N-) glycosylated, the database also integrated with the predicted canonical N-glycosylation sites (Asn-x-Ser/Thr (x ≠ Pro)).

**Table 1: Number of LRR-RLK genes in 17 represented land plants** 

|Species|Five Digit Code|Number of LRR-RLKs|
|-------|---------------|------------------|
|Amborella trichopoda|AMBTR|108|
|Arabidopsis lyrata|ARALY|207|
|Arabidopsis thaliana|ARATH|213|
|Brachypodium distachyon|BRADI|222|
|Brassica rapa|BRARA|273|
|Glycine max|GLYMA|424|
|Marchantia polymorpha|MARPL|102|
|Medicago truncatula|MEDTR|316|
|Oryza sativa ssp. Indica|ORYSI|287|
|Oryza sativa ssp. Japonica|ORYSJ|281|
|Phoenix dactylifera|PHODC|149|
|Physcomitrella patens|PHYPA|162|
|Populus trichocarpa|POPTR|373|
|Selaginella moellendorffii|SELML|138|
|Solanum lycopersicum|SOLLC|213|
|Solanum tuberosum|SOLTU|291|
|Zea mays|MAIZE|228|
|Total| |3987|


