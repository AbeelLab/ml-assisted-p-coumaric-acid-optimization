# p-Coumaric acid project
Exploration and exploitation in large combinatorial design spaces

Project structure


## Data
1. Promoter strength screening: [PromoterScreen_FCS_Elif](data/raw/PromoterScreen_FCS_Elif)
2. Data for choosing the genes for optimization: [GeneTargets](data/raw/GeneTargets)
- The Brenda database for brendapyrser
- String-DB Enrichment analysis for the genes coming out of FBA (one seperate file for PAL/C4H 
since these are natively from a different organism.)
1. Library clones (in 4tu repo)
2. GFP expression: [quantified_promoter_strength_rewritten.csv](results/GFP_PromoterStrengths/quantified_promoter_strength_rewritten.csv)

## Scripts
#### Sequence pool analysis
1. Analysis of the pooled sequencing of the four library designs: [SequencePoolAnalysis(Irsan).py](scripts/SequencePoolAnalysis(Irsan))

#### Preprocessing for ML
Sequencing data is stored in the [4tu-repo](https://data.4tu.nl/private_datasets/10.4121/fa0782ab-760d-4fa7-babf-09bdaab0f509). We also provide the count matrix and the numeric matrix used 
for machine learning. Gff files are made available upon request.

1. Processing fcs files to (mean) promoter strengths: [PromoterScreenAnalysis.py](scripts/PromoterScreenAnalysis.py)
2. Finding gene targets from GEM + additional info on Thermodynamics and Literature enrichment: 
[FindGeneTargets_YeastGEM.py](scripts/FindGeneTargets_YeastGEM.py)
3. Processing gff files to [pro-orf-ter](scripts/ConstructCountMatrix.py) count matrix
4. Merging previous WUR round with TUD round (strain count matrix)
[MergeTUDWURdata.py](scripts/MergeTUDWURdata.py). Outputs the file batch_corrected_screening_results_integrated.csv, which is used for training.

5. Convert to a numeric matrix: [ConstructNumericMatrix.py](scripts/ConstructNumericMatrix.py)
6. Batch-correction script for WUR data [BatchCorrectWUR.py](scripts/BatchCorrectWUR.py)
7. Batch correction and some plots of TUD round: [BatchCorrection.py](scripts/BatchCorrection.py)

8. Sampling the designs for sequencing, after screening [SampleDesignsForSequencing.py](scripts/SampleDesignsForSequencing.py).

#### Rescreening and rebuild
1. Rescreened top 86: [RescreeningTopProducers](scripts/RescreeningTopProducers.py)
2. Rebuild

#### Library visualizations
10. Gene count fraction: [AssessGeneContent.py](scripts/AssessGeneContent.py)

#### ML + Feature importance
11. Dense Weight implementation + XGBoost training [DenseWeightTraining.py](scripts/DenseWeightTraining.py). 
Remeasured top 100 strains were also integrated here. These were lower than the first (n=1) measurement round.
12. TODO Feature importance scripts need to be rewritten into one py file
From old project: 1704_combinatorial_interrogation.py
010724_individualcontributions_genes.py, 2504_combination_analysis.py

#### Round 2 (DoE)
13. analysis of new designs: [AnalysisValidationRound.py](scripts/AnalysisValidationRound.py). 


## Results

