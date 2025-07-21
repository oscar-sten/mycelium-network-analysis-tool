# mycelium-network-analysis-tool
A repository with the code and data for reproducing the results in: A Ridge-Based Detection Algorithm with Filament Overlap Identification for 2D Mycelium Network Analysis.

Article authors: Oscar Sten, Emanuela del Dottore, Nicola Pugno, 
and Barbara Mazzolai.

DOI: https://doi.org/10.1016/j.ecoinf.2024.102670

Code author: Oscar Sten, unless otherwise indicated.

# Folders:
## Src
Contains:

**fungal_network_2D_qualitative.m** for reproducing the qualitative
 results.

**branch_crossing_identification_test.m** for reproducing results in
 Table 3.

**branch_crossing_identification_test_2.m** for reproducing results in
 Table 4.

**branch_crossing_matching_test.m** for reproducing results in Table 5.

**comparison_length_measurement.m** for reproducing results in Figure 11.

**comparison_morph_graphs.m** for reproducing the results presented in 
Section 3.4.  

**create_GT_skeletons_for_benchmark_on_Cardini_et_al_images.m** converts 
the manual annotations in `data\lbl_length_Cardini_et_al_images\` to the 
skeletons found in `ground_truth_skeletons_Cardini_et_al_images`.

**test_detection_benchmark_on_Cardini_et_al_images** performs benchmark 
test on the image set in `/data/Cardini_et_al_images`.

**test_detection_benchmark_on_GUFI_1_images** performs benchmark test on
the image set in `/data/Ghent_University_Fungal_Images_I`.

**test_detection_benchmark_on_Vidal_et_al_images** performs benchmark 
test on the image set in `/data/Vidal_et_al_images`.

`fcn/` contains all necessary functions.

`fcn/external/` functions originally authored by external authors.



## result_tables
Contains full results of the bench marking test results displayed in
 Table 2.

## data
Contains all data including images and annotations.

### Cardini_et_al_images
Images of *R. Irregularis* samples downloaded from 
[Repository](https://gitlab.iit.it/EDelDottore/hylength/-/tree/master/HyLength_validation_images/Rhizophagus%20irregularis_in%20vivo%20hyphae).\
Images #1 and #2 displayed in Figure 9 are here called
 `s.jpg` and `o.jpg` respectively.

Original publication: Cardini, A., Pellegrino, E., Del Dottore, E. et al. 
HyLength: a semi-automated digital image analysis tool for measuring the 
length of roots and fungal hyphae of dense mycelia. Mycorrhiza 30, 
229–242 (2020). https://doi.org/10.1007/s00572-020-00956-w 


### Ghent_University_Fungal_Images_I
Images in the GUFI-1 image set downloaded from:
 http://dx.doi.org/10.13140/RG.2.1.4441.1607 

Full documentation is available at source.

### Vidal_et_al_images
Contains images of *A. Oligospora* aquired by Vidal et al. (2019).

Source: https://github.com/hsueh-lab/FFT/tree/master/GT_Mycelium 

The images not available in this repository were sent to us by the authors.

Original publication: Vidal-Diez de Ulzurrun G, Huang TY, Chang CW, Lin HC,
 Hsueh YP. Fungal feature tracker (FFT): A tool for quantitatively
 characterizing the morphology and growth of filamentous fungi.
 PLoS Comput Biol. 2019 Oct 31;15(10):e1007428.
 doi: 10.1371/journal.pcbi.1007428. 

### lbl_length_Cardini_et_al_images
Contains line annotations of the images from `Cardini_et_al_images`, 
the annotations were perfromed with the MATLAB image labler tool 
by Oscar Sten in 2023.

### lbl_intersections_test
Contains categorical box annotations of the images `s.jpg` and `o.jpg`, 
the annotations were perfromed with the MATLAB image labler tool
 by Oscar Sten in 2023.

### lbl_matching_test
Contains pixel annotations of the image `s.jpg`, the annotations were
 perfromed with the MATLAB image labler tool by Oscar Sten in 2023.

### ground_truth_skeletons_Cardini_et_al_images
Contains ground truth skeletons for the images in `Cardini_et_al_images`, 
the annotations were perfromed with the MATLAB image labler
 tool by Oscar Sten in 2023.


## Remarks
The branch/crossing detection and branch/crossing matching tests are
not fully automated, since a couple of cases were too complicated to 
annotate. In these cases the person performing the test must determine 
visually is the algorithm's assessment is correct or not.

1) In image `s.jpg` there is a crossing in the bottom left corner whose
 matching must be assessed manually. Running 
`branch_crossing_matching_test.m` with the parameters reported in the
 manuscript makes this crossing wrongly matched, so one "fault"
 must be added.

2) In image`o.jpg`, a filament in the top right corner traverses
 multiple structures and must be assessed maunally.  Running 
`branch_crossing_identification_test_2.m` with the parameters 
reported in the manuscript makes this crossing not detected, so one
 "false negative" must be added.

