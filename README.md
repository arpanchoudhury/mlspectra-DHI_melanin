# mlspectra-DHI_melanin

Designing of virtual chemical space of DHI melanin and machine learning (ML) prediction of their entire UV-visible spectra and thermodynamic stability. 
The step-wise procedure is described below.

## 1. Structure generation using combinatorial algorithm
```
python generator_dimer.py
python generator_trimer.py
python generator_tetramer.py
```
The codes need to be run in the above order. This will generate all possible unique dimers, trimers and tetramers structures from the monomers. 
In our study, we considered only the tetramers.
## 2. ML training data 
A subset of generated tetramers are first relaxed with the UFF force field and then optimized using the B3LYP/6-31G(d) level. The final optimized geometries are contained in ``` ML/optimized_geometry/```. Their TD-CAM-B3LYP/6-31G(d) electronic spectra of the lowest 70 singlet excitation energies and the corresponding oscillator strengths are contained in ``` ML/spectra/```.
```
python kmeans_clusters.py
```
for each of the DHICA, DKICA and MKICA.
Subsequently, subclustering of OH dihedral angles within each cluster can be performed based on the clustering result by running 
```
python kmeans_subclusters.py
```
### Spectrum binning and averaging
The next steps are (i) to bin the spectra within a given range and (ii) to calculate the mean spectra and mean structures of every subclusters. Step (ii) is required to discard the redundancies in the dataset. It ensures that we take only one entry from each subcluster instead of all similar/like entires. 
This can be done by running 
```
python subclusters_averaging.py
```
It will create two .csv files; one for geometries and one for spectra.
### Screening important clusters
In a given spectral range, instead of taking all the clusters, we can take most important few clusters which have higher intensity values than others. 
```
python screen_imp_clusters.py
```
will print the desired number of important clusters.
### Final dataset generation
Final dataset generation for ML regarding the important clusters, can be done by running 
```
python make_ML_input-output.py
```
This will create two binary files which contain ML input spectra and output geometries.
## Machine learning training and prediction
Final ML training and predictions can be done by running the following command, 
```
python run_KRR-ML.py --minrange 290 --maxrange 300 --Ntrain 10000
```
`--minrange` and `--maxrange` specify the spectral range and `--Ntrain` specifies the training set size.
## Requirements
The codes were tested with Python 3.8.10 on HP-OMEN [Intel(R) Core(TM) i5-10300H CPU @ 2.50GHz] as well as with Python 2.7.5 on linux workstation [Intel(R) Xeon(R) CPU E5-1650 v4 @ 3.60GHz]. Below are the versions of python modules used on the HP-OMEN
```
numpy 1.22.3
pandas 1.4.2
scipy 1.8.0
scikit-learn 1.0.2
joblib 1.1.0
```
