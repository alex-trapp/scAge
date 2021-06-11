![image](https://user-images.githubusercontent.com/83630286/121631269-bfda1d80-ca4c-11eb-830b-6c8378d18160.png)

# scAge

scAge is a probabilistic framework for epigenetic age profiling at single-cell resolution.
To learn more about the underlying algorithms driving scAge, consult our [preprint on bioRxiv](https://www.biorxiv.org/content/10.1101/2021.03.13.435247v1).

## Installation
To install scAge, clone the GitHub repository. All required functions are stored within scAge.py.

`git clone https://github.com/alex-trapp/scAge`

## Usage
scAge is a workflow that enables epigenetic age prediction in single cells using a combination of linear models to estimate age.

### Training
CpG-specific linear models are first calculated using a reference bulk dataset, which may contain some missing values.
Using single tissue or single cell-type datasets is preferred to improve prediction accuracy, although multi-tissue datasets
can also be used for training.

We provide a few pre-computed training datasets inside of the **training** directory.

In order to train a custom set of linear models from a given DNAm matrix, run
```
construct_reference(training_DNAm_matrix,
                    output_path,
                    n_cores = 30,
                    chunksize = 100)
```
This function takes as input a pandas dataframe DNAm matrix, with rows as samples and columns as CpGs (in the form chr9_85324737). 
Methylation values must be in the range from 0 (fully unmethylated) to 1 (fully methylated)
This dataframe must contain an "Age" column, which is used to compute correlations and linear regressions. <br>
`training_DNAm_matrix` --> the input bulk reference matrix <br>
`output_path` --> desired full path <br>
`n_cores` -->  the number of cores that scAge should use via parallel processing <br>
`chunksize` --> the number of individual CpG-age series that should be distributed at once to each worker <br>

### Loading single-cell methylomes
scAge requires binary methylation matrices as input for the epigenetic age profiling algorithm. These binary matrices can be obtained
by processing existing .cov files produced by [Bismark](https://www.bioinformatics.babraham.ac.uk/projects/bismark/). Depending on the tool used,
final single-cell methylome files may have slightly different formats. This may require slightly modifying the function provided. 
In the end, scAge requires a .tsv file or pandas datraframe with two columns:
provide additional info on Bismark files

ChrPos | MetLev
------------ | -------------
chr1_3037802 | 1
chr19_61305429 | 0
... | ...

To load in .cov files processed by Bismark, run
```
load_cov(cov_directory, 
         n_cores = 10,
         max_met = 100,
         split = ".",
         chunksize = 1,
         binarization = "hard",
         output_path = None):
```

`cov_directory` --> the path to the directory where .cov single-cell methylation files are stored
`n_cores` --> the number of cores that should be used to simulateneously load and process .cov files
`maxmet` --> the maximum methylation value (
