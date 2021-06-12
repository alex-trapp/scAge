![image](https://user-images.githubusercontent.com/83630286/121631269-bfda1d80-ca4c-11eb-830b-6c8378d18160.png)

# scAge

scAge is a probabilistic framework for epigenetic age profiling at single-cell resolution, developed in Python. <br>
This tool leverages the linear relationship of methylation and chronological age at some CpGs in the DNA, and uses these models to predict
epigenetic age in single cell data. <br>
To learn more about the underlying algorithms driving scAge, consult our [preprint on bioRxiv](https://www.biorxiv.org/content/10.1101/2021.03.13.435247v1).

## Installation
To install scAge and associated data, clone the GitHub repository. All required functions are stored within scAge.py.

`git clone https://github.com/alex-trapp/scAge`

To run scAge functions, import the module into a Python script or Jupyter notebook

`import scAge`

In order to use scAge, the following packages need to be installed:

`numpy` (developed with v1.20.2)
`pandas` (developed with v1.2.4)
`scipy` (developed with v1.6.3)
`sklearn` (developed with v0.24.2) 
`tqdm` (developed with v4.60.0)

## Usage
scAge is a workflow that enables epigenetic age prediction in single cells using a combination of linear models to estimate age.

3 example Jupyter notebooks are provided in the `notebooks` directory:
`example_process_coverage` --> processing .cov files from Bismark into processed binary methylation matrices
`example_construct_reference` --> constructing a reference set of linear models from a bulk methylation matrix
`example_run_scAge` --> predicting epigenetic age in single cells

These notebooks use a sample of cells from the [Gravina et al. study](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-016-1011-3), shown in Figure 2 of the manuscript. Required data to run these example scripts is provided:
* Raw .cov files are provided in `sc_data_raw`
* Processed binary methylation matrices are provided in `sc_data_processed`
* Raw bulk data used to  is provided in `bulk`
* Processed reference matrix is provided in `train`

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
* `training_DNAm_matrix` --> the input bulk reference matrix <br>
* `output_path` --> desired full path <br>
* `n_cores` -->  the number of cores that scAge should use via parallel processing <br>
* `chunksize` --> the number of individual CpG-age series that should be distributed at once to each worker <br>

### Loading single-cell methylomes
scAge requires binary methylation matrices as input for the epigenetic age profiling algorithm. These binary matrices can be obtained
by processing existing .cov files produced by [Bismark](https://www.bioinformatics.babraham.ac.uk/projects/bismark/). Depending on the tool used,
final single-cell methylome files may have slightly different formats. The Bismark .cov file format is the following:

Chromosome | Position 1 | Position 2 | Methylation level | Methylated counts | Unmethylated counts
:---: | :---: | :---: | :---: | :---: | :---: 
11 | 3100225 | 3100225 | 100 | 1 | 0
11 | 3101286 | 3101286 | 0 | 0 | 2

This may require slightly modifying the function provided. 
In the end, scAge requires a .tsv file or pandas datraframe with two columns:

ChrPos | MetLev
:---: | :---:
chr1_3037802 | 1
chr19_61305429 | 0
... | ...

To load and process .cov files generated by Bismark, run
```
process_cov(cov_directory, 
            n_cores = 10,
            max_met = 100,
            split = ".",
            chunksize = 1,
            binarization = "hard",
            output_path = None):
```

* `cov_directory` --> the path to the directory where .cov single-cell methylation files are stored <br>
* `n_cores` --> the number of cores to use for parallel processing <br>
* `maxmet` --> the maximum methylation value (normally, methylation ratios from Bismark range from 0 to 100) <br>
* `split` -->  desired string to split the file name on for single-cell name generation (i.e. "SRR3136624.cov" --> "SRR3136624") <br>
* `chunksize` --> number of coverage files that will be fed into a single worker process at a time <br>
* `binarization` --> choice of hard vs. soft.
                   Both methods involve dropping methylation values of 0.5.
                   "hard binarization" rounds remaining non-binary values to 0 or 1 depending on proximity, while
                   "soft binarization" discards remaining non-binary values <br>
* `output_path` --> the output directory to which processed .tsv binary matrices should be written to
If `output_path` is set to `None`, named binary methylation matrices are returned <br>

### Predicting epigenetic age
The core of scAge is run_scAge, a function that enables epigenetic age predictions from a previously computed set of binarized single-cell
methylome profiles and a reference linear regression matrix generated from bulk data. To run this function:
```
scAge.run_scAge(single_cell_dir_or_dict,
                single_cell_set_name,
                reference_data,
                selection_mode = "percentile",
                CpG_parameter = 1,
                zero_met_replacement = 0.001,
                one_met_replacement = 0.999,
                min_age = -20,
                max_age = 60,
                age_step = 0.1,
                n_cores = 3,
                uncertainty = 1,
                output_path = None,
                chunksize = 1)
```

* `single_cell_dir_or_dict` --> full path to a directory or a dictionary of processed single-cell methylomes, created with `process_cov`
* `single_cell_set_name` --> desired name of the scAge run (i.e. 'dataset_x)
* `reference_data` --> full path to the reference matrix created with `create_reference`
* `selection_mode` --> one of `percentile`, `numCpGs`, `cutoff`, which determines which CpG selection mode should be used in the algorithm
* `CpG_parameter` --> parameter accompanying selection mode (1 in percentile mode --> top 1% age-correlated CpGs)
* `zero_met_replacement1` --> if the linear model goes below 0, this value replaces the probability
* `one_met_replacement` --> if the linear model goes above 1, this value replaces the probability
* `min_age` --> minimum age for probability computations
* `max_age` --> maximum age for probability computations
* `age_step` --> step (in months) that probability computations should be performed at (age_step = 0.1 --> (min_age, min_age + age_step, ..., max_age)
* `n_cores` --> number of cores that should be used for parallelization
* `uncertainty` --> the uncertainty metric that should be used to compute upper and lower bounds (higher value --> wider interval)
* `output_path` --> full path to the directory where predictions and the report file should be written to
* `chunksize` --> number of single-cell methylomes that should be passed to each parallel worker at once

The output of run_scAge is a .tsv matrix, which contains a number of columns which detail the internal computations of the algorithm.
The leftmost columns are the most crucial, including: 
* `PredictedAge` --> scDNAm epigenetic age predictions
* `MeanMet` --> mean global methylation of the single cell
* `CellCoverage` --> number of CpGs covered in the single cell

Cell | PredictedAge | MeanMet | CellCoverage
:---: | :---: | :---: | :---: |
SRR3136627 | 0.5 | 0.663169 | 3914949 |
SRR3136659 | 4.0 | 0.683454 | 799350 |
SRR3136628 | 25.0 | 0.695256 | 2511084 |
