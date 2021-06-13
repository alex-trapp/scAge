![image](images/scAge_schematic_GitHub.jpg)
# scAge

scAge is a probabilistic framework for profiling epigenetic age at single-cell resolution, developed in Python. <br> <br>
This tool leverages the relationship between DNA methylation in bulk samples and chronological age to perform
probability-based profiling of epigenetic age in intrinsically sparse and binarized single-cell data. <br> <br>
This approach constitutes the first available single-cell clock, and is both highly scalable and flexible. <br> <br> 
It can be used on any number of cells and can be trained on any methylation-age dataset. <br> <br>
To learn more about the underlying algorithms driving scAge, please consult our [preprint on bioRxiv](https://www.biorxiv.org/content/10.1101/2021.03.13.435247v1).

## Installation
To install scAge and associated data, please clone the GitHub repository:

`git clone https://github.com/alex-trapp/scAge.git`

This will download all required data to use and test the software. <br>

For ease of use, all functions needed to run the full scAge pipeline are included within scAge.py

## Usage

To run scAge, first import the module into a Python script or Jupyter notebook:

`import scAge`

In order to use the functions provided in scAge, the following packages need to be installed:

`numpy` (developed with version 1.20.2) <br>
`pandas` (developed with version 1.2.4) <br>
`scipy` (developed with version 1.6.3) <br>
`sklearn` (developed with version 0.24.2) <br>
`tqdm` (developed with version 4.60.0) <br>

We also recommend `seaborn` and `matplotlib` to visualized epigenetic age predictions in Python. If desired, predicted epigenetic age output files
can also be analyzed in other programming languages (i.e. R).

scAge is a workflow that enables epigenetic age prediction in single cells using a combination of linear models to estimate age.

3 example Jupyter notebooks are provided in the `examples/notebooks` directory: <br>
* `process_coverage_notebook.ipynb` --> processing .cov files from Bismark into processed binary methylation matrices <br>
* `construct_reference_notebook.ipynb` --> constructing a reference set of linear models from a bulk methylation matrix <br>
* `example_run_scAge_notebook.ipynb` --> predicting epigenetic age in single cells <br>

These notebooks use a sample of cells from the [Gravina et al. study](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-016-1011-3),
described in Figure 2 of our manuscript. <br>
Required data to run these example scripts is provided:
* Raw .cov files of single-cell methylomes are located in `sc_data_raw`
* Processed binary methylation matrices are located in `sc_data_processed`
* Raw bulk data for liver used to construct reference models is located in `bulk`
* Processed reference matrix for liver is located in `train`

The functions driving scAge are explained and documented below. Default values for parameters are shows (i.e. `cores = 30`)

## Training
CpG-specific linear models are first calculated using a reference bulk methylation matrix, which may contain some missing values.
From our findings, using single tissue or single cell-type datasets is preferred to improve prediction accuracy,
although multi-tissue datasets may also be used for training.

I provide a pre-computed training reference dataset for C57BL/6J mice livers inside of the `train` directory, and I will add more shortly.

In order to train a custom set of linear models from a given DNAm matrix, run:

```
construct_reference(training_DNAm_matrix,
                    output_path,
                    n_cores = 30,
                    chunksize = 100)
```

where
* `training_DNAm_matrix` --> input bulk methylation matrix <br>
* `output_path` --> desired full path to the output reference matrix (/path/to/file.tsv) <br>
* `n_cores` --> number of cores to use via parallel processing <br>
* `chunksize` --> number of individual CpG methylation series to distribute at once to each worker <br>

This function takes as input a pandas dataframe DNAm matrix, with rows as samples and columns as CpGs (in the form chr9_85324737). <br>
Methylation values must be in the range from 0 (fully unmethylated) to 1 (fully methylated). <br>
This dataframe must contain an "Age" column, which is used to compute correlations and linear regressions. <br>
An example bulk matrix of C57BL/6J mice livers is provided in the `bulk` directory <br>


## Loading single-cell methylomes
scAge requires binary methylation matrices as input for the core epigenetic age profiling algorithm. These binary matrices can be obtained
by processing existing .cov files produced by [Bismark](https://www.bioinformatics.babraham.ac.uk/projects/bismark/) with our function
`process_coverage`. However, there are other tools to process FASTQ bisulfite sequencing data, such as [BSSeeker](https://github.com/BSSeeker/BSseeker2) and [methylpy](https://github.com/yupenghe/methylpy). Final single-cell methylome files may have slightly different formats. The Bismark .cov file format is the following
(the columns are not named in the files):

Chromosome | Position 1 | Position 2 | Methylation level | Methylated counts | Unmethylated counts
:---: | :---: | :---: | :---: | :---: | :---: 
11 | 3100225 | 3100225 | 100 | 1 | 0
11 | 3101286 | 3101286 | 0 | 0 | 2

Native support for a variety of input formats will be added shortly. For now, I recommend processing single-cell methylation data with [Bismark](https://www.bioinformatics.babraham.ac.uk/projects/bismark/) or converting existing methylation data to a Bismark format (shown above). <br> <br>

Most of the publicly available data that we used for our study was in the .cov format, including data from [Hernando-Herraez et al](https://www.nature.com/articles/s41467-019-12293-4), [Angermueller et al.](https://www.nature.com/articles/nmeth.3728) and [Smallwood et al.](https://www.nature.com/articles/nmeth.3035). The sequencing data from the [Gravina et al.](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-016-1011-3) study was also processed via the Bismark pipeline. Processed methylation data from the [Argelaguet et al.](https://www.nature.com/articles/s41586-019-1825-8) 
contained the same information as .cov files, but the columns were labeled and in a different order. These files were handled separately to produce
the required format for scAge, described below <br><br>

Ultimately, scAge requires a .tsv file or pandas datraframe with two columns, shown below. <br><br>
`ChrPos` should be in the form chr_position (i.e. chr1_3037802), while `MetLev` should be binary (0 or 1).

ChrPos | MetLev
:---: | :---:
chr1_3037802 | 1
chr19_61305429 | 0
... | ...

To load and process .cov files generated by Bismark, run:

```
process_coverage(cov_directory, 
                 n_cores = 10,
                 max_met = 100,
                 split = ".",
                 chunksize = 1,
                 binarization = "round",
                 output_path = None):
```

where
* `cov_directory` --> path to the directory where .cov single-cell methylation files are stored <br>
* `n_cores` --> number of cores to use for parallel processing <br>
* `maxmet` --> maximum methylation value (normally, methylation ratios of Bismark-generated files range from 0 to 100) <br>
* `split` --> desired string to split the file name on for single-cell name generation (if split is ".", then "SRR3136624.cov" becomes "SRR3136624") <br>
* `chunksize` --> number of coverage files that will be fed into a single worker process at a time <br>
* `binarization` --> choice of `round` vs. `discard`
                   Both methods involve dropping methylation values of 0.5
                   `round` rounds remaining non-binary values to 0 or 1 (this is the default), while
                   `discard` discards remaining non-binary values <br>
* `output_path` --> the full path to the output directory in which processed .tsv binary matrices should be written in
If `output_path` is set to `None`, named binary methylation matrices are returned in the form of a dictionary of pandas DataFrames <br>

### Predicting epigenetic age
The core of scAge is run_scAge, a function that enables epigenetic age predictions from a previously processed set of binarized single-cell
methylome profiles and a reference regression matrix generated from bulk data (ideally from the same tissue). To run this function:

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

where
* `single_cell_dir_or_dict` --> full path to a directory OR a dictionary of processed single-cell methylomes, created with `process_coverage`
* `single_cell_set_name` --> desired name of the scAge run (i.e. 'dataset_x')
* `reference_data` --> full path to the reference regression matrix created with `create_reference`
* `selection_mode` --> one of `percentile` (default), `numCpGs`, `cutoff`, which determines which CpG selection mode should be used in the algorithm.
I recommend `percentile` mode, which selects the most highly-correlated CpGs while account for differential coverage between cells. Alternatively,
a defined number of CpGs can be chosen using `numCpGs` mode, or a cutoff based on Pearson correlation can be set with `cutoff` 
(where only CpGs with |r| ≥ threshold are selected).
* `CpG_parameter` --> parameter accompanying selection mode 
(1 in percentile mode --> top 1% age-correlated CpGs) <br>
(500 in numCpGs mode --> top 500 age-correlated CpGs) <br>
(0.7 in cutoff mode --> top age-correlated CpGs with |r| ≥ 0.7) <br>
* `zero_met_replacement1` --> if the linear model built on bulk data goes below 0 for a given age, this value replaces the probability
* `one_met_replacement` --> if the linear model built on bulk data goes below 0 for a given age, this value replaces the probability
* `min_age` --> minimum age for probability computations
* `max_age` --> maximum age for probability computations
* `age_step` --> step for probability computations (i.e. if age_step == 1, likelihoods will be calculated for every 1 month between `min_age` and `max_age`)
* `n_cores` --> number of cores that should be used for parallelization
* `uncertainty` --> the uncertainty metric that should be used to compute upper and lower bounds (higher value results in a wider interval)
* `output_path` --> full path to the directory where predictions and the report file should be written to
* `chunksize` --> number of single-cell methylomes that should be passed to each parallel worker at once

The output of run_scAge is a .tsv matrix, which contains a number of columns detailing the internal computations of the algorithm.
The leftmost columns are the most critical, including: 
* `PredictedAge` --> scDNAm epigenetic age predictions in the unit of the training data (in our case, months)
* `MeanMet` --> mean global methylation of the single cell (which we show is related to epigenetic age in some cell types)
* `CellCoverage` --> number of CpGs covered in the single cell (which can be used to filter predictions)

An example output file is shown below:

Cell | PredictedAge | MeanMet | CellCoverage | ... |
:---: | :---: | :---: | :---: | :---: |
SRR3136627 | 0.5 | 0.663169 | 3914949 | ... |
SRR3136659 | 4.0 | 0.683454 | 799350 | ... |
SRR3136628 | 25.0 | 0.695256 | 2511084 | ... |

## Troubleshooting
If you encounter any issue when trying to run scAge, please feel free to contact me by email: alexandre.trapp1@gmail.com

## Information and acknowledgments
This software was developed by Alexandre Trapp, Technical Research Assistant in the Gladyshev Lab 
at Harvard Medical School and Brigham and Women's Hospital. I want to acknowledge all the members
of the Gladyshev Lab for their input, particularly Csaba Kerepesi and Vadim Gladyshev for their
contributions to the project, as well as Tiamat Fox and Adit Ganguly for their help with schematic
design.

This algorithm is covered by a provisional patent application filed by Brigham and Women's Hospital
on 3/12/2021 which names Alexandre Trapp, Csaba Kerepesi, and Vadim Gladyshev as inventors.
