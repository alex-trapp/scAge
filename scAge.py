'''
scAge v1.0 (7/28/2021)

scAge is a flexible framework for epigenetic age predictions in single cells.
For more information on the algorithm, please consult Trapp et al., bioRxiv (2021).
If you use this software, we ask that you please cite our work.

The scAge pipeline consists of three key steps:
    1) Computing reference linear models for each CpG within a bulk DNAm matrix
    2) Loading in and processing single-cell methylation coverage files
    3) Predicting age of single cells from processed binary methylation matrices
   
and can be executed with the following functions:
    1) construct_reference
    2) process_coverage
    3) run_scAge

For more details on the algorithm and how to run the functions, 
please consult the GitHub page @ https://github.com/alex-trapp/scAge/

Developed by Alexandre Trapp in the Gladyshev Lab
Copyright, Alexandre Trapp, 2021
'''

# import packages and check dependencies
try:
    import numpy as np
    import pandas as pd
    import os
    import time
    from datetime import datetime
    import scipy.stats as ss
    from multiprocessing import Pool
    from sklearn.linear_model import LinearRegression
    import tqdm
    from tqdm.contrib.concurrent import process_map
    import warnings
    warnings.filterwarnings("ignore")
except:
    print("One or more required packages is not installed.")
    print("Please verify dependencies and try again.")
    exit()

def commas(value):
    '''
    Summary:
    ----------
    This function formats an integer or a float into a comma-separated string
    (i.e. 1000 -> '1,000')
    
    Parameters
    ----------
    value : int or float value
    
    Returns
    ----------
    value_w_comma : string of comma-separated number (i.e. 1000 -> '1,000')
    '''
    value_w_comma = f'{value:,}' 
    return value_w_comma

def get_range(value_list):
    '''
    Summary:
    ----------
    This functions return the range (minimum value, maximum value) of a given list or array
    
    Parameters
    ----------
    value_list: list of floats or ints
    
    Returns
    ----------
    min_max_tuple: tuple of the form (min_value, max_value)
    '''
    min_max_tuple = (min(value_list), max(value_list))
    return min_max_tuple

def closest(lst, K):
    '''
    Summary:
    ----------
    This function returns the closest item to "K" in a given list "lst"
    
    Parameters
    ----------
    lst: list of numerical values
    K: number to query within list
    
    Returns
    ----------
    closest_value: the closest value (absolute) from K within lst
    '''
    closest_value = lst[min(range(len(lst)), key = lambda i: abs(lst[i]-K))]
    return closest_value

def load_cov_file(args):
    '''
    Summary:
    ----------
    This function acts as the workhorse for the parallelization function process_coverage.
    It takes as input a tuple of arguments from process_coverage and either returns a 
    processed methylation matrix or writes it to a specified output path.
    
    Refer to process_coverage for more information.
    
    Parameters:
    ----------
    args: a tuple of arguments supplied in parallelization function load_cov
          in the form (file, cov_directory, max_met, split, binarization, write_path)
              file -- the name of the .cov/.cov.gz file
              cov_directory -- the path to the directory containing coverage files
              max_met -- the maximum methylation rate in the .cov file (1 or 100, depending on processing)
              split -- an optional argument dictating how the file name should be split to name cells
                       i.e. if split is ".", then "SRR3136624.cov" becomes SRR3136624 
              binarization -- "round" or "discard": both remove methylation values of 0.5,
                               --> "round" rounds remanining values to 0 or 1
                               --> "discard" removes remaining values that are not 0 or 1
              write_path -- the full path of the directory in which to store processed .tsv files
    Returns:
    ----------
    if write_path == None:
        (cell, cov): a tuple containing
                    cell -- the name/identifier of a cell
                    cov -- the processed binary methylation matrix 
                           (col 1 = genomic position, col 2 = binary methylation)
    else:
        cov (the processed binary methylation matrix) is written as a .tsv file in write_path
        
    '''
    # load arguments
    file = args[0]
    cov_directory = args[1]
    max_met = args[2]
    split = args[3]
    binarization = args[4]
    write_path = args[5]

    # list of autosomes for filtering
    autosome_list = [str(x) for x in range(1, 20)]
    
    # split cell name from file name with desired split
    cell = file.split(split)[0]
    
    # read in coverage file
    cov = pd.read_csv(cov_directory + file,
                      sep = "\t", header = None,
                      names = ["Chr", "Pos1", "Pos2", "MetLev", "Met", "Unmet"],
                      dtype={'Chr': 'str',
                             'Pos1' : 'str'})
    
    if "chr" in cov.iloc[0, 0]: # check for if column is labeled "chr15" instead of "15"
        cov['Chr'] = cov['Chr'].str.replace('chr', '')
    
    
    # filter autosomes
    cov = cov[cov["Chr"].isin(autosome_list)]
    
    # create ChrPos column
    cov["ChrPos"] = "chr" + cov["Chr"] + "_" + cov["Pos1"]
    
    # set index to ChrPos
    cov = cov.set_index("ChrPos")
    
    # sort genomic positions
    cov[["Chr", "Pos1"]] = cov[["Chr", "Pos1"]].astype("int")
    cov = cov.sort_values(["Chr", "Pos1"])
    
    # check maximum methylation
    if max_met == 100: 
        # methylation must be adjusted to [0, 1] for scAge
        cov['MetLev'] = cov["MetLev"] / 100
    else:
        pass
    
    # remove methylation levels of 0.5
    cov = cov[cov['MetLev'] != 0.5]
    
    # binarization method:
    # 'round' refers to rounding remaining values to 0 or 1
    # 'discard' refers to removing values that are not 0 or 1
    if binarization == "round":
        cov["MetLev"] = cov["MetLev"].round().astype("int8")
    elif binarization == "discard":
        cov = cov[cov["MetLev"].isin([0, 1])]
        cov["MetLev"] = cov["MetLev"].astype("int8")
        
    # minimizing the processed dataframe to only necessary columns
    cov = cov.drop(["Chr", "Pos1", "Pos2", "Met", "Unmet"], axis = 1)
    
    # remove duplicate indeces if there are any (there should not be any)
    cov = cov[~cov.index.duplicated(keep='first')]
    
    # determining whether to return a tuple or write data to a .tsv
    if write_path == None:
        return (cell, cov) # return tuple of cell name and coverage dataframe
    else:
        # write to file
        if binarization == "round":  
            cov.to_csv(write_path + cell + ".tsv", sep = "\t")
        elif binarization == "discard":  
            cov.to_csv(write_path + cell + "-dis" + ".tsv", sep = "\t")
    del cov

def process_coverage(cov_directory,
                     output_path = "./sc_data_processed/",
                     n_cores = 1,
                     max_met = 100,
                     split = ".",
                     chunksize = 1,
                     binarization = "round"):
    
    '''
    Summary:
    ----------
    This function is a parallelization tool that internally relies on load_cov_file. It uses
    multi-core processing to load .cov or .cov.gz files and process them into what scAge
    needs as input (binarized methylation matrices). This function either returns
    a dictionary of single-cell methylome matrices, or writes them to a specified directory.
    
    Parameters:
    ----------
    cov_directory: str, the path to the directory containing .cov/.cov.gz files
    n_cores: int, the number of cores to run in parallel
    max_met: int, the maximum methylation in coverage files (usually 100, sometimes 1)
    split: str, the symbol/letter/number to split by when generating single-cell names from files
    chunksize: int, number of elements to feed to each worker during parallel processing
    binarization: "round" or "discard". 
                  Both methods involve dropping methylation values of 0.5.
                  "round" rounds remaining non-binary values to 0 or 1 (default)
                  "discard" discards remaining non-binary values
    output_path: str, the path to the output directory in which to write .tsv files
              
    Returns:
    ----------
    if output_path == None:
        sc_dict: dict, with cell names/identifiers as keys and 
                 binary methylation matrices as values
    else:
        binary methylation matrices are written to .tsv files inside of output_path directory
        
    '''
    
    print("process_coverage function starting!\n")
    
    start_time = time.time()
    
    print("----------------------------------------------------------")
    # get list of files in directory
    print("Loading .cov files from '%s'" % cov_directory)
    cov_list = sorted(os.listdir(cov_directory))
    print("Number of Bismark .cov files = %s" % len(cov_list))
    for file in cov_list:
        if file == ".ipynb_checkpoints":
            cov_list.remove(file)
    print("First .cov file name: '%s'" % cov_list[0])
    print("----------------------------------------------------------\n")

    # create tuple of arguments for load_cov_file
    file_tuples = []
    for file in cov_list:
        file_tuples.append((file, cov_directory, max_met, split, binarization, output_path))
    
    print("----------------------------------------------------------")
    print("Starting parallel loading and processing of .cov files...")
    # parallelization function with tqdm progress bar
    results = process_map(load_cov_file,
                          file_tuples,
                          max_workers = n_cores,
                          chunksize = chunksize,
                          desc = "Single-cell loading progress ",
                          unit = " cell methylomes")
    print("\nParallel loading complete!")
    
    # determine whether to write data or return a dictionary
    if output_path == None: # return dictionary of cell binary matrices
        sc_dict = {}
        for result_tuple in results:
            sc_dict[result_tuple[0]] = result_tuple[1]
        print("Returning a dictionary, as no output path was given.")
        return sc_dict
    else: # cell matrices are written to output_path inside of load_cov_file
        print("Processed binary methylation matrices written to '%s'" % output_path)    
    print("----------------------------------------------------------\n") 
          
    end_time = time.time()
    elapsed_time = end_time - start_time
    print("Time elapsed to process coverage files = %.3f seconds" % elapsed_time)
    print("\nprocess_coverage run complete!")

def compute_linear_relationship(args):
    
    '''
    Summary:
    ----------
    This function computes correlation and regression metrics for a single CpG from a
    bulk reference dataset. It is integrated into construct_reference, which encapsulates
    it for efficient multi-core parallel processing.
    
    Refer to the documentation for construct_reference for additional information.
    
    Parameters:
    ----------
    args: tuple, of the form (CpG, metlev_series, age_series) where
          CpG -- str, genomic position of the CpG, in the form chr5_1234567
          metlev_series -- pandas Series, bulk methylation levels in the range [0, 1]
          age_series -- pandas Series, bulk sample ages
    
    Note: the unit in which scAge is trained (days, weeks, months) is the unit in
    which the epigenetic age prediction outputs will be.
              
    Returns:
    ----------
    ref_output: tuple, of the form (CpG, Pearson_r, Pearson_p, coef, intercept) where
                    CpG: str, genomic position of the CpG
                    Pearson_r: float, Pearson correlation coefficient of methylation level and age
                    Pearson_p: float, p-value associated with the Pearson correlation coefficient
                    coef: float, linear regression coefficient
                    intercept: float, linear regression intercept
    Notes: 
        Pearson metrics are calculated with scipy.stats.pearsonr
        Regression metrics are calculated with sklearn.linear_model.LinearRegression
        Missing values are dropped in both computations
        
    '''
    
    # load arguments
    CpG = args[0]
    metlev_series = args[1]
    age_series = args[2]
    
    # isolate valid values (i.e. metlev/age non-NaN pairs)
    valid_indeces = metlev_series.isna()
    valid_metlevs = np.array(metlev_series[~valid_indeces])
    valid_ages = np.array(age_series[~valid_indeces])
    
    # get Pearson correlation metrics
    Pearson_results = ss.pearsonr(valid_metlevs, valid_ages)
    Pearson_r = Pearson_results[0]
    Pearson_p = Pearson_results[1]

    # calculate linear regression (with traditional ordinary least squares (OLS))
    reg = LinearRegression(n_jobs = 1).fit(valid_ages.reshape(-1, 1),
                                           valid_metlevs)
    
    # isolate coefficient (slope) and intercept
    coef = reg.coef_[0]
    intercept = reg.intercept_
    
    # return output tuple
    ref_output = (CpG, Pearson_r, Pearson_p, coef, intercept)
    return ref_output

def construct_reference(training_DNAm_matrix,
                        output_path,
                        n_cores = 1,
                        chunksize = 100):
    
    '''
    Summary:
    ----------
    This function parallelizes the worker function compute_linear_relationship.
    It takes as input a training matrix of methylation values for bulk samples 
    of different ages, and computes linear associations between age and methylation
    for all CpGs in the matrix. Methylation values must range between 0 or 1, and the
    matrix must be have samples as rows and CpGs as columns. Addditionally, there must
    be at least one numerical column labeled "Age", which is used for computing
    correlation and regression metrics.
    
    Parameters:
    ----------
    training_DNAm_matrix: pandas dataframe of samples (rows) and CpG sites (columns)
                          with some additional metadata columns (at least "Age")
    output_path: the full path to the desired output reference file
    n_cores: the number of cores to use for parallel processing
    chunksize: the number of elements to feed to each worker at once
              
    Returns:
    ----------
    No output, the reference matrix is simply saved to output_path as a .tsv file
    
    Notes: 
        Pearson metrics are calculated with scipy.stats.pearsonr
        Regression metrics are calculated with sklearn.linear_model.LinearRegression
        Missing values are dropped in both computations
        
    '''
    
    start_time = time.time()
    print("construct_reference function starting!\n")
    
    # get age data
    age_series = training_DNAm_matrix.loc[:, "Age"]

    # get samples
    bulk_sample_names = list(training_DNAm_matrix.index)
    
    # get CpGs
    bulk_CpG_names = [x for x in training_DNAm_matrix.columns if "chr" in x]
    
    print("----------------------------------------------------------")
    print("Number of samples = %s" % len(bulk_sample_names))
    print("Number of CpGs = %s" % commas(len(bulk_CpG_names)))
    print("----------------------------------------------------------\n\n")
    
    # construct list of arguments for process_map parallel function
    print("----------------------------------------------------------")
    print("Constructing list of arguments for parallel processing...")
    list_of_arguments_linear = []
    
    for CpG in tqdm.auto.tqdm(bulk_CpG_names,
                              desc = "Reference progress (1/2) ",
                              unit = " CpGs"):
        metlev_series = training_DNAm_matrix[CpG]
        list_of_arguments_linear.append((CpG, metlev_series, age_series))
    print("Argument list constructed!")
    print("----------------------------------------------------------\n\n")
    
    # parallel processing with a tqdm progress bar
    print("----------------------------------------------------------")
    print("Starting parallel processing with %s cores..." % n_cores) 
    results = process_map(compute_linear_relationship,
                          list_of_arguments_linear,
                          max_workers = n_cores,
                          chunksize = chunksize,
                          desc = "Reference progress (2/2) ",
                          unit = " CpGs")
    
    # get reference matrix
    results_df = pd.DataFrame(results,
                              columns = ["ChrPos",
                                         "PearsonR",
                                         "PearsonP",
                                         "Coef",
                                         "Intercept"]).set_index("ChrPos")
    
    # remove missing values (correlations of NaN, i.e. CpGs that do not change with age)
    results_df = results_df.dropna(axis = 0)
    
    # write reference to a .tsv file
    results_df.to_csv(output_path, sep = "\t")
    print("\nReference model dataset written to '%s'" % output_path)

    # write report file detailing input matrix and distributions
    output_file_name = output_path.split("/")[-1][:-4]
    with open('%s.report.txt' % output_path[:-4], 'w') as writer:
        writer.write("scAge reference report for %s\n" % output_file_name)
        now = datetime.now()
        write_datetime = now.strftime("%m/%d/%Y %H:%M:%S")
        writer.write("Reference file created: %s\n\n" % write_datetime)
        writer.write("Number of input samples = %s\n" % len(bulk_sample_names))
        writer.write("Number of input CpGs = %s\n" % commas(len(bulk_CpG_names)))
        writer.write("Number of output CpGs (after dropping NA) = %s\n\n" % commas(len(results_df)))
        
        # if some metadata was included, it is used
        # currently, the function supports 'Tissue', 'Strain', and 'Gender',
        # which are the metadata provided in the Thompson et al. (2018) study
        if "Tissue" in training_DNAm_matrix.columns:
            writer.write("Tissue(s)\n%s\n\n" % training_DNAm_matrix["Tissue"].value_counts())
        if "Strain" in training_DNAm_matrix.columns:
            writer.write("Strain(s)\n%s\n\n" % training_DNAm_matrix["Strain"].value_counts())
        if "Gender" in training_DNAm_matrix.columns:
            writer.write("Sex\n%s\n\n" % training_DNAm_matrix["Gender"].value_counts())
            
    print("Report file generated at '%s.report.txt'" % output_path[:-4])
    print("----------------------------------------------------------\n\n")
    
    end_time = time.time()
    elapsed_time = end_time - start_time
    print("\nTime to run construct_reference: %0.3f seconds" % elapsed_time)
    
    print("\nconstruct_reference run complete!")
    
def compute_probabilities(args):
    '''
    Summary:
    ----------
    This function is the core epigenetic age profiling worker tool of the scAge framework.
    It takes as input a tuple of arguments from the parallelization function run_scAge,
    and outputs a tuple containing a variety of information regarding the
    single-cell methylome profile, most importantly the predicted epigenetic age.
    
    Refer to run_scAge for additional information.
    
    Parameters:
    ----------
    args: tuple, in the form (single_cell_name, single_cell_met, corr_regr_data, selection_mode,
                              CpG_parameter, zero_met_replacement, one_met_replacement,
                              min_age, max_age, age_step, uncertainty), where
          single_cell_name: str, the name of the cell
          single_cell_met: str or pandas dataframe, the processed binary methylation matrix
                           You can supply either the input path of the file as a string
                           or directly a pandas dataframe structure of the matrix
          reference_data: str, the full file path to the desired reference data
          selection_mode: str, one of [numCpGs, percentile, cutoff] where
                          numCpGs -- selects a defined number of age-associated CpGs per cell
                          percentile -- selects the top x% age-associated CpGs per cell
                          cutoff -- selects only CpGs with a Pearson correlation ≥ cutoff
          CpG_parameter: float, the parameter to feed in given a specific selection mode
                         ex1: selection_mode == numCpGs --> CpG_parameter = 1000 (1000 CpGs/cell)
                         ex2: selection_mode == percentile --> CpG_parameter = 1 (Top 1% CpG)
                         ex3: selection_mode == cutoff --> CpG_parameter = 0.7 (Only CpGs with r ≥ 0.7)
          zero_met_replacement: float, the lower bound (for when met ≤ 0) of bulk methylation level
                                predictions based on the linear models generated in construct_reference
          one_met_replacement: float, the upper bound (for when met ≥ 1) of bulk methylation level
                                predictions based on the linear models generated in construct_reference
          min_age: float, the minimum age for which to build a probability profile
          max_age: float, the maximum age for which to build a probability profile
          age_step: float, the step value for computing probability profiles
                    (i.e. if age_step == 1, likelihoods will be calculated for every 1 month)
          uncertainty: float, the width of the uncertainty metric (+/- uncertainty) to provide
                       a confidence interval of predictions
              
    Returns:
    ----------
    probabilities_output: tuple, in the form 
                          (single_cell_name, max_probability_age,
                          age_probability_df, corr_regr_singlecell_df,
                          mean_met, coverage, num_intersections,
                          lower_bound_age, upper_bound_age) where:
                          
                             single_cell_name: str, the name of the cell
                             max_probability_age: float, maximum likelihood age estimate
                             age_probability_df: pandas DataFrame, with log-likelihoods for each age
                             corr_regr_singlecell_df: pandas DataFrame, with information regarding
                                                      selected CpGs (Pearson correlations,
                                                      regression coefficient, regression intercepts,
                                                      and binary methylation levels)
                             mean_met: float, mean global methylation of all binary states in the cell
                             coverage: int, the number of individual CpGs covered (either strand)
                             num_intersections: int, the number of CpGs that intersect with a given
                                                reference training set
                             lower_bound_age: float, the lower age bound of the uncertainty estimate
                             upper_bound_age: float, the upper age bound of the uncertainty estimate
    
    '''
    
    # get arguments
    single_cell_name = args[0]
    single_cell_met = args[1]
    corr_regr_data = args[2]
    selection_mode = args[3]
    CpG_parameter = args[4]
    zero_met_replacement = args[5]
    one_met_replacement = args[6]
    min_age = args[7]
    max_age = args[8]
    age_step = args[9]
    uncertainty = args[10]
    
    # determine whether input is pandas dataframe or string (file path)
    input_type = str(type(single_cell_met))
    if input_type == "<class 'str'>":
        single_cell_met = pd.read_csv(single_cell_met, sep = "\t", index_col = 0)
    elif input_type == "<class 'pandas.core.frame.DataFrame'>":
        pass
        
    start = time.time()
    
    # intersecting CpGs between bulk and reference dataset
    
    ref_sc_intersect_df = pd.concat([corr_regr_data, single_cell_met], axis = 1,
                                     join = "inner")

    # profiling mode selection
    if selection_mode == "percentile": # ex: top 1% age-associated CpGs per cell
        quantile = ref_sc_intersect_df['PearsonR'].abs().quantile(q = CpG_parameter)
        PearsonR_abs_top = ref_sc_intersect_df[ref_sc_intersect_df["PearsonR"].abs() >= quantile]
    elif selection_mode == "numCpGs": # ex: top 500 age-associated CpGs per cell
        PearsonR_abs_top = ref_sc_intersect_df['PearsonR'].abs().nlargest(CpG_parameter)
    elif selection_mode == "cutoff": # ex: CpGs above a cutoff of r ≥ 0.7 per cell
        PearsonR_abs_top = ref_sc_intersect_df[ref_sc_intersect_df["PearsonR"].abs() >= CpG_parameter]
    
    # subset dataframe to chosen highly age-correlated CpGs
    ref_sc_intersect_df_subset = ref_sc_intersect_df.loc[PearsonR_abs_top.index, :]
    
    # isolate selected sites
    selected_sites = list(ref_sc_intersect_df_subset.index)
    
    # get age steps from min_age to max_age (inclusive of both)
    age_steps = np.arange(min_age, max_age + age_step, age_step)

    # create list to store probability profiles
    list_of_profile_probabilities_per_age = []
    
    # loop through each age step
    for age in age_steps:
        # create list to hold probability for all chosen CpGs for a given age
        probability_list_one_age = []
        # loop through each site
        for site in selected_sites:
            # isolate slope
            slope = ref_sc_intersect_df_subset.loc[site, "Coef"]
            # isolate intercept
            intercept = ref_sc_intersect_df_subset.loc[site, "Intercept"]
            # compute methylation probability from reference data
            # using slope, intercept, and the current age-step
            methylation_probability = slope * age + intercept # compute methylation probability

            # methylation must be inherently bounded between 0 (fully unmethylated)
            # and 1 (fully methyalted). Hence:
            # if methylation_probability is 1 or above (based on linear regression)
            # replace with "one_met_replacement"
            if methylation_probability >= 1: 
                methylation_probability = one_met_replacement
                
            # if methylation_probability is 0 or below (based on linear regression)
            # replace with "one_met_replacement"
            elif methylation_probability <= 0: 
                methylation_probability = zero_met_replacement
                
            # in most cases, methylation probability stays untouched
            else: 
                methylation_probability = methylation_probability

            # get single cell binary methylation level
            site_methylation_sc = ref_sc_intersect_df_subset.loc[site, "MetLev"] 
            
            # if the CpG is methylated, append ln(methylation_probability)
            if site_methylation_sc == 1:
                probability_list_one_age.append(np.log(methylation_probability))
                
            # if the CpG is methylated, append ln(1 - methylation_probability)
            elif site_methylation_sc == 0:
                probability_list_one_age.append(np.log(1 - methylation_probability))
                
            else:
                raise ValueError("Encountered a non-binary methylation state!\n" + \
                                 "Please verify input single-cell data is properly" + \
                                 "binarized before running run_scAge.")
                
        # compute log-likelihood sum and appened to list
        # this is equivalent to computing the product of probabilities
        # but neatly avoid underflow errors that result when multiplying
        # many small numbers together
        list_of_profile_probabilities_per_age.append(np.sum(probability_list_one_age)) 
        
    # transform into dataframe with age steps
    age_probability_df = pd.DataFrame({"Pr" : list_of_profile_probabilities_per_age},
                                      index = age_steps)
    
    # compute highest likelihood age among age steps
    max_probability_age = round(float(age_probability_df.idxmax()), 2) 
    
    # compute maximum probability
    max_probability = float(age_probability_df['Pr'].max())
    
    # get likelihood based on uncertainty parameter
    likelihood_uncertainty_bound = max_probability - uncertainty
    
    # isolate probability curve below max point
    age_pred_df_below = age_probability_df[age_probability_df.index.astype("float") < max_probability_age]
    
    # isolate probability curve above max point
    age_pred_df_above = age_probability_df[age_probability_df.index.astype("float") > max_probability_age]
    
    # get age bounds for uncertainty confidence interval
    try:
        # get closest log-likelihood below the maximum based on provided uncertainty parameter
        closest_num_below = closest(list(age_pred_df_below['Pr']),
                                         likelihood_uncertainty_bound)
        lower_bound_age = round(float(age_pred_df_below[age_pred_df_below["Pr"] \
                                                        == closest_num_below].index.values), 1)

        # get closest log-likelihood above the maximum based on provided uncertainty parameter
        closest_num_above = closest(list(age_pred_df_above['Pr']),
                                        likelihood_uncertainty_bound)
        upper_bound_age = round(float(age_pred_df_above[age_pred_df_above["Pr"] \
                                                        == closest_num_above].index.values), 1)
        
    # if this throws an error, return NaN
    except:
        lower_bound_age = np.nan
        upper_bound_age = np.nan
    
    end = time.time()
    
    # compute single-cell characteristics
    mean_met = single_cell_met["MetLev"].mean()
    coverage = len(single_cell_met)
    num_intersections = len(ref_sc_intersect_df)
    
    # return tuple output
    probabilities_output = (single_cell_name, max_probability_age,
                             age_probability_df, ref_sc_intersect_df_subset,
                             mean_met, coverage, num_intersections,
                             lower_bound_age, upper_bound_age)
    
    return probabilities_output

def run_scAge(single_cell_dir_or_dict,
              single_cell_set_name,
              reference_data,
              output_path = "./predictions/",
              selection_mode = "percentile",
              CpG_parameter = 1,
              zero_met_replacement = 0.001,
              one_met_replacement = 0.999,
              min_age = -20,
              max_age = 60,
              age_step = 0.1,
              n_cores = 1,
              uncertainty = 1,
              chunksize = 5):
    
    '''
    Summary:
    ----------
    This is the main function of the pipeline. It parallelizes compute_probabilities,
    enabling rapid and scalable estimation of epigenetic age across many cells simultaneously.
    It takes as input a single-cell profile and a reference matrix, and returns
    a dataframe with predicted epigenetic age (scDNAm age), as well as abundant 
    information regarding single-cell characteristics and age-correlated CpGs 
    chosen as part of the algorithm. These additional columns are useful to conduct
    downstream analyses on covariates and other techincal or biological factors.
    
    Parameters:
    ----------
    single_cell_dir_or_dict: str or dict, either the directory containing processed methylation
                             data as .tsv/.tsv.gz files (i.e. generated by process_coverage),
                             or a dictionary of labeled methylation matrices
    single_cell_set_name: str, the desired name of the single cell data
                          this is used in setting the name of the output file
    reference_data: str, the full file path to the desired reference data
    selection_mode: str, one of [numCpGs, percentile, cutoff] where
                         percentile -- selects the top x% age-associated CpGs per cell
                         numCpGs -- selects a defined number of age-associated CpGs per cell
                         cutoff -- selects only CpGs with a Pearson correlation ≥ cutoff
    CpG_parameter: float, the parameter to feed in given a specific selection mode
                          ex1: selection_mode == percentile --> CpG_parameter = 1 (Top 1% percentile)
                          ex2: selection_mode == numCpGs --> CpG_parameter = 1000 (1000 CpGs/cell)
                          ex3: selection_mode == cutoff --> CpG_parameter = 0.7 (Only CpGs with r ≥ 0.7)
    zero_met_replacement: float, the lower bound (methylation level ≤ 0) of bulk methylation level
                          predictions based on the linear models generated in construct_reference
    one_met_replacement: float, the upper bound (methylation level ≥ 1) of bulk methylation level
                         predictions based on the linear models generated in construct_reference
    min_age: float, the minimum age for which to build a probability profile
    max_age: float, the maximum age for which to build a probability profile
    age_step: float, the step value for computing probability profiles
              (i.e. if age_step == 1, likelihoods will be calculated for every 1 month)
    n_cores: int, the number of cores to use for parallel processing
    uncertainty: float, the width of the uncertainty metric (+/- uncertainty) to provide
                 a confidence interval of predictions
              
    Returns:
    ----------
    2 files are created: 1 .report.txt file containing the parameters used in running the algorithm
                         and 1 .tsv file containing the results and predictions for each single cell
    
    '''
    
    start = time.time()
    print("scAge algorithm starting!\n")
    
    print("----------------------------------------------------------")
    print("Profiling epigenetic age in '%s' single-cell data..." % single_cell_set_name)
    
    # determine whether input data is a path to a directory or a dictionary of data
    input_type = str(type(single_cell_dir_or_dict))
    if input_type == "<class 'str'>":
        single_cell_cov_dir = single_cell_dir_or_dict
        print("Loading processed binary methylation files from '%s'..." % single_cell_cov_dir)
        single_cell_files = sorted(os.listdir(single_cell_cov_dir))
        for file in single_cell_files:     
            if file == ".ipynb_checkpoints":
                single_cell_files.remove(file)
        # check if files are gzipped
        if ".gz" in single_cell_files[0]:
            add_gz = True
        else:
            add_gz = False
        single_cells = [cell.split(".tsv")[0] for cell in single_cell_files]
    elif input_type == "<class 'dict'>":
        print("Using cell stored in dictionary...")
        single_cells = list(single_cell_dir_or_dict.keys())

    print("Number of single cells to analyze: %s" % len(single_cells))
    print("----------------------------------------------------------")
    
    # get name of the reference dataset from the path input
    training_dataset_name = reference_data.split("/")[-1].split(".tsv")[0]
    
    print("\nscAge parameters:")
    print("----------------------------------------------------------")
    print("Using reference training data: %s" % training_dataset_name)
    
    # if the reference file is there, loads it in
    try:
        corr_regr_data = pd.read_csv(reference_data,
                                     sep = "\t",
                                     index_col = 0)
        
    # if the reference file cannot be found, an error is thrown
    except:
        raise NameError("Reference training set not found, please verify input directory.")
        
    print("Shape of reference matrix: {} CpGs, {} metric columns".format(commas(corr_regr_data.shape[0]),
                                                                                corr_regr_data.shape[1]))
    print("\n")
    print("Using %s cores with chunksize of %s" % (n_cores, chunksize))
    print("\n")
    print("Setting minimum age to %s month(s)" % min_age)
    print("Setting maximum age to %s month(s)" % max_age)
    print("Using age step of %s month(s)" % age_step)
    print("\n")
    print("Replacing modeled bulk methylation ≤ 0 with %s" % zero_met_replacement)
    print("Replacing modeled bulk methylation ≥ 1 with %s" % one_met_replacement)
    print("\n")
    print("Using profiling mode: %s" % selection_mode)
    if selection_mode == "percentile":
        print("--> Profiling top %s%% age-related CpGs by absolute Pearson correlation " % \
              (str(CpG_parameter)))
        # for example, providing a value of 1 means selecting the top 1% (absolute highest)
        # age-correlated CpGs
        CpG_parameter_num = (100 - CpG_parameter) / 100
    elif selection_mode == "numCpGs":
        # for example, providing a value of 1000 means selecting the top 1000 (absolute highest)
        # age-correlated CpGs
        print("--> Profiling top %s age-related CpGs by absolute Pearson correlation" % CpG_parameter)
        CpG_parameter_num = CpG_parameter
    elif selection_mode == "cutoff":
        # for example, providing a value of 0.7 means CpGs with an absolute age-correlation 
        # greater than or equal to 0.7 
        print("--> Profiling top age-related CpGs above an absolute correlation cutoff of %s" % CpG_parameter)
        CpG_parameter_num = CpG_parameter
    else:
        raise ValueError("Incorrect selection mode, must be one of ['percentile', 'numCpGs', 'cutoff']")
           
    print("\nUsing a prediction uncertainty metric of +/- %s " % uncertainty + \
          "for confidence interval computation")
    print("----------------------------------------------------------")
    
    # create tuple of arguments for parallel processing
    list_of_arguments_parallel_scAge = []
    
    # if the full path path to processed single cell methylation files is given
    if input_type == "<class 'str'>":
        for cell in single_cells:
            if add_gz == True:
                cell_path = single_cell_dir_or_dict + cell + ".tsv.gz"
            elif add_gz == False:
                cell_path = single_cell_dir_or_dict + cell + ".tsv"
            list_of_arguments_parallel_scAge.append((cell,
                                                     cell_path,
                                                     corr_regr_data,
                                                     selection_mode,
                                                     CpG_parameter_num,
                                                     zero_met_replacement,
                                                     one_met_replacement,
                                                     min_age,
                                                     max_age,
                                                     age_step,
                                                     uncertainty))
        
    # or if single-cell data is provided as a labeled dictionary
    elif input_type == "<class 'dict'>":
        for cell in single_cell_dir_or_dict:
            list_of_arguments_parallel_scAge.append((cell,
                                                     single_cell_dir_or_dict[cell],
                                                     corr_regr_data,
                                                     selection_mode,
                                                     CpG_parameter_num,
                                                     zero_met_replacement,
                                                     one_met_replacement,
                                                     min_age,
                                                     max_age,
                                                     age_step,
                                                     uncertainty))
    print("\n\n----------------------------------------------------------")
    print("Starting parallel processing of all cells with %s cores!\n" % n_cores)
        
    # compute probabilities using parallel processing
    # with a progress bar using tqdm
    results = process_map(compute_probabilities,
                          list_of_arguments_parallel_scAge,
                          max_workers = n_cores,
                          chunksize = chunksize,
                          desc = "scAge progress ",
                          unit = " age predictions")

    # process output data into a final dataframe
    cell_data_dict = {}
    for cell in results:
        # get name of the cell
        cell_name = cell[0]
        
        # get predicted age of the cell
        cell_age = cell[1]

        # get the list of age steps that were tested
        ages_tested = list(np.around(cell[2].index.values.astype("float"), 2))
        
        # get the likelihood for each age step
        ages_likelihoods = list(cell[2]['Pr'])
        
        # get CpGs chosen by the ranking algorithm
        CpGs_chosen = list(cell[3].index)
        
        # get the number of CpGs that were selected
        numCpGs_selected = len(CpGs_chosen)
        
        # isolate Pearson correlations of chosen CpGs
        correlations = list(cell[3]["PearsonR"])
        
        # isolate the linear regression coefficient of chosen CpGs
        regression_coefs = list(cell[3]["Coef"])
        
        # isolate the linear regression intercept of chosen CpGs
        regression_intercepts = list(cell[3]["Intercept"])
        
        # isolate the binary methylation value of chosen CpGs
        methylation_values = list(cell[3]["MetLev"])
        
        # get mean methylation of the cell
        mean_met = cell[4]
        
        # get CpG coverage of the cell
        coverage = cell[5]
        
        # get the number of intersections between single-cell and reference data
        num_intersections = cell[6]
        
        # get lower and upper bounds for probabilistic confidence interval
        lower_age_bound = cell[7]
        upper_age_bound = cell[8]

        # combine all the data into a list and save to dictionary
        cell_data_dict[cell_name] = [cell_age,
                                     mean_met,
                                     coverage,
                                     num_intersections,
                                     ages_tested,
                                     ages_likelihoods,
                                     CpGs_chosen,
                                     numCpGs_selected,
                                     correlations,
                                     regression_coefs,
                                     regression_intercepts,
                                     methylation_values,
                                     lower_age_bound,
                                     upper_age_bound]
        
    # create dataframe from dictionary
    cell_predictions_df = pd.DataFrame.from_dict(cell_data_dict,
                                                 columns = ["PredictedAge",
                                                            "MeanMet",
                                                            "CellCoverage",
                                                            "Intersections",
                                                            "AgesTested",
                                                            "AgeLikelihood",
                                                            "SelectedCpGs",
                                                            "NumberCpGs",
                                                            "Correlations",
                                                            "RegressionCoefs",
                                                            "RegressionIntercepts",
                                                            "MethylationValues",
                                                            "LowerBound", 
                                                            "UpperBound"],
                                                 orient = "index")
    
    # create descriptive name for output file
    # ex: name-train(Thompson_Liver_BL6)-mode(percentile)-param(top_1_pct).tsv
    # the most crucial parameters (the training data, the selection mode,
    # and the selection parameter) are automatically encoded in the output file name
    # additional data about the run is written to a .report.txt file
    base_output_name = output_path + single_cell_set_name + "-train(" + \
                       training_dataset_name + ")-mode(" + selection_mode
    
    if selection_mode == "percentile":
        output_file_name =  base_output_name + ")-param(top_%s_pct).tsv" % CpG_parameter
    if selection_mode == "numCpGs":
        output_file_name =  base_output_name + ")-param(%sCpGs).tsv" % CpG_parameter
    if selection_mode == "cutoff":
        output_file_name =  base_output_name + ")-param(above_%s_cutoff).tsv" % CpG_parameter
        
    # check if output path directory exists, otherwise creates it
    print("\n")
    try:
        os.listdir(output_path)
    except:
        print("Output path does not exist, creating directory...")
        os.makedirs(output_path)
        
    # save prediction dataframe to a .tsv file
    cell_predictions_df = cell_predictions_df.rename_axis(index='Cell')
    cell_predictions_df.to_csv(output_file_name,
                               sep = "\t")
    
    end = time.time()
    print("\nPredictions stored in '%s'" % output_path)
    print("----------------------------------------------------------")
    print("\nTime elapsed to generate scAge results = %0.3f seconds\n" % (end-start))
    print("scAge run complete!")
    
    # write a report file containing the parameters that were used
    # in the current scAge run
    with open('%s.report.txt' % output_file_name[:-4], 'w') as writer:
        writer.write("scAge report for '%s'\n" % output_file_name)
        now = datetime.now()
        write_datetime = now.strftime("%m/%d/%Y %H:%M:%S")
        writer.write("Files created: %s\n\n" % write_datetime)
        if input_type == "<class 'dict'>":
            writer.write("Methylation data loaded in from dictionary\n")
        if input_type == "<class 'str'>":
            writer.write("Methylation data loaded in from '%s'\n" % single_cell_dir_or_dict)
        writer.write("Training dataset used: %s\n" % training_dataset_name)
        writer.write("Selection mode: %s\n" % selection_mode)
        if selection_mode == "percentile":
            writer.write("CpG parameter: top %s%% age-associated CpGs\n" % CpG_parameter)
        elif selection_mode == "numCpGs":
            writer.write("CpG parameter: top %s age-associated CpGs\n" % commas(CpG_parameter))
        elif selection_mode == "cutoff":
            writer.write("CpG parameter: CpGs with Pearson correlation ≥ %s\n" % commas(CpG_parameter))
        writer.write("Minimum age: %s\n" % min_age)
        writer.write("Maximum age: %s\n" % max_age)
        writer.write("Age step: %s\n" % age_step)
        writer.write("Replacement for modeled bulk methylation below 0: %s\n" % zero_met_replacement)
        writer.write("Replacement for modeled bulk methylation above 1: %s\n" % one_met_replacement)
        writer.write("Uncertainty parameter for confidence interval: +/- %s \n\n" % uncertainty)
        writer.write("Time to generate results: %0.3f seconds" % (end-start))
