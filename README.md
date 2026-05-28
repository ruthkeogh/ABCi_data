# ABCi: A benchmark causal inference data set

This repository provides the ABCi data set, data generation code, and code for illustrating the use of the data, as described in the draft paper (See ABCi_data_paper_draft_August2025.pdf for the current draft):

### "A benchmark causal inference (ABCi) data set: simulated data on time-varying treatments, confounders and outcomes based on patients with type 2 diabetes." Ruth Keogh, Nan van Geloven, Daniala Weir.

# The data
The ABCi data set is provided as an R data file, data_ABCi.RData. 
The code for generating the data is in sim_ABCi_data_generate.R. 
Descriptive information about the data can be generated using data_descriptives.R.

# Illustration
An illustration of the use of the data, as described in the draft paper for Example Question 1 (Table 2). The files for performing data set-up and analyses, and plotting and tabulating results, are in the folder Q1_files. All other files are called from the file _master_Q1.R. 

Similar code for Example Questions 2,3,4 appears in the folders Q2_files, Q3_files, Q4_files. 


