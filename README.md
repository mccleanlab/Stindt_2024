## The following scripts were used to analyze data for Stindt et. al. 2023, and are organized here based on analysis function (which corresponds roughly to paper section). These utilize data that can be found at the Dryad repository for this paper:  


# These following scripts were directly used to generate all figures in Stindt et. al. 2023:

"All_Plots.m" was used to generate all plotted data (besides modeling), using a master compiled data table "__Consolidated_Data.mat". Sections of the code correspond to Figure numbers.
"Conjugation_TKCsweep.m" was used to generate Figure 4A using "__Conjugation_Clump_Data.mat" and "AllFits.mat".
"Clumping_ParamSweeps.m" was used to generate Figures 4B,C using "__Conjugation_Clump_Data.mat" and "AllFits.mat".
"Rescue_ParamSweeps.m" was used to generate Figure 6C using "__Conjugation_Clump_Data.mat" and "AllFits.mat".



# The following scripts were used to analyze raw data, based on analysis type:

# Flow cytometry data:

"FlowAnalyze.m" imports raw fcs files and sample info, allows drawing/editing gates, calculates grouping statistics, consolidates data from separate voltages, and analyzes samples (e.g. normalization, IDC frequencies, etc.). It requires the following supporting functions:

  - add_gate.m
  - apply_gate.m
  - clean_grpstats.m
  - draw_gate.m
  - fca_readfcs.m
  - format_fcsdat.m
  - load_fcs.m
  - plot_gate.m


# Fluorimetry data:

"TecAnalyzeSparkMulti.m" imports consolidated raw Excel data outputs from Tecan Spark plate reader, organizes it, applies sample info, and analyzes sample data. 


# Colony image ("expansion assay") data:

"ExpAnalyze.m" imports .tif files, calculates colony-wide channel data including average intensities and Li colocalization, and calculates radial intensity info. It requires the following supporting functions:

  - ColID.m: determines colony definitions based on outermost fluorescent circle, flags problematic images
  - LiColoc.m: calculates Li Colocalization for entire colony
  - RadialScan.m: performs intensity scans for each fluorescent channel, from center of colony to each point along circumference, corrects for oversampling near center,     
    calculates metrics per radial distance


# Clumping image data:

"ClumpAnalyze.m" imports .tif files, filters images by event counts, calculates clump sizes and bacterial coincidence, imputes number clumps and clumped cells for use in modeling.


# Modeling:

Modeling scripts were written and run piecewise, in order of complexity, while building out the ODE model and adding parameters. Parameters were determined in the following order: rates and capacities -> Monod terms and amino acid feeding -> competition / ecological niche overlap -> conjugation -> clumping -> rescue. Most of these were used simply as parameter exploration before developing the "Free", "Clumped", and "Rescue" models.

To run "free cell" model:

  - Conjugation_LHS.m: run latin hypercube sampling model with user-defined number of guesses (paper uses 100,000 iteratively)
  - Conjugation_LHS_Verify.m: runs model with 1 guess, constrained to the top parameter set. This is used to quickly generate top fit (or saved fit) plots
  - Conjugation_TKCsweep.m: sweeps IDC rate term with saved parameters
  - Conjugation_ODEs_wDeath.m: child function for model running, contains ODEs

To run "clumped cell" model:

  - Clumping_LHS.m: run latin hypercube sampling model with user-defined number of guesses (paper uses 100,000 iteratively)
  - Clumping_LHS_Verify.m: runs model with 1 guess, constrained to the top parameter set. This is used to quickly generate top fit (or saved fit) plots
  - Clumping_ParamSweeps.m: sweeps IDC rate term and proximity terms with saved parameters
  - Clump_ODEs.m: child function for model running, contains ODEs

To run "rescue" model:

  - Rescue_ParamSweeps.m: sweeps leucine, uracil, and histidine, with clump model parameters modified for rescue assay conditions. Generates Figure 6C
  - Rescue_ODEs.m: child function for rescue model, contains ODEs

Additional model, parameter exploration:

  - fit_growths.m: Gompertz fits for monoculture fluorescence traces, for obtaining initial rate and capacity terms
  - Nutrient_Dynamics_Info.m: Includes sources for amino acid concentrations and secretion values, performs appropriate conversions
  - Monod.m: Calculates initial Monod terms
  - Competition_LHS.m: performs model fitting (free-cell) for initial competition (aka niche overlap) terms
  - Competition_ODEs.m: child function for previous model fit script, contains ODEs
  - Conjugation_DataPrep.m: converts fluorescence data for use in free-cell model fitting
  - Clumping_DataPrep.m: converts fluorescence data for use in clumped-cell model fitting
 
