rm(list = ls())
library(MetaboAnalystR)
library(OptiLCMS)



## MS1 DATA
# 1.Extract an Region of Interest (ROI) for parameters optimization.
unzip("malaria_raw.zip", exdir = "upload")
file_list <- list.files("upload/malaria_LCMS/", pattern = ".zip", full.names = TRUE)
for (file in file_list){
  unzip(file, exdir = "upload/")
}
DataFiles <- list.files("upload/", pattern = "QC", full.names = TRUE)

mSet <- PerformROIExtraction(datapath = DataFiles, rt.idx = 0.9, rmConts = TRUE)
# rt.idx refers to the retained percentage of chromatogram in retention time dimension. Suggested to use more than 50% (0.5) to include larger ROI and improve the accuracy of following parameters’ optimization.
# rmConts means “remove contamination” before ROI extraction.

# 2.Auto-optimization of parameters
best_params <- PerformParamsOptimization(mSet, param = SetPeakParam(platform = "UPLC-Q/E"), ncore = 1) #windows only recommend 1 core.

# 3.Importing example data (both metadata file and sub folders can be used)
# Plotting functions can be enabled with parameter
mSet <- ImportRawMSData(path = c("upload"), metadata = "metadata.txt", plotSettings = SetPlotParam(Plot = T))

# 4.Raw spectral data processing
mSet <- PerformPeakProfiling(mSet, Params = best_params, plotSettings = SetPlotParam(Plot=TRUE))

# 5. Feature annotation.  (identify real peaks, and to clarify the relationships among them, using CAMERA40 and CliqueMS41)
# 'polarity' is required, can be either 'negative' or 'positive';
# 'perc_fwhm' is used to set the percentage of the width of the FWHM for peak grouping. 
#              Default is set to 0.6;
# 'mz_abs_iso' is used to set the allowed variance for the search (for isotope annotation). 
#              The default is set to 0.005;
# 'max_charge' is set the maximum number of the isotope charge. 
#              For example, the default is 2, therefore the max isotope charge is 2+/-;
# 'max_iso' is used to set the maximum number of isotope peaks.
#              For example, the default is 2, therefore the max number of isotopes per peaks is 2;
# 'corr_eic_th' is used to set the threshold for intensity correlations across samples. 
#              Default is set to 0.85.
# 'mz_abs_add' is used to set the allowed variance for the search (for adduct annotation). 
#              Default is set to 0.001.
# 'adducts' is used to specify the adducts based on your instrument settings.
annParams <- SetAnnotationParam(polarity = 'positive', mz_abs_add = 0.015)

# "annParams" is the parameters used for annotation
mSet <- PerformPeakAnnotation(mSet, annParams)

# 6.Feature table generation
# Here we format and filter the peak list for following analysis with MetaboAnalystR

# annParams, is the object created using the SetAnnotationParam function above;
# filtIso, is used to decide to filter out all isotopes (TRUE) or not (FALSE);
# filtAdducts, is used to decide to filter out all adducts (TRUE) or not (FALSE);
# missPercent, specify the threshold to remove features missing in a certain percentage of all samples in a group.

mSet <- FormatPeakList(mSet, annParams, filtIso = FALSE, filtAdducts = FALSE, missPercent = 1)

# Export annotation results, the annotation will be save as "annotated_peaklist.csv";
Export.Annotation(mSet)

# Export complete feature table. It will be saved as "metaboanalyst_input.csv";
Export.PeakTable(mSet)

# Export a summary table (peak_result_summary.txt) to summarize the information of all peaks in
# different samples. The columns are sample names, groups, retention time range, m/z range of all peaks,
# number of all peaks and percentage of missing features.
Export.PeakSummary(mSet)



## DDA MS/MS data
# 1.prepare raw data
unzip("ms2_dda_example.zip", exdir = "ms2_dda")

# 2.Download MS/MS spectra reference database (https://www.metaboanalyst.ca/resources/vignettes/LCMSMS_Raw_Spectral_Processing.html)
download.file("https://rest.xialab.ca/api/download/metaboanalyst/MS2ID_Bio_v09102023.zip",
              destfile = "MS2ID_Bio.zip",
              method = "curl")
unzip("MS2ID_Bio.zip", exdir = "ms2_db")

ft_dt <- qs::qread("ms2_dda/ms2_dda_example/ft_dt.qs")
# See the previous section for automatic optimization and the following section for your own settings

mSet <- PerformMSnImport(filesPath = c(list.files("ms2_dda/ms2_dda_example/",
                                                  pattern = ".mzML",
                                                  full.names = T, recursive = T)),
                         targetFeatures = ft_dt,
                         acquisitionMode = "DDA")

# 3.DDA MS/MS spectra Deconvolution
system.time(mSet <- PerformDDADeconvolution(mSet,
                                            ppm1 = 5,
                                            ppm2 = 10,
                                            sn = 12,
                                            filtering = 0,
                                            window_size = 1.5,
                                            intensity_thresh = 1.6e5,
                                            database_path = "ms2_db/MS2ID_Bio_v09102023.sqlite",
                                            ncores = 2L))
# 4.Spectrum consensus of replicates
mSet <- PerformSpectrumConsenus (mSet,
                                 ppm2 = 15,
                                 concensus_fraction = 0.2,
                                 database_path = "",
                                 use_rt = FALSE,
                                 user_dbCorrection = FALSE)
# 5.Reference database searching
mSet <- PerformDBSearchingBatch (mSet,
                                 ppm1 = 10,
                                 ppm2 = 25,
                                 rt_tol = 5,
                                 database_path = "ms2_db/MS2ID_Bio_v09102023.sqlite",
                                 use_rt = FALSE,
                                 enableNL = FALSE,
                                 ncores = 2L)
# 6.Results export
mSet <- PerformResultsExport (mSet,
                              type = 0L,
                              topN = 10L,
                              ncores = 2L)
# 2nd argument, TopN, is the argument used to specifiy the number of compounds to export

dtx <- FormatMSnAnnotation(mSet, 5L, F)
# If the dataset is a lipidomics dataset, please set 3rd argument as "TRUE" to extract the lipid classification information

# mzmin, minimum m/z value for this feature;
# mzmax, maximum m/z value for this feature;
# rtmin, minimum retention time value for this feature;
# rtmax, maximum retention time value for this feature;
# Compound_N, Compound name of the Nth identification. The identification results are sorted based on the scores (decreasing);
# InChiKey_N, InChiKey of the Nth identification;
# Formula_N, Formula of the Nth identification;
# Score_N, Similarity score of the Nth identification. 100 is the perfect match. 0 is the negative match;
# Database_N, the MS/MS reference library source of the Nth identification;
# SuperClass_N, Superclass of the Nth identification (only for lipidomics database).
# MainClass_N, Mainclass of the Nth identification (only for lipidomics database).
# SubClass_N, Subclass of the Nth identification (only for lipidomics database).

# 7.Matching results exploration
# Here we are downloading biology MS/MS reference database.
download.file("https://rest.xialab.ca/api/download/metaboanalyst/FragsAnnotateDB_v02042024.zip",
              destfile = "FragsAnnotateDB.zip",
              method = "curl")
unzip("FragsAnnotateDB.zip", exdir = "ms2_db")

save(mSet, file = "mSet_raw_dda.rda") # you are suggested to save this, in case of lossing data

# Here, user need to convert raw spec mSet object into regular analysis object;
# 1. peak_idx: Index of the peak. For example, 11 refers to the 11th peak in the target peak list;
# 2. sub_idx: Index of the identification result. For example, 1 refers to the first identified compound (the one with the highest matching score.)
# 3. interative: can be TRUE or FALSE. If TRUE, will plot an interactive figure.
# 4. ppm, the parameter used to determine the mactching results.
mSet <- Convert2AnalObject(mSet, "raw", "spec", F)

mSet <- PerformMirrorPlotting(mSetObj = mSet, 
                              fragDB_path = "ms2_db/FragsAnnotateDB_v02042024.sqlite", 
                              peak_idx = 11, sub_idx = 1, 
                              interactive = T, ppm = 25, 
                              dpi = 72, format = "png", width = 8, height = 6)



## SWATH-DIA MS/MS data
# 1.data preparing
download.file("https://rest.xialab.ca/api/download/metaboanalyst/ms2_dia_example.zip",
              destfile = "ms2_dia_example.zip",
              method = "curl")

unzip("ms2_dia_example.zip", exdir = "ms2_dia")

# Construct meta data
meta_dt <- data.frame(samples = c("210210-SWATH-NEG-Covid-Cov-16.mzML", 
                                  "210210-SWATH-NEG-Covid-Cov-17.mzML", 
                                  "210210-SWATH-NEG-Covid-Cov-18.mzML",
                                  "210210-SWATH-NEG-Covid-Ct-1.mzML",
                                  "210210-SWATH-NEG-Covid-Ct-2.mzML",
                                  "210210-SWATH-NEG-Covid-Ct-3.mzML"),
                      groups = c(rep("COVID",3), rep("Control",3)))

# Import raw MS1 data
mSet1 <- ImportRawMSData(path = "ms2_dia/ms2_dia_example/", metadata = meta_dt)

# Process MS1 peaks by using customized parameters
mSet1 <- PerformPeakProfiling(mSet1, 
                              Params = SetPeakParam(ppm = 25,
                                                    bw = 3,
                                                    mzdiff = 0.001,
                                                    max_peakwidth = 35,
                                                    min_peakwidth = 5,
                                                    noise = 200, 
                                                    minFraction = 0.2),
                              ncore = 6,
                              plotSettings = SetPlotParam(Plot = F))

# Annotation
annParams <- SetAnnotationParam(polarity = 'negative',
                                mz_abs_add = 0.035);

mSet1 <- PerformPeakAnnotation(mSet = mSet1,
                               annotaParam = annParams,
                               ncore =1)

# Format MS1 feature table
mSet1 <- FormatPeakList(mSet = mSet1,
                        annParams,
                        filtIso =FALSE,
                        filtAdducts = FALSE,
                        missPercent = 1)
# export MS1
# Export.Annotation(mSet1)
# Export.PeakTable(mSet1)
# Export.PeakSummary(mSet1)

mSet <- PerformMSnImport(mSet = mSet1,
                         filesPath = c(list.files("ms2_dia/ms2_dia_example/", pattern = ".mzML", full.names = T, recursive = T)),
                         acquisitionMode = "DIA",
                         SWATH_file = "ms2_dia/ms2_dia_example/DIA_SWATH_MS_experiment_file_neg.txt")

# Deconvolution
mSet <- PerformDIADeconvolution(mSet,
                                min_width = 5,
                                span = 0.3,
                                ppm2 = 30,
                                sn = 12,
                                filtering = 0,
                                ncores = 6L)

# Consenus
mSet <- PerformSpectrumConsenus (mSet,
                                 ppm2 = 30,
                                 concensus_fraction = 0.25,
                                 database_path = "",
                                 use_rt = FALSE,
                                 user_dbCorrection = FALSE)

# Reference database searching
mSet <- PerformDBSearchingBatch (mSet,
                                 ppm1 = 15,
                                 ppm2 = 30,
                                 rt_tol = 5,
                                 database_path = "ms2_db/MS2ID_Bio_v09102023.sqlite",
                                 use_rt = FALSE,
                                 enableNL = FALSE,
                                 ncores = 4L)
# Export
mSet <- PerformResultsExport (mSet,
                              type = 0L,
                              topN = 10L,
                              ncores = 4L)

dtx2 <- FormatMSnAnnotation(mSet, 5L, F)

# convert raw spec mSet object into regular analysis object
save(mSet, file = "mSet_raw_dia.rda")
mSet <- Convert2AnalObject(mSet, "raw", "spec", F)
mSet <- PerformMirrorPlotting(mSetObj = mSet, 
                              fragDB_path = "ms2_db/FragsAnnotateDB_v02042024.sqlite", 
                              peak_idx = 5, sub_idx = 1, 
                              interactive = T, ppm = 25, 
                              dpi = 72, format = "png", width = 8, height = 6)