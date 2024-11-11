library(MetaboAnalystR)
rm(list = ls())

## When input is a list
# Create vector consisting of compounds for enrichment analysis 
tmp.vec <- c("Acetoacetic acid", "Beta-Alanine", "Creatine", "Dimethylglycine", "Fumaric acid", "Glycine", "Homocysteine", "L-Cysteine", "L-Isolucine", "L-Phenylalanine", "L-Serine", "L-Threonine", "L-Tyrosine", "L-Valine", "Phenylpyruvic acid", "Propionic acid", "Pyruvic acid", "Sarcosine", "Arsenic", "Benzene", "Caffeic acid", "Cotinine", "Cadmium", "Lead", "Thiocyanate")
# Create mSetObj
mSet<-InitDataObjects("conc", "msetora", FALSE)
#Set up mSetObj with the list of compounds
mSet<-Setup.MapData(mSet, tmp.vec)
# Cross reference list of compounds against libraries (hmdb, pubchem, chebi, kegg, metlin)
mSet<-CrossReferencing(mSet, "name")
# Example compound name map
mSet$name.map 
# Create the mapping results table
mSet<-CreateMappingResultTable(mSet)
# Input the name of the compound without any matches 
mSet<-PerformDetailMatch(mSet, "L-Isolucine")
# Create list of candidates to replace the compound
mSet <- GetCandidateList(mSet)
# Identify the name of the compound to replace
mSet<-SetCandidate(mSet, "L-Isolucine", "L-Isoleucine")
# Set the metabolite filter
mSet<-SetMetabolomeFilter(mSet, F)
# Select metabolite set library, refer to 
mSet<-SetCurrentMsetLib(mSet, "smpdb_pathway", 0)
# Calculate hypergeometric score, results table generated in your working directory
mSet<-CalculateHyperScore(mSet)
# Plot the ORA, bar-graph
mSet<-PlotORA(mSet, "ora_0_", "bar", "png", 72, width=NA)



# Create mSetObj
mSet<-InitDataObjects("conc", "msetqea", FALSE)
# Read in data table
mSet<-Read.TextData(mSet, "https://rest.xialab.ca/api/download/metaboanalyst/human_cachexia.csv", "rowu", "disc")
# Perform cross-referencing of compound names
mSet<-CrossReferencing(mSet, "name")
# Create mapping results table
mSet<-CreateMappingResultTable(mSet)
# Mandatory check of data 
mSet<-SanityCheckData(mSet)
# Replace missing values with minimum concentration levels
mSet<-ReplaceMin(mSet)
# Perform no normalization
mSet<-PreparePrenormData(mSet)
mSet<-Normalization(mSet, "NULL", "NULL", "NULL", "PIF_178", ratio=FALSE, ratioNum=20)
# Plot normalization
mSet<-PlotNormSummary(mSet, "norm_0_", "png", 72, width=NA)
# Plot sample-wise normalization
mSet<-PlotSampleNormSummary(mSet, "snorm_0_", "png", 72, width=NA)
# Set the metabolome filter
mSet<-SetMetabolomeFilter(mSet, F)
# Set the metabolite set library to pathway
mSet<-SetCurrentMsetLib(mSet, "smpdb_pathway", 0)
# Calculate the global test score
mSet<-CalculateGlobalTestScore(mSet)
# Plot the QEA
mSet<-PlotQEA.Overview(mSet, "qea_0_", "bar", "png", 72, width=NA)



## path way
# Create vector consisting of compounds for enrichment analysis 
tmp.vec <- c("Acetoacetic acid", "Beta-Alanine", "Creatine", "Dimethylglycine", "Fumaric acid", "Glycine", "Homocysteine", "L-Cysteine", "L-Isolucine", "L-Phenylalanine", "L-Serine", "L-Threonine", "L-Tyrosine", "L-Valine", "Phenylpyruvic acid", "Propionic acid", "Pyruvic acid", "Sarcosine")
# Create mSetObj for storing objects created during your analysis
mSet<-InitDataObjects("conc", "pathora", FALSE)
# Set up mSetObj with the list of compounds
mSet<-Setup.MapData(mSet, tmp.vec)
# Cross reference list of compounds against libraries (hmdb, pubchem, chebi, kegg, metlin)
mSet<-CrossReferencing(mSet, "name")
mSet<-CreateMappingResultTable(mSet)
# Perform matching against compound with out matches
mSet<-PerformDetailMatch(mSet, "L-Isolucine")
# Get list of candidates for matching
mSet<-GetCandidateList(mSet)
# Replace selected compound
mSet<-SetCandidate(mSet, "L-Isolucine", "Isoleucine")
# Select the pathway library, ranging from mammals to prokaryotes
mSet<-SetKEGG.PathLib(mSet, "hsa", "current")
# Set the metabolite filter
mSet<-SetMetabolomeFilter(mSet, F)
# Calculate the over representation analysis score, here we selected to use the hypergeometric test (alternative is Fisher's exact test)
# A results table "pathway_results.csv" will be created and found within your working directory
mSet<-CalculateOraScore(mSet, "rbc", "hyperg")
# Plot of the Pathway Analysis Overview 
mSet<-PlotPathSummary(mSet, T, "path_view_0_", "png", 72, width=NA)
# Plot a specific metabolic pathway, in this case "Glycine, serine and threonine metabolism"
mSet<-PlotKEGGPath(mSet, "Glycine, serine and threonine metabolism",576, 480, "png", NULL)
