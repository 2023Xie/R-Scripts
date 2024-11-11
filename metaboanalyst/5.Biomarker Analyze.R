library(MetaboAnalystR)
rm(list = ls())

# Create objects for storing processed data from biomarker analysis
mSet<-InitDataObjects("conc", "roc", FALSE)
# Read in data and fill in the dataSet list
mSet<-Read.TextData(mSet, "http://www.metaboanalyst.ca/MetaboAnalyst/resources/data/plasma_nmr_new.csv", "rowu", "disc")
# Sanity check, replace missing values, check if the sample size is too small
mSet<-SanityCheckData(mSet)
mSet<-ReplaceMin(mSet)
mSet<-IsSmallSmplSize(mSet)
mSet<-PreparePrenormData(mSet)

### *** OPTION 1 FOR NORMALIZATION
# Perform no normalization, no ratio calculation
mSet<-Normalization(mSet, "NULL", "NULL", "NULL", ref=NULL, ratio=FALSE, ratioNum=20)
### To perform normalization with ratio calculation use Option 2. . .
### *** OPTION 2 FOR NORMALIZATION
# No normalization, and computeS metabolite ratios and includeS the top 20
mSet<-Normalization(mSet, "NULL", "NULL", "NULL", "C01", ratio=TRUE, ratioNum=20)
# If ratio = TRUE: view the normalized dataset including the top ranked ratios
# The ratios will be towards the end of the matrix
mSet$dataSet$norm
#If ratio = TRUE: view just the top ranked included ratios
mSet$dataSet$ratio


# Set the biomarker analysis mode to perform Classical ROC curve analysis ("univ")
mSet<-SetAnalysisMode(mSet, "univ")
# Prepare data for biomarker analysis
mSet<-PrepareROCData(mSet)
### OPTION 1 Perform univariate ROC curve analysis ###
mSet<-Perform.UnivROC(mSet, feat.nm = "Isoleucine", imgName = "Isoleucine", "png", dpi=300, isAUC=F, isOpt=T, optMethod="closest.topleft", isPartial=F, measure="sp", cutoff=0.2)
# Create box plot showing the concentrations of the selected compound between the groups
mSet<-PlotBoxPlot(mSet, "Isoleucine", "Isoleucineboxplot_0_", "png", 72, T, FALSE)
# Perform calculation of feature importance (AUC, p value, fold change)
mSet<-CalculateFeatureRanking(mSet)
### OPTION 2 Perform univariate ROC curve analysis, resulting in a partial AUC with a 95% CI band ###
mSet<-Perform.UnivROC(mSet, feat.nm = "Valine", imgName = "Valine", "png", dpi=300, isAUC=T, isOpt=T, optMethod="closest.topleft", isPartial=T, measure="se", cutoff=0.2)
### OPTION 3 Perform univariate ROC curve analysis on a metabolite ratio pair, note that you cannot save an image with a "\" in the name ###
mSet<-Perform.UnivROC(mSet, feat.nm = "Isoleucine/Valine", imgName = "IsoleucineValine", "png", dpi=300, isAUC=T, isOpt=T, optMethod="closest.topleft", isPartial=T, measure="se", cutoff=0.2)



# Set the biomarker analysis mode to perform Multivariate exploratory ROC curve analysis ("explore")
mSet<-SetAnalysisMode(mSet, "explore")
# Prepare data for biomarker analysis
mSet<-PrepareROCData(mSet)
# Perform multivariate ROC curve analysis, using SVM classification and ranking
mSet<-PerformCV.explore(mSet, cls.method = "svm", rank.method = "svm", lvNum = 2)
### OPTION 1 Comparison plot of ROC curves of all models ###
mSet<-PlotROC(mSet, imgName = "ROC_all_models", format = "png", dpi = 300, mdl.inx= 0, avg.method = "threshold", show.conf = 0, show.holdout = 0, focus="fpr", cutoff=0.5)
# Plot predicted class probabilities for each sample for a selected model, not showing labels of wrongly classified samples
mSet<-PlotProbView(mSet, imgName = "multi_roc_prob", format = "png", dpi = 300, mdl.inx = -1, show = 0, showPred = 0)
# Plot the predictive accuracy of models with increasing number of features
mSet<-PlotAccuracy(mSet, imgName = "multi_roc_accuracy", format = "png", dpi = 300)
# Plot the most important features of a selected model ranked from most to least important
mSet<-PlotImpVars(mSet, imgName = "multi_roc_impvar", format="png", dpi=300, mdl.inx = -1, measure="freq", feat.num=15)



### The above workflow plotted ROC curves for all models, to plot the ROC curve for a single model instead use:
### OPTION 2 Plot the ROC curve of a single selected model, in this case model 1 and display the confidence interval ###
mSet<-PlotROC(mSet, imgName = "ROC_model1", format = "png", dpi = 300, mdl.inx = 1, avg.method = "threshold", show.conf = 1, 0, "fpr", 0.2)
# To view the top-ranked 10 variables use the code below.
imp.feats <- GetImpFeatureMat(mSet, mSet$analSet$multiROC$imp.cv, 10)
head(imp.feats)