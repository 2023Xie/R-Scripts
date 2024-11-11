# Load the R package
library(readxl)
library(caret)
library(randomForest)
library(ggplot2)


# Get all sheet names
sheets <- excel_sheets("Tables.S2.xlsx")

# Create data storage list
data <- list()

# Loop reading worksheet data
for (shee in sheets) {
        # Read worksheet data
        da <- read_excel("Tables.S2.xlsx", sheet = shee, skip = 1, col_names = TRUE)
        # Convert group column to factor type
        da$groups <- as.factor(da$groups)
        # Store data in a list and name it after the worksheet name
        data[[shee]] <- da
        # Delete da
        rm(da)
}
# Delete NA value (there is no NA value in this experimental data)
data1<-na.omit(data)

# Setting the random number seed ensures that the result is repeatable
set.seed(430)

# Using serum data as an example for random forest analysis
# Create training and testing sets using the "caret" package
partition <- createDataPartition(data1[["serum"]]$groups, p = 0.7, list = FALSE)
train_set <- data1[["serum"]][partition, ]
test_set <- data1[["serum"]][-partition, ]

# Use the "randomForest" package and test set for testing analysis
test.rf = randomForest(groups ~ ., data=test_set, ntree=200, importance=TRUE, proximity=TRUE)
plot(test.rf)

# Select the number of decision trees based on test results and use the training set to create a random forest model
randomForest.rf = randomForest(groups ~ ., data=train_set, ntree=100, importance=TRUE, proximity=TRUE)
plot(randomForest.rf)

# Use the test set to evaluate predictive performance
test_pred<-predict(randomForest.rf, newdata=test_set)
table(test_pred, test_set$groups)

# Using the training set to evaluate prediction accuracy
train_pred<-predict(randomForest.rf, newdata=train_set)
table(train_pred, train_set$groups)

# Evaluate the importance of predictive variables
importance <- randomForest.rf$importance
head(importance)

# Visualize the "importance" data
varImpPlot(randomForest.rf, n.var = min(30, nrow(randomForest.rf$importance)),
           main = "Variable importance")

# Ten-fold cross validation
train.cv <- replicate(5, rfcv(train_set[-ncol(train_set)], train_set$groups, cv.fold = 10, step = 1.5), simplify = FALSE)
train.cv
# Extract the validation results
train.cv <- data.frame(sapply(train.cv, '[[', 'error.cv'))
train.cv$Metabolites <- rownames(otu_train.cv)
train.cv <- reshape2::melt(train.cv, id = 'Metabolites')
train.cv$Metabolites <- as.numeric(as.character(train.cv$Metabolites))
train.cv.mean <- aggregate(train.cv$value, by = list(train.cv$Metabolites), FUN = mean)
head(train.cv.mean, 10)

# Plot
ggplot(train.cv.mean, aes(Group.1, x)) +
          geom_line() +
          theme(panel.grid = element_blank(), panel.background = element_rect(color = 'black', fill = 'transparent')) +
          labs(title = '',x = 'Metabolites', y = 'Cross-validation serror')
