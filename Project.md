Machine learning using RNAseq data
================
Alsu
February-March 2024

# Data preparation

``` r
library(janitor)
```

    ## 
    ## Attaching package: 'janitor'

    ## The following objects are masked from 'package:stats':
    ## 
    ##     chisq.test, fisher.test

``` r
library(tidyverse)
```

    ## ── Attaching core tidyverse packages ──────────────────────── tidyverse 2.0.0 ──
    ## ✔ dplyr     1.1.4     ✔ readr     2.1.5
    ## ✔ forcats   1.0.0     ✔ stringr   1.5.1
    ## ✔ ggplot2   3.5.0     ✔ tibble    3.2.1
    ## ✔ lubridate 1.9.3     ✔ tidyr     1.3.1
    ## ✔ purrr     1.0.2

    ## ── Conflicts ────────────────────────────────────────── tidyverse_conflicts() ──
    ## ✖ dplyr::filter() masks stats::filter()
    ## ✖ dplyr::lag()    masks stats::lag()
    ## ℹ Use the conflicted package (<http://conflicted.r-lib.org/>) to force all conflicts to become errors

``` r
setwd("/cloud/project/data/")  

normalized.counts.ibd <- read.table(file="Normalized_Count.csv",
                             sep="",
                             header=T,
                             fill=T,
                             check.names=F)

# Transpose data frame
data <- t(normalized.counts.ibd)

# Move the gene names as headers, add column for status
data <- data %>%
  row_to_names(row_number = 1) %>%
  as.data.frame() # %>%
```

``` r
#install.packages('rsample')

library(rsample)

# Fix random numbers by setting the seed 
# Enables analysis to be reproducible when random numbers are used 
set.seed(123)
# Put 80% of the data into the training set 
data_split <- initial_split(data, prop = 0.70)
# Create data frames for the two sets:
train_data <- training(data_split)
test_data  <- testing(data_split)
```

# Logistic regression after PCA

### Dimensionality reduction (PCA)

``` r
#convert character dataframe to numeric
train_data <- as.data.frame(sapply( train_data , as.numeric))

# Identify constant or zero columns
train_data_transposed_cons <- sapply(train_data, function(x) is.atomic(x) && length(unique(x)) == 1)

# Remove constant or zero columns
train_data_transposed_no <- train_data[, !train_data_transposed_cons]

# Run PCA
pca_train_result <- prcomp(train_data_transposed_no, scale. = TRUE)
```

``` r
# Plot PCA results
plot(pca_train_result$x[,1], pca_train_result$x[,2], 
     xlab = "PC1", ylab = "PC2", 
     main = "PCA of Normalized Counts for the training data")
```

![](Project_files/figure-gfm/visualize-pca-1.png)<!-- -->

``` r
# 3D scatterplot

#install.packages("scatterplot3d") #install as needed
library(scatterplot3d)

scatterplot3d(pca_train_result$x[,1], pca_train_result$x[,2], pca_train_result$x[,3], 
              xlab = "PC1", ylab = "PC2", zlab = "PC3", 
              main = "3D PCA Plot")
```

![](Project_files/figure-gfm/visualize-pca-2.png)<!-- -->

``` r
#Screeplot
screeplot(pca_train_result, type = "lines")
```

![](Project_files/figure-gfm/visualize-pca-3.png)<!-- -->

``` r
#heatmap
train_matrix <- as.matrix(train_data[, -1])   
train_matrix <- matrix(as.numeric(unlist(train_matrix)),nrow=nrow(train_matrix))

#heatmap(train_matrix, scale = "row", add.expr = TRUE) #PositCloud does not have enought computing power for this one
```

### Logistic regression

``` r
# Select first few principal components (based on screeplot)
num_components <- 6
selected_pcs <- pca_train_result$x[, 1:num_components]

#Add binary variable (parkinson's vs control). Condition is written for the seed(123)
condition <- c(rep("Control", 1), rep("Parkinsons", 2), rep("Control", 4), rep("Parkinsons", 2))

# Assuming 'condition' is the binary outcome variable
data_train_selected_pcs <- cbind(selected_pcs, condition)
```

``` r
library(broom)
library(parsnip)

data_train_selected_pcs <- as.data.frame(data_train_selected_pcs)

#Fitting the logistic regression model
model_fit <- logistic_reg() %>%
  set_engine("glm") %>%
  fit(as.factor(condition) ~ ., data = data_train_selected_pcs, family = "binomial")

tidy(model_fit)
```

    ## # A tibble: 49 × 5
    ##    term                  estimate std.error statistic p.value
    ##    <chr>                    <dbl>     <dbl>     <dbl>   <dbl>
    ##  1 (Intercept)          -2.46e+ 1   131011. -1.88e- 4    1.00
    ##  2 PC1-100.669753612644  4.91e+ 1   185277.  2.65e- 4    1.00
    ##  3 PC1-2.37775419621684 -2.69e-14   185277. -1.45e-19    1   
    ##  4 PC1-39.1392375268399 -4.99e-14   185277. -2.69e-19    1   
    ##  5 PC1-9.62955839720086 -5.83e-15   185277. -3.15e-20    1   
    ##  6 PC119.0287120199584   4.91e+ 1   185277.  2.65e- 4    1.00
    ##  7 PC12.73379426679209  -1.44e-14   185277. -7.75e-20    1   
    ##  8 PC132.5721926050154   4.91e+ 1   185277.  2.65e- 4    1.00
    ##  9 PC198.4609141258824   4.91e+ 1   185277.  2.65e- 4    1.00
    ## 10 PC2-27.1763035517658 NA              NA  NA          NA   
    ## # ℹ 39 more rows

``` r
#Assessing the goodness of fit of the logistic regression model
summary(model_fit)
```

    ##              Length Class        Mode     
    ## lvl           2     -none-       character
    ## spec          8     logistic_reg list     
    ## fit          30     glm          list     
    ## preproc       1     -none-       list     
    ## elapsed       1     -none-       list     
    ## censor_probs  0     -none-       list

``` r
# Get coefficient estimates and their significance
#coef(summary(model_fit))
```

# Logistic regression without PCA

``` r
#install.packages("glmnet", "parsnip") #install as needed

library(glmnet)
```

    ## Loading required package: Matrix

    ## 
    ## Attaching package: 'Matrix'

    ## The following objects are masked from 'package:tidyr':
    ## 
    ##     expand, pack, unpack

    ## Loaded glmnet 4.1-8

``` r
library(parsnip)

lr_mod <- 
  logistic_reg(penalty = tune(), mixture = 1) %>% 
  set_engine("glmnet")
```
