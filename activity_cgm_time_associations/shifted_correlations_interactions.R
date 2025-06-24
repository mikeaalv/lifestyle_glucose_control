# load data and the needed libraries 
rm(list = ls())
library(plyr)
library(dplyr)
library(tidyverse)
library(ggplot2)
library(reshape)
library(reshape2); 
library(grid)
library(gridExtra)
library(scales)
library(RColorBrewer)
library(extrafont)
library(ggpubr)
library(factoextra)
library(matrixStats)
library("writexl")
library(psycho)
library(lubridate)
library("rjson")
library(readtext)
library(readr)
library(stringr)
library(viridis)
library(hms)
library(anytime)
loadfonts()
detach("package:plyr", unload = TRUE)

# load the formatted daily zones data (for activity and CGM data)
daily_zones_data = read_csv("daily_zones_cgm_physical_activity_diet.csv")
daily_zones_data <- daily_zones_data %>% select(-X1)

#################################################################################################
#correlation of steps with cgm in different time zones with sampling and permutation
#################################################################################################
#for each metabolic subphenotype
#sample sample_size consecutive days for each participant
#combine data for all participants of that metabolic subphenotype
#perform a shifted pearson correlation analysis between cgm and steps
#repeat perm_num permutations of the shifted correlation analysis
#calculate average correlation coefficients from perm_num permutations
#calculate p-value of the perm_num permutations (test if correlation coefficients from permutations are significantly different from 0 as if there was no correlation)
####################################################################################


###################################################
###################################################
#steps v. cgm next 0-24 h
###################################################
###################################################

# initialize parameters
set.seed(5) 
sample_size = 25
perm_num = 8 # previously 10
perm_num_minus = perm_num - 1
i_sub = 0
plots<- list()

# data.sub: subset of data with the u_sub subgroup status 
data.sub <- daily_zones_data

if (nrow(data.sub) <= sample_size){next} # insufficient data, given sample size

# initialize the results to be stored across the permutations
test <- list()
test_interaction <- list()
test_SDinteraction <- list()
test_pvalinteraction <- list()
all_predictions <- list()

df <- data.frame(matrix("", ncol = 7, nrow = 7))  
x <- c("steps.zone.corr1", "steps.zone.corr2", "steps.zone.corr3", "steps.zone.corr4", "steps.zone.corr5", "steps.zone.corr6", "steps.zone.corr7")
colnames(df) <- x

df <- as.data.frame(sapply(df, as.numeric))

x <- c("CGM.zone1", "CGM.zone2", "CGM.zone3", "CGM.zone4", "CGM.zone5", "CGM.zone6", "CGM.zone7")
rownames(df) <- x

# iterate through the permutations 
for (t in 1:perm_num)
{
  
  sampling <- data.sub[1,]
  sampling <- sampling[-1,]
  
  for (u in (unique(data.sub$ID))) # iterate through the participants 
    
  {
    data.temp <- data.sub %>% subset(ID == u) # data for the individual
    if (nrow(data.temp) < sample_size){next} # insufficient data, given sample size
    # extract sample of the data 
    sampling_index = nrow(data.temp) - sample_size + 1
    start_index <- sample(1:sampling_index, 1)
    end_index <- start_index + sample_size - 1
    sampling <- rbind(sampling,data.temp[start_index:end_index, ]) # sample of data for individual 
  }
  
  # initialize for the permutation
  test[[t]] <- df
  test_interaction[[t]] <- df
  test_SDinteraction[[t]] <- df
  test_pvalinteraction[[t]] <- df
  
  num_cgm_start = 6 # based on the data structure 
  num_steps_start = 13
  
  # initialize an empty dataframe for storing all predictions
  all_predictions[[t]] <- data.frame()
  j_labels <- c("steps: 12am-5am", "steps: 5am-8am", "steps: 8am-11am", "steps: 11am-2pm", "steps: 2pm-5pm", "steps: 5pm-9pm", "steps: 9pm-12am")
  i_labels <- c("12am-5am", "5am-8am", "8am-11am", "11am-2pm", "2pm-5pm", "5pm-9pm", "9pm-12am")
  
  for (i in 1:7)
  {
    #i for CGM
    for (j in 1:7)
    {
      #j for steps
      temp <- sampling
      #skip if more than half of the cgm or carb values are zero for that time-window
      if (sum(temp[,i+num_cgm_start] == 0) > nrow(temp)/4*3){next}
      if (sum(temp[,j+num_steps_start] == 0) > nrow(temp)/4*3){next}
      cgm.temp <- temp[,i+num_cgm_start]
      cgm.temp <- cgm.temp[!is.na(cgm.temp)]
      steps.temp <- temp[,j+num_steps_start]
      steps.temp <- steps.temp[!is.na(steps.temp)]
      res <- cor.test(cgm.temp, steps.temp, method = "pearson")
      test[[t]][i,j] <- res$estimate
      
      # set up the linear model
      numeric_indices <- c(i+num_cgm_start, j+num_steps_start)
      selected_columns <- temp[, numeric_indices]
      selected_columns <- cbind(selected_columns, temp[,"sspg_status"])
      names(selected_columns)[1] <- "CGM"
      names(selected_columns)[2] <- "steps"
      
      # fit the model, with steps*sspg interaction
      model <- lm(CGM ~ steps + sspg_status + steps:sspg_status, data = selected_columns)
      model_summary <- summary(model)
      coefficients <- as.data.frame(model_summary$coefficients)
      test_interaction[[t]][i,j] <- coefficients[4,1] # estimate for interaction term
      test_SDinteraction[[t]][i,j] <- coefficients[4,2] # std. error for interaction
      test_pvalinteraction[[t]][i,j] <- coefficients[4,4] # p-value for interaction
      
      # create a new data frame for prediction
      steps_range <- seq(0, 0.8, length.out = 100)
      sspg_status_levels <- na.omit(unique(selected_columns$sspg_status))
      prediction_data <- expand.grid(steps = steps_range, sspg_status = sspg_status_levels)
      prediction_data$CGM_predicted <- predict(model, newdata = prediction_data)
      
      prediction_data$i_index <- i
      prediction_data$j_index <- j
      
      # add custom labels for the current iteration
      prediction_data$i_label <- i_labels[i]
      prediction_data$j_label <- j_labels[j]
      
      # combine with the all_predictions dataframe
      all_predictions[[t]] <- rbind(all_predictions[[t]], prediction_data)
      
    }}
}

combined_df <- do.call(rbind, all_predictions)

# now calculate the mean of the predictions for each combination
average_predictions_0_24 <- combined_df %>%
  group_by(steps, sspg_status, i_index, j_index, i_label, j_label) %>%
  summarize(CGM_predicted = mean(CGM_predicted, na.rm = TRUE))

j_labels_reversed = rev(j_labels)
average_predictions_0_24$i_label <- factor(average_predictions_0_24$i_label, levels = i_labels)
average_predictions_0_24$j_label <- factor(average_predictions_0_24$j_label, levels = j_labels_reversed)

# combine correlation matrices from permutations
test.combined <- test[[1]][1,]
test.combined <- test.combined[-1,]

test_interaction.combined <- test_interaction[[1]][1,]
test_interaction.combined <- test_interaction.combined[-1,]

test_SDinteraction.combined <- test_SDinteraction[[1]][1,]
test_SDinteraction.combined <- test_SDinteraction.combined[-1,]

test_pvalinteraction.combined <- test_pvalinteraction[[1]][1,]
test_pvalinteraction.combined <- test_pvalinteraction.combined[-1,]

# combine the results across permutations 
for (i in 1:perm_num)
{
  test.combined <- rbind(test.combined,test[[i]])
  test_interaction.combined <- rbind(test_interaction.combined,test_interaction[[i]])
  test_SDinteraction.combined <- rbind(test_SDinteraction.combined,test_SDinteraction[[i]])
  test_pvalinteraction.combined <- rbind(test_pvalinteraction.combined,test_pvalinteraction[[i]])
}


# initialize avg correlation matrices from permutations
initialize_matrix <- function(base_name) {
  mat <- matrix(NA, nrow = 7, ncol = 7)
  rownames(mat) <- paste0("CGM.zone", 1:7)
  colnames(mat) <- paste0(base_name, 1:7)
  as.data.frame(mat)
}

test.avg <- initialize_matrix("steps.zone.corr")
test_interaction.avg <- initialize_matrix("steps.zone.corr")
test_SDinteraction.avg <- initialize_matrix("steps.zone.corr")
test_pvalinteraction.avg <- initialize_matrix("steps.zone.corr")

test.p <- initialize_matrix("steps.pvalue.corr")
test_interaction.p <- initialize_matrix("steps.pvalue.corr")
test_SDinteraction.p <- initialize_matrix("steps.pvalue.corr")
test_pvalinteraction.p <- initialize_matrix("steps.pvalue.corr")

# compute average correlation values and p-values
for (i in 1:7)
{
  for (j in 1:7)
  {
    ii <- (7*c(0:perm_num_minus))+i
    test.avg[i,j] <- mean(test.combined[ii,j])
    tt <- t.test(test.combined[ii,j], mu = 0)
    test.p[i,j] <- tt$p.value
    
    test_interaction.avg[i,j] <- mean(test_interaction.combined[ii,j])
    tt <- t.test(test_interaction.combined[ii,j], mu = 0)
    test_interaction.p[i,j] <- tt$p.value
    
  }
}

# corr coeffs
temp1 <- test.avg
temp1_interaction <- test_interaction.avg

# upper triangle removal (avoid symmetric entries)
tri.temp <- upper.tri(temp1, diag = FALSE)
tri_interaction.temp <- upper.tri(temp1_interaction, diag = FALSE)
temp1[tri.temp] <- NA
temp1_interaction[tri_interaction.temp] <- NA

# label rows and columns
x <- c("12am-5am","5am-8am","8am-11am","11am-2pm", "2pm-5pm", "5pm-9pm", "9pm-12am")
colnames(temp1) <- x
colnames(temp1_interaction) <- x
x <- c("12am-5am","5am-8am","8am-11am","11am-2pm", "2pm-5pm", "5pm-9pm", "9pm-12am")
rownames(temp1) <- x
rownames(temp1_interaction) <- x
temp1$var <- rownames(temp1)
temp1_interaction$var <- rownames(temp1_interaction)

# melt into long format and rename melted columns 
melt.temp1 = melt(temp1, id = "var")
melt_interaction.temp1 = melt(temp1_interaction, id = "var")
x <- c("CGM", "steps", "corr_coeff")
colnames(melt.temp1) <- x
colnames(melt_interaction.temp1) <- x

# factorize time bins to preserve order
melt.temp1$CGM <- factor(melt.temp1$CGM, levels=c("12am-5am","5am-8am","8am-11am","11am-2pm", "2pm-5pm", "5pm-9pm", "9pm-12am"))
melt_interaction.temp1$CGM <- factor(melt_interaction.temp1$CGM, levels=c("12am-5am","5am-8am","8am-11am","11am-2pm", "2pm-5pm", "5pm-9pm", "9pm-12am"))


# p-values
temp2 <- test.p
temp2_interaction <- test_interaction.p

# upper triangle removal (avoid symmetric entries)
tri.temp <- upper.tri(temp2, diag = FALSE)
tri_interaction.temp <- upper.tri(temp2_interaction, diag = FALSE)
temp2[tri.temp] <- NA
temp2_interaction[tri_interaction.temp] <- NA

# label rows and columns
x <- c("12am-5am","5am-8am","8am-11am","11am-2pm", "2pm-5pm", "5pm-9pm", "9pm-12am")
colnames(temp2) <- x
colnames(temp2_interaction) <- x
x <- c("12am-5am","5am-8am","8am-11am","11am-2pm", "2pm-5pm", "5pm-9pm", "9pm-12am")
rownames(temp2) <- x
rownames(temp2_interaction) <- x
temp2$var <- rownames(temp2)
temp2_interaction$var <- rownames(temp2_interaction)

# melt into long format and rename melted columns 
melt.temp2 = melt(temp2, id = "var")
melt_interaction.temp2 = melt(temp2_interaction, id = "var")
x <- c("CGM", "steps", "perm.p.value")
colnames(melt.temp2) <- x
colnames(melt_interaction.temp2) <- x

# factorize time bins to preserve order
melt.temp2$CGM <- factor(melt.temp2$CGM, levels=c("12am-5am","5am-8am","8am-11am","11am-2pm", "2pm-5pm", "5pm-9pm", "9pm-12am"))
melt_interaction.temp2$CGM <- factor(melt_interaction.temp2$CGM, levels=c("12am-5am","5am-8am","8am-11am","11am-2pm", "2pm-5pm", "5pm-9pm", "9pm-12am"))

# combine corr coef and p-value dataframes and format 
melt.cgm.steps.perm = left_join(melt.temp1,melt.temp2,by=c("CGM"="CGM", "steps"="steps"))
melt.cgm.steps.perm$corr_coeff <- round(melt.cgm.steps.perm$corr_coeff,digits=2) 
melt.cgm.steps.perm$perm.p.value[melt.cgm.steps.perm$perm.p.value > 0.001] <- NA
melt.cgm.steps.perm$perm.p.value[melt.cgm.steps.perm$perm.p.value < 0.001] <- "*"
melt.cgm.steps.perm$perm.p.value[!is.na(melt.cgm.steps.perm$perm.p.value)] <- formatC(melt.cgm.steps.perm$perm.p.value[!is.na(melt.cgm.steps.perm$perm.p.value)], format = "e", digits = 0)

# interaction now 
melt_interaction.cgm.steps.perm = left_join(melt_interaction.temp1,melt_interaction.temp2,by=c("CGM"="CGM", "steps"="steps"))
melt_interaction.cgm.steps.perm$corr_coeff <- round(melt_interaction.cgm.steps.perm$corr_coeff,digits=2)
melt_interaction.cgm.steps.perm$perm.p.value[melt_interaction.cgm.steps.perm$perm.p.value > 0.001] <- NA
melt_interaction.cgm.steps.perm$perm.p.value[melt_interaction.cgm.steps.perm$perm.p.value < 0.001] <- "*"
melt_interaction.cgm.steps.perm$perm.p.value[!is.na(melt_interaction.cgm.steps.perm$perm.p.value)] <- formatC(melt_interaction.cgm.steps.perm$perm.p.value[!is.na(melt_interaction.cgm.steps.perm$perm.p.value)], format = "e", digits = 0)

###################################################
###################################################
#steps v. cgm next 24-48h
###################################################
###################################################

# initialize lists
test <- list() 
test_interaction <- list()
all_predictions <- list()

df <- data.frame(matrix("", ncol = 7, nrow = 7))  
x <- c("steps.zone.corr1", "steps.zone.corr2", "steps.zone.corr3", "steps.zone.corr4", "steps.zone.corr5", "steps.zone.corr6", "steps.zone.corr7")
colnames(df) <- x

df <- as.data.frame(sapply(df, as.numeric))

x <- c("CGM.zone1", "CGM.zone2", "CGM.zone3", "CGM.zone4", "CGM.zone5", "CGM.zone6", "CGM.zone7")
rownames(df) <- x

# iterate through the permuations 
for (t in 1:perm_num)
{
  sampling <- data.sub[1,]
  sampling <- sampling[-1,]
  for (u in (unique(data.sub$ID))) # iterate through the participants 
    
  {
    data.temp <- data.sub %>% subset(ID == u) # data for the individual
    if (nrow(data.temp) < sample_size){next} # insufficient data, given sample size
    # extract sample of the data 
    sampling_index = nrow(data.temp) - sample_size + 1
    start_index <- sample(1:sampling_index, 1)
    end_index <- start_index + sample_size - 1
    sampling <- rbind(sampling,data.temp[start_index:end_index, ]) # sample of data for the individual
  }
  
  # initialize for the permutation
  test[[t]] <- df
  test_interaction[[t]] <- df
  
  a <- nrow(sampling)/sample_size
  ii <- sample_size*c(0:(a-1)) + 1
  jj <- sample_size*c(1:a)
  
  # initialize an empty dataframe for storing all predictions
  all_predictions[[t]] <- data.frame()

  j_labels <- c("steps: 12am-5am", "steps: 5am-8am", "steps: 8am-11am", "steps: 11am-2pm", "steps: 2pm-5pm", "steps: 5pm-9pm", "steps: 9pm-12am")
  i_labels <- c("12am-5am+24h", "5am-8am+24h", "8am-11am+24h", "11am-2pm+24h", "2pm-5pm+24h", "5pm-9pm+24h", "9pm-12am+24h")
  
  for (i in 1:7)
  {
    #i for CGM
    for (j in 1:7)
    {
      #j for steps
      temp <- sampling
      #skip if more than half of the cgm or carb values are zero for that time-window
      if (sum(temp[,i+num_cgm_start] == 0) > nrow(temp)/4*3){next}
      if (sum(temp[,j+num_steps_start] == 0) > nrow(temp)/4*3){next}
      temp[ii,i+num_cgm_start] <- NA
      temp[jj,j+num_steps_start] <- NA
      cgm.temp <- temp[,i+num_cgm_start]
      cgm.temp <- cgm.temp[!is.na(cgm.temp)]
      steps.temp <- temp[,j+num_steps_start]
      steps.temp <- steps.temp[!is.na(steps.temp)]
      res <- cor.test(cgm.temp, steps.temp, method = "pearson")
      test[[t]][i,j] <- res$estimate
      
      # set up the linear model
      numeric_indices <- c(i+num_cgm_start, j+num_steps_start)
      selected_columns <- temp[, numeric_indices]
      selected_columns <- cbind(selected_columns, temp[,"sspg_status"])
      names(selected_columns)[1] <- "CGM"
      names(selected_columns)[2] <- "steps"
      
      # fit the model 
      model <- lm(CGM ~ steps + sspg_status + steps:sspg_status, data = selected_columns)
      model_summary <- summary(model)
      coefficients <- as.data.frame(model_summary$coefficients)
      test_interaction[[t]][i,j] <- coefficients[4,1] # estimate for interaction term
      #test_SDinteraction[[t]][i,j] <- coefficients[4,2] # std. error for interaction
      #test_pvalinteraction[[t]][i,j] <- coefficients[4,4] # p-value for interaction
      
      # create a new data frame for prediction
      steps_range <- seq(0, 0.8, length.out = 100)
      sspg_status_levels <- na.omit(unique(selected_columns$sspg_status))
      prediction_data <- expand.grid(steps = steps_range, sspg_status = sspg_status_levels)
      prediction_data$CGM_predicted <- predict(model, newdata = prediction_data)
      
      prediction_data$i_index <- i
      prediction_data$j_index <- j
      
      # add custom labels for the current iteration
      prediction_data$i_label <- i_labels[i]
      prediction_data$j_label <- j_labels[j]
      
      # combine with the all_predictions dataframe
      all_predictions[[t]] <- rbind(all_predictions[[t]], prediction_data)
      
    }}
}

combined_df <- do.call(rbind, all_predictions)

# now calculate the mean of the predictions for each combination
average_predictions_24_48 <- combined_df %>%
  group_by(steps, sspg_status, i_index, j_index, i_label, j_label) %>%
  summarize(CGM_predicted = mean(CGM_predicted, na.rm = TRUE))

j_labels_reversed = rev(j_labels)
average_predictions_24_48$i_label <- factor(average_predictions_24_48$i_label, levels = i_labels)
average_predictions_24_48$j_label <- factor(average_predictions_24_48$j_label, levels = j_labels_reversed)

# combine correlation matrices from permutations
test.combined <- test[[1]][1,]
test.combined <- test.combined[-1,]

test_interaction.combined <- test_interaction[[1]][1,]
test_interaction.combined <- test_interaction.combined[-1,]

# combine the results across permutations 
for (i in 1:perm_num)
{
  test.combined <- rbind(test.combined,test[[i]])
  test_interaction.combined <- rbind(test_interaction.combined,test_interaction[[i]])
}

# avg correlation matrices from permutations
test.avg <- data.frame(matrix("", ncol = 7, nrow = 7)) 
test_interaction.avg <- data.frame(matrix("", ncol = 7, nrow = 7)) 
x <- c("steps.zone.corr1", "steps.zone.corr2", "steps.zone.corr3", "steps.zone.corr4", "steps.zone.corr5", "steps.zone.corr6", "steps.zone.corr7")
colnames(test.avg) <- x
colnames(test_interaction.avg) <- x
test.avg <- as.data.frame(sapply(test.avg, as.numeric))
test_interaction.avg <- as.data.frame(sapply(test_interaction.avg, as.numeric))
x <- c("CGM.zone1", "CGM.zone2", "CGM.zone3", "CGM.zone4", "CGM.zone5", "CGM.zone6", "CGM.zone7")
rownames(test.avg) <- x
rownames(test_interaction.avg) <- x

test.p <- data.frame(matrix("", ncol = 7, nrow = 7))  
test_interaction.p <- data.frame(matrix("", ncol = 7, nrow = 7))  
x <- c("steps.pvalue.corr1", "steps.pvalue.corr2", "steps.pvalue.corr3", "steps.pvalue.corr4", "steps.pvalue.corr5", "steps.pvalue.corr6", "steps.pvalue.corr7")
colnames(test.p) <- x
colnames(test_interaction.p) <- x
test.p <- as.data.frame(sapply(test.p, as.numeric))
test_interaction.p <- as.data.frame(sapply(test_interaction.p, as.numeric))
x <- c("CGM.zone1", "CGM.zone2", "CGM.zone3", "CGM.zone4", "CGM.zone5", "CGM.zone6", "CGM.zone7")
rownames(test.p) <- x
rownames(test_interaction.p) <- x

# compute average correlation values and p-values
for (i in 1:7)
{
  for (j in 1:7)
  {
    ii <- (7*c(0:perm_num_minus))+i
    test.avg[i,j] <- mean(test.combined[ii,j])
    tt <- t.test(test.combined[ii,j], mu = 0)
    test.p[i,j] <- tt$p.value
    
    test_interaction.avg[i,j] <- mean(test_interaction.combined[ii,j])
    tt <- t.test(test_interaction.combined[ii,j], mu = 0)
    test_interaction.p[i,j] <- tt$p.value
  }
}

# corr coeffs
temp1 <- test.avg
temp1_interaction <- test_interaction.avg
x <- c("12am-5am","5am-8am","8am-11am","11am-2pm", "2pm-5pm", "5pm-9pm", "9pm-12am")
colnames(temp1) <- x
colnames(temp1_interaction) <- x
x <- c("12am-5am+24h","5am-8am+24h","8am-11am+24h","11am-2pm+24h", "2pm-5pm+24h", "5pm-9pm+24h", "9pm-12am+24h")
rownames(temp1) <- x
rownames(temp1_interaction) <- x
temp1$var <- rownames(temp1)
temp1_interaction$var <- rownames(temp1_interaction)

# melt into long format and rename melted columns 
melt.temp1 = melt(temp1, id = "var")
melt_interaction.temp1 = melt(temp1_interaction, id = "var")
x <- c("CGM", "steps", "corr_coeff")
colnames(melt.temp1) <- x
colnames(melt_interaction.temp1) <- x

# factorize time bins to preserve order
melt.temp1$CGM <- factor(melt.temp1$CGM, levels=c("12am-5am+24h","5am-8am+24h","8am-11am+24h","11am-2pm+24h", "2pm-5pm+24h", "5pm-9pm+24h", "9pm-12am+24h"))
melt_interaction.temp1$CGM <- factor(melt_interaction.temp1$CGM, levels=c("12am-5am+24h","5am-8am+24h","8am-11am+24h","11am-2pm+24h", "2pm-5pm+24h", "5pm-9pm+24h", "9pm-12am+24h"))


# p-values
temp2 <- test.p
temp2_interaction <- test_interaction.p
x <- c("12am-5am","5am-8am","8am-11am","11am-2pm", "2pm-5pm", "5pm-9pm", "9pm-12am")
colnames(temp2) <- x
colnames(temp2_interaction) <- x
x <- c("12am-5am+24h","5am-8am+24h","8am-11am+24h","11am-2pm+24h", "2pm-5pm+24h", "5pm-9pm+24h", "9pm-12am+24h")
rownames(temp2) <- x
rownames(temp2_interaction) <- x
temp2$var <- rownames(temp2)
temp2_interaction$var <- rownames(temp2_interaction)

# melt into long format and rename melted columns 
melt.temp2 = melt(temp2, id = "var")
melt_interaction.temp2 = melt(temp2_interaction, id = "var")
x <- c("CGM", "steps", "perm.p.value")
colnames(melt.temp2) <- x
colnames(melt_interaction.temp2) <- x

# factorize time bins to preserve order
melt.temp2$CGM <- factor(melt.temp2$CGM, levels=c("12am-5am+24h","5am-8am+24h","8am-11am+24h","11am-2pm+24h", "2pm-5pm+24h", "5pm-9pm+24h", "9pm-12am+24h"))
melt_interaction.temp2$CGM <- factor(melt_interaction.temp2$CGM, levels=c("12am-5am+24h","5am-8am+24h","8am-11am+24h","11am-2pm+24h", "2pm-5pm+24h", "5pm-9pm+24h", "9pm-12am+24h"))

# combine corr coef and p-value dataframes and format 
melt.cgm.steps.perm.next = left_join(melt.temp1,melt.temp2,by=c("CGM"="CGM", "steps"="steps"))
melt.cgm.steps.perm.next$corr_coeff <- round(melt.cgm.steps.perm.next$corr_coeff,digits=2)
melt.cgm.steps.perm.next$perm.p.value[melt.cgm.steps.perm.next$perm.p.value > 0.001] <- NA
melt.cgm.steps.perm.next$perm.p.value[melt.cgm.steps.perm.next$perm.p.value < 0.001] <- "*"
melt.cgm.steps.perm.next$perm.p.value[!is.na(melt.cgm.steps.perm.next$perm.p.value)] <- formatC(melt.cgm.steps.perm.next$perm.p.value[!is.na(melt.cgm.steps.perm.next$perm.p.value)], format = "e", digits = 0)

# interaction
melt_interaction.cgm.steps.perm.next = left_join(melt_interaction.temp1,melt_interaction.temp2,by=c("CGM"="CGM", "steps"="steps"))
melt_interaction.cgm.steps.perm.next$corr_coeff <- round(melt_interaction.cgm.steps.perm.next$corr_coeff,digits=2)
melt_interaction.cgm.steps.perm.next$perm.p.value[melt_interaction.cgm.steps.perm.next$perm.p.value > 0.001] <- NA
melt_interaction.cgm.steps.perm.next$perm.p.value[melt_interaction.cgm.steps.perm.next$perm.p.value < 0.001] <- "*"
melt_interaction.cgm.steps.perm.next$perm.p.value[!is.na(melt_interaction.cgm.steps.perm.next$perm.p.value)] <- formatC(melt_interaction.cgm.steps.perm.next$perm.p.value[!is.na(melt_interaction.cgm.steps.perm.next$perm.p.value)], format = "e", digits = 0)

###################################################
###################################################
#steps v. cgm next 48-72h
###################################################
###################################################

test <- list() 
test_interaction <- list()
all_predictions <- list()

df <- data.frame(matrix("", ncol = 7, nrow = 7))  
x <- c("steps.zone.corr1", "steps.zone.corr2", "steps.zone.corr3", "steps.zone.corr4", "steps.zone.corr5", "steps.zone.corr6", "steps.zone.corr7")
colnames(df) <- x

df <- as.data.frame(sapply(df, as.numeric))

x <- c("CGM.zone1", "CGM.zone2", "CGM.zone3", "CGM.zone4", "CGM.zone5", "CGM.zone6", "CGM.zone7")
rownames(df) <- x

# iterate through the permutations
for (t in 1:perm_num)
{
  sampling <- data.sub[1,]
  sampling <- sampling[-1,]
  for (u in (unique(data.sub$ID))) # iterate through the participants 
    
  {
    data.temp <- data.sub %>% subset(ID == u) # data for the individual
    if (nrow(data.temp) <= sample_size){next} # insufficient data, given sample size
    # extract sample of the data 
    sampling_index = nrow(data.temp) - sample_size + 1
    start_index <- sample(1:sampling_index, 1)
    end_index <- start_index + sample_size - 1
    sampling <- rbind(sampling,data.temp[start_index:end_index, ]) # sample of data for the individual 
  }
  
  # initialize for the permutation
  test[[t]] <- df
  test_interaction[[t]] <- df
  
  a <- nrow(sampling)/sample_size
  ii <- sample_size*c(0:(a-1)) + 1
  ii2 <- sample_size*c(0:(a-1)) + 2
  ii3 <- as.numeric(rbind(ii2,ii))
  jj <- sample_size*c(1:a)
  jj2 <- sample_size*c(1:a) - 1 
  jj3 <- as.numeric(rbind(jj2,jj))
  
  # initialize an empty dataframe for storing all predictions
  all_predictions[[t]] <- data.frame()
  j_labels <- c("steps: 12am-5am", "steps: 5am-8am", "steps: 8am-11am", "steps: 11am-2pm", "steps: 2pm-5pm", "steps: 5pm-9pm", "steps: 9pm-12am")
  i_labels <- c("12am-5am+48h", "5am-8am+48h", "8am-11am+48h", "11am-2pm+48h", "2pm-5pm+48h", "5pm-9pm+48h", "9pm-12am+48h")
  
  for (i in 1:7)
  {
    #i for CGM
    for (j in 1:7)
    {
      #j for steps
      temp <- sampling
      #skip if more than half of the cgm or carb values are zero for that time-window
      if (sum(temp[,i+num_cgm_start] == 0) > nrow(temp)/4*3){next}
      if (sum(temp[,j+num_steps_start] == 0) > nrow(temp)/4*3){next}
      temp[ii3,i+num_cgm_start] <- NA
      temp[jj3,j+num_steps_start] <- NA
      cgm.temp <- temp[,i+num_cgm_start]
      cgm.temp <- cgm.temp[!is.na(cgm.temp)]
      steps.temp <- temp[,j+num_steps_start]
      steps.temp <- steps.temp[!is.na(steps.temp)]
      res <- cor.test(cgm.temp, steps.temp, method = "pearson")
      test[[t]][i,j] <- res$estimate
      
      # fit the linear model
      numeric_indices <- c(i+num_cgm_start, j+num_steps_start)
      selected_columns <- temp[, numeric_indices]
      selected_columns <- cbind(selected_columns, temp[,"sspg_status"])
      names(selected_columns)[1] <- "CGM"
      names(selected_columns)[2] <- "steps"
      
      # fit the model 
      model <- lm(CGM ~ steps + sspg_status + steps:sspg_status, data = selected_columns)
      model_summary <- summary(model)
      coefficients <- as.data.frame(model_summary$coefficients)
      test_interaction[[t]][i,j] <- coefficients[4,1] # estimate for interaction term
      #test_SDinteraction[[t]][i,j] <- coefficients[4,2] # std. error for interaction
      #test_pvalinteraction[[t]][i,j] <- coefficients[4,4] # p-value for interaction
      
      # create a new data frame for prediction
      steps_range <- seq(0, 0.8, length.out = 100)
      sspg_status_levels <- na.omit(unique(selected_columns$sspg_status))
      prediction_data <- expand.grid(steps = steps_range, sspg_status = sspg_status_levels)
      prediction_data$CGM_predicted <- predict(model, newdata = prediction_data)
      
      prediction_data$i_index <- i
      prediction_data$j_index <- j
      
      # add custom labels for the current iteration
      prediction_data$i_label <- i_labels[i]
      prediction_data$j_label <- j_labels[j]
      
      # combine with the all_predictions dataframe
      all_predictions[[t]] <- rbind(all_predictions[[t]], prediction_data)
      
    }}
  
}

combined_df <- do.call(rbind, all_predictions)

# now calculate the mean of the predictions for each combination
average_predictions_48_72 <- combined_df %>%
  group_by(steps, sspg_status, i_index, j_index, i_label, j_label) %>%
  summarize(CGM_predicted = mean(CGM_predicted, na.rm = TRUE))

j_labels_reversed = rev(j_labels)
average_predictions_48_72$i_label <- factor(average_predictions_48_72$i_label, levels = i_labels)
average_predictions_48_72$j_label <- factor(average_predictions_48_72$j_label, levels = j_labels_reversed)

# combine correlation matrices from permutations
test.combined <- test[[1]][1,]
test.combined <- test.combined[-1,]

test_interaction.combined <- test_interaction[[1]][1,]
test_interaction.combined <- test_interaction.combined[-1,]

# combine the results across permutations 
for (i in 1:perm_num)
{
  test.combined <- rbind(test.combined,test[[i]])
  test_interaction.combined <- rbind(test_interaction.combined,test_interaction[[i]])
}

# avg correlation matrices from permutations
test.avg <- data.frame(matrix("", ncol = 7, nrow = 7)) 
test_interaction.avg <- data.frame(matrix("", ncol = 7, nrow = 7)) 
x <- c("steps.zone.corr1", "steps.zone.corr2", "steps.zone.corr3", "steps.zone.corr4", "steps.zone.corr5", "steps.zone.corr6", "steps.zone.corr7")
colnames(test.avg) <- x
colnames(test_interaction.avg) <- x
test.avg <- as.data.frame(sapply(test.avg, as.numeric))
test_interaction.avg <- as.data.frame(sapply(test_interaction.avg, as.numeric))
x <- c("CGM.zone1", "CGM.zone2", "CGM.zone3", "CGM.zone4", "CGM.zone5", "CGM.zone6", "CGM.zone7")
rownames(test.avg) <- x
rownames(test_interaction.avg) <- x

test.p <- data.frame(matrix("", ncol = 7, nrow = 7))  
test_interaction.p <- data.frame(matrix("", ncol = 7, nrow = 7))  
x <- c("steps.pvalue.corr1", "steps.pvalue.corr2", "steps.pvalue.corr3", "steps.pvalue.corr4", "steps.pvalue.corr5", "steps.pvalue.corr6", "steps.pvalue.corr7")
colnames(test.p) <- x
colnames(test_interaction.p) <- x
test.p <- as.data.frame(sapply(test.p, as.numeric))
test_interaction.p <- as.data.frame(sapply(test_interaction.p, as.numeric))
x <- c("CGM.zone1", "CGM.zone2", "CGM.zone3", "CGM.zone4", "CGM.zone5", "CGM.zone6", "CGM.zone7")
rownames(test.p) <- x
rownames(test_interaction.p) <- x

# compute average correlation values and p-values
for (i in 1:7)
{
  for (j in 1:7)
  {
    ii <- (7*c(0:perm_num_minus))+i
    test.avg[i,j] <- mean(test.combined[ii,j])
    tt <- t.test(test.combined[ii,j], mu = 0)
    test.p[i,j] <- tt$p.value
    
    test_interaction.avg[i,j] <- mean(test_interaction.combined[ii,j])
    tt <- t.test(test_interaction.combined[ii,j], mu = 0)
    test_interaction.p[i,j] <- tt$p.value
  }
}

# corr coeffs
temp1 <- test.avg
temp1_interaction <- test_interaction.avg

# upper triangle removal (avoid symmetric entries)
tri.temp <- lower.tri(temp1, diag = FALSE)
tri_interaction.temp <- lower.tri(temp1_interaction, diag = FALSE)
temp1[tri.temp] <- NA
temp1_interaction[tri_interaction.temp] <- NA

# label rows and columns
x <- c("12am-5am","5am-8am","8am-11am","11am-2pm", "2pm-5pm", "5pm-9pm", "9pm-12am")
colnames(temp1) <- x
colnames(temp1_interaction) <- x
x <- c("12am-5am+48h","5am-8am+48h","8am-11am+48h","11am-2pm+48h", "2pm-5pm+48h", "5pm-9pm+48h", "9pm-12am+48h")
rownames(temp1) <- x
rownames(temp1_interaction) <- x
temp1$var <- rownames(temp1)
temp1_interaction$var <- rownames(temp1_interaction)

# melt into long format and rename melted columns 
melt.temp1 = melt(temp1, id = "var")
melt_interaction.temp1 = melt(temp1_interaction, id = "var")
x <- c("CGM", "steps", "corr_coeff")
colnames(melt.temp1) <- x
colnames(melt_interaction.temp1) <- x

# factorize time bins to preserve order
melt.temp1$CGM <- factor(melt.temp1$CGM, levels=c("12am-5am+48h","5am-8am+48h","8am-11am+48h","11am-2pm+48h", "2pm-5pm+48h", "5pm-9pm+48h", "9pm-12am+48h"))
melt_interaction.temp1$CGM <- factor(melt_interaction.temp1$CGM, levels=c("12am-5am+48h","5am-8am+48h","8am-11am+48h","11am-2pm+48h", "2pm-5pm+48h", "5pm-9pm+48h", "9pm-12am+48h"))

# p-values
temp2 <- test.p
temp2_interaction <- test_interaction.p
tri.temp <- lower.tri(temp2, diag = FALSE)
tri_interaction.temp <- lower.tri(temp2_interaction, diag = FALSE)
temp2[tri.temp] <- NA
temp2_interaction[tri_interaction.temp] <- NA

# label rows and columns
x <- c("12am-5am","5am-8am","8am-11am","11am-2pm", "2pm-5pm", "5pm-9pm", "9pm-12am")
colnames(temp2) <- x
colnames(temp2_interaction) <- x
x <- c("12am-5am+48h","5am-8am+48h","8am-11am+48h","11am-2pm+48h", "2pm-5pm+48h", "5pm-9pm+48h", "9pm-12am+48h")
rownames(temp2) <- x
rownames(temp2_interaction) <- x
temp2$var <- rownames(temp2)
temp2_interaction$var <- rownames(temp2_interaction)

# melt into long format and rename melted columns 
melt.temp2 = melt(temp2, id = "var")
melt_interaction.temp2 = melt(temp2_interaction, id = "var")
x <- c("CGM", "steps", "perm.p.value")
colnames(melt.temp2) <- x
colnames(melt_interaction.temp2) <- x

# factorize time bins to preserve order
melt.temp2$CGM <- factor(melt.temp2$CGM, levels=c("12am-5am+48h","5am-8am+48h","8am-11am+48h","11am-2pm+48h", "2pm-5pm+48h", "5pm-9pm+48h", "9pm-12am+48h"))
melt_interaction.temp2$CGM <- factor(melt_interaction.temp2$CGM, levels=c("12am-5am+48h","5am-8am+48h","8am-11am+48h","11am-2pm+48h", "2pm-5pm+48h", "5pm-9pm+48h", "9pm-12am+48h"))

# combine corr coef and p-value dataframes and format 
melt.cgm.steps.perm.next2 = left_join(melt.temp1,melt.temp2,by=c("CGM"="CGM", "steps"="steps"))
melt.cgm.steps.perm.next2$corr_coeff <- round(melt.cgm.steps.perm.next2$corr_coeff,digits=2)
melt.cgm.steps.perm.next2$perm.p.value[melt.cgm.steps.perm.next2$perm.p.value > 0.001] <- NA
melt.cgm.steps.perm.next2$perm.p.value[melt.cgm.steps.perm.next2$perm.p.value < 0.001] <- "*"
melt.cgm.steps.perm.next2$perm.p.value[!is.na(melt.cgm.steps.perm.next2$perm.p.value)] <- formatC(melt.cgm.steps.perm.next2$perm.p.value[!is.na(melt.cgm.steps.perm.next2$perm.p.value)], format = "e", digits = 0)
melt.cgm.steps.perm.next2.comb <- melt.cgm.steps.perm.next2
melt.cgm.steps.perm.next2.comb$perm.p.value[is.na(melt.cgm.steps.perm.next2.comb$corr_coeff)] <- NA
melt.cgm.steps.perm.next2.comb <- rbind(melt.cgm.steps.perm,melt.cgm.steps.perm.next,melt.cgm.steps.perm.next2)

# factorize time bins to preserve order
melt.cgm.steps.perm.next2.comb$CGM <- factor(melt.cgm.steps.perm.next2.comb$CGM, levels=c("12am-5am","5am-8am","8am-11am","11am-2pm", "2pm-5pm", "5pm-9pm", "9pm-12am", "12am-5am+24h","5am-8am+24h","8am-11am+24h","11am-2pm+24h", "2pm-5pm+24h", "5pm-9pm+24h", "9pm-12am+24h", "12am-5am+48h","5am-8am+48h","8am-11am+48h","11am-2pm+48h", "2pm-5pm+48h", "5pm-9pm+48h", "9pm-12am+48h"))
melt.cgm.steps.perm.next2.comb <- melt.cgm.steps.perm.next2.comb %>%
  mutate(Steps = steps, correlation_coefficient = corr_coeff) %>%
  select(-steps, -corr_coeff)

# interaction 
melt_interaction.cgm.steps.perm.next2 = left_join(melt_interaction.temp1,melt_interaction.temp2,by=c("CGM"="CGM", "steps"="steps"))
melt_interaction.cgm.steps.perm.next2$corr_coeff <- round(melt_interaction.cgm.steps.perm.next2$corr_coeff,digits=2)
melt_interaction.cgm.steps.perm.next2$perm.p.value[melt_interaction.cgm.steps.perm.next2$perm.p.value > 0.001] <- NA
melt_interaction.cgm.steps.perm.next2$perm.p.value[melt_interaction.cgm.steps.perm.next2$perm.p.value < 0.001] <- "*"
melt_interaction.cgm.steps.perm.next2$perm.p.value[!is.na(melt_interaction.cgm.steps.perm.next2$perm.p.value)] <- formatC(melt_interaction.cgm.steps.perm.next2$perm.p.value[!is.na(melt_interaction.cgm.steps.perm.next2$perm.p.value)], format = "e", digits = 0)
melt_interaction.cgm.steps.perm.next2.comb <- melt_interaction.cgm.steps.perm.next2
melt_interaction.cgm.steps.perm.next2.comb$perm.p.value[is.na(melt_interaction.cgm.steps.perm.next2.comb$corr_coeff)] <- NA
melt_interaction.cgm.steps.perm.next2.comb <- rbind(melt_interaction.cgm.steps.perm,melt_interaction.cgm.steps.perm.next,melt_interaction.cgm.steps.perm.next2)

# factorize time bins to preserve order
melt_interaction.cgm.steps.perm.next2.comb$CGM <- factor(melt_interaction.cgm.steps.perm.next2.comb$CGM, levels=c("12am-5am","5am-8am","8am-11am","11am-2pm", "2pm-5pm", "5pm-9pm", "9pm-12am", "12am-5am+24h","5am-8am+24h","8am-11am+24h","11am-2pm+24h", "2pm-5pm+24h", "5pm-9pm+24h", "9pm-12am+24h", "12am-5am+48h","5am-8am+48h","8am-11am+48h","11am-2pm+48h", "2pm-5pm+48h", "5pm-9pm+48h", "9pm-12am+48h"))
melt_interaction.cgm.steps.perm.next2.comb <- melt_interaction.cgm.steps.perm.next2.comb %>%
  mutate(Steps = steps, correlation_coefficient = corr_coeff) %>%
  select(-steps, -corr_coeff)

num_cases = nrow(sampling) / sample_size

###################################################
###################################################
# visualizations
###################################################
###################################################

p_title <- paste("all cohort")

# plot heatmap
dodge <- position_dodge(width=0.9)

plts = ggplot(melt.cgm.steps.perm.next2.comb, aes(Steps, CGM, fill= correlation_coefficient)) + 
  geom_tile() +
  coord_equal() +
  geom_text(aes(label=perm.p.value), size = 6, fontface = "bold") +
  scale_fill_gradient2(low="blue", mid="white", high="red", 
                       midpoint=0,    
                       limits=c(-1, 1)) + #limits=c(-0.5, 0.5)) +
  ylab("CGM (mean)") +
  xlab("Step Count") +
  ggtitle(p_title) +
  theme(panel.border = element_rect(colour = "black", fill=NA, size=0.5), 
        panel.background = element_rect(fill="white", colour="white", size=0.1,linetype="solid", color="black"),
        panel.grid.major = element_line(size = 0.1, linetype = 'solid', colour = "grey"),
        panel.grid.minor = element_line(size = 0.1, linetype = 'solid', colour = "grey"), 
        axis.ticks.x = element_line(colour="black"),
        axis.ticks.y = element_line(colour="black"),
        axis.line = element_line(colour = 'black', size = 0.1), 
        axis.title=element_text(size=20,face="bold"), 
        plot.title = element_text(size=21, face="bold"),
        axis.text.x = element_text(colour = "black", size=16, angle = 90, vjust = 0.5, hjust = 1, face = "bold"), 
        axis.text.y = element_text(colour = "black",size=16, face = "bold"),
        legend.title = element_text(size=14,face = "bold"),
        legend.text = element_text(size=14,face = "bold")) 

names(melt_interaction.cgm.steps.perm.next2.comb)[names(melt_interaction.cgm.steps.perm.next2.comb) == "correlation_coefficient"] <- "interaction_estimate"
plts_interaction <- ggplot(melt_interaction.cgm.steps.perm.next2.comb, aes(Steps, CGM, fill = interaction_estimate)) + 
  geom_tile() +
  coord_equal() +
  geom_text(aes(label=perm.p.value), size = 6, fontface = "bold") +
  scale_fill_gradient2(low="blue", mid="white", high="red", midpoint=0) +
  ylab("CGM (mean)") +
  xlab("Step Count") +
  ggtitle(p_title) +
  theme(panel.border = element_rect(colour = "black", fill=NA, size=0.5), 
        panel.background = element_rect(fill="white", colour="white", size=0.1, linetype="solid", color="black"),
        panel.grid.major = element_line(size = 0.1, linetype = 'solid', colour = "grey"),
        panel.grid.minor = element_line(size = 0.1, linetype = 'solid', colour = "grey"), 
        axis.ticks.x = element_line(colour="black"),
        axis.ticks.y = element_line(colour="black"),
        axis.line = element_line(colour = 'black', size = 0.1), 
        axis.title=element_text(size=20,face="bold"), 
        plot.title = element_text(size=21, face="bold"),
        axis.text.x = element_text(colour = "black", size=16, angle = 90, vjust = 0.5, hjust = 1, face = "bold"), 
        axis.text.y = element_text(colour = "black",size=16, face = "bold"),
        legend.title = element_text(size=14,face = "bold"),
        legend.text = element_text(size=14,face = "bold"))


# create the plot
combined_predictions <- rbind(average_predictions_0_24, average_predictions_24_48, average_predictions_48_72)
# add in significance
combined_predictions <- combined_predictions %>%
  mutate(perm.p.value = NA)
i_label_column = combined_predictions['i_label']
i_label_column <- i_label_column %>% mutate(new_label = str_replace(i_label, "CGM: ", ""))
j_label_column = combined_predictions['j_label']
j_label_column <- j_label_column %>%mutate(new_label = str_replace(j_label, "steps: ", ""))
CGM_new_label <- i_label_column$new_label
Steps_new_label <- j_label_column$new_label
combined_predictions$CGM <- CGM_new_label
combined_predictions$Steps <- Steps_new_label
combined_predictions <- left_join(combined_predictions, 
                                  melt_interaction.cgm.steps.perm.next2.comb[, c("CGM", "Steps", "perm.p.value")], 
                                  by = c("CGM", "Steps"))
combined_predictions <- combined_predictions %>%
  mutate(perm.p.value.y = ifelse(round(steps, digits = 2) == 0.50 & sspg_status == 'IR', perm.p.value.y, NA))


#excluded_indices <- data.frame(i_label = c("CGM: 12am-5am", "CGM: 12am-5am", "CGM: 12am-5am", "CGM: 12am-5am","CGM: #12am-5am","CGM: 12am-5am","CGM: 5am-8am","CGM: 5am-8am","CGM: 5am-8am","CGM: 5am-8am","CGM: 5am-8am", "CGM: #8am-11am","CGM: 8am-11am","CGM: 8am-11am","CGM: 8am-11am", "CGM: 11am-2pm","CGM: 11am-2pm","CGM: 11am-2pm","CGM: #2pm-5pm","CGM: 2pm-5pm","CGM: 5pm-9pm", "CGM: 5am-8am+48h", "CGM: 8am-11am+48h","CGM: 8am-11am+48h","CGM: #11am-2pm+48h","CGM: 11am-2pm+48h","CGM: 11am-2pm+48h","CGM: 2pm-5pm+48h","CGM: 2pm-5pm+48h","CGM: 2pm-5pm+48h","CGM: #2pm-5pm+48h","CGM: 5pm-9pm+48h","CGM: 5pm-9pm+48h","CGM: 5pm-9pm+48h","CGM: 5pm-9pm+48h","CGM: 5pm-9pm+48h","CGM: #9pm-12am+48h","CGM: 9pm-12am+48h","CGM: 9pm-12am+48h","CGM: 9pm-12am+48h","CGM: 9pm-12am+48h","CGM: 9pm-12am+48h"), #j_label = c("steps: 5am-8am","steps: 8am-11am","steps: 11am-2pm","steps: 2pm-5pm","steps: 5pm-9pm","steps: #9pm-12am","steps: 8am-11am","steps: 11am-2pm","steps: 2pm-5pm","steps: 5pm-9pm","steps: 9pm-12am","steps: #11am-2pm","steps: 2pm-5pm","steps: 5pm-9pm","steps: 9pm-12am","steps: 2pm-5pm","steps: 5pm-9pm","steps: #9pm-12am","steps: 5pm-9pm","steps: 9pm-12am","steps: 9pm-12am","steps: 12am-5am","steps: 12am-5am","steps: #5am-8am","steps: 12am-5am","steps: 5am-8am","steps: 8am-11am","steps: 12am-5am","steps: 5am-8am","steps: #8am-11am","steps: 11am-2pm","steps: 12am-5am","steps: 5am-8am","steps: 8am-11am","steps: 11am-2pm","steps: #2pm-5pm","steps: 12am-5am","steps: 5am-8am","steps: 8am-11am","steps: 11am-2pm","steps: 2pm-5pm", "steps: 5pm-9pm"))

excluded_indices <- data.frame(i_label = c("12am-5am", "12am-5am", "12am-5am", "12am-5am","12am-5am","12am-5am","5am-8am","5am-8am","5am-8am","5am-8am","5am-8am", "8am-11am","8am-11am","8am-11am","8am-11am", "11am-2pm","11am-2pm","11am-2pm","2pm-5pm","2pm-5pm","5pm-9pm", "5am-8am+48h", "8am-11am+48h","8am-11am+48h","11am-2pm+48h","11am-2pm+48h","11am-2pm+48h","2pm-5pm+48h","2pm-5pm+48h","2pm-5pm+48h","2pm-5pm+48h","5pm-9pm+48h","5pm-9pm+48h","5pm-9pm+48h","5pm-9pm+48h","5pm-9pm+48h","9pm-12am+48h","9pm-12am+48h","9pm-12am+48h","9pm-12am+48h","9pm-12am+48h","9pm-12am+48h"), j_label = c("steps: 5am-8am","steps: 8am-11am","steps: 11am-2pm","steps: 2pm-5pm","steps: 5pm-9pm","steps: 9pm-12am","steps: 8am-11am","steps: 11am-2pm","steps: 2pm-5pm","steps: 5pm-9pm","steps: 9pm-12am","steps: #11am-2pm","steps: 2pm-5pm","steps: 5pm-9pm","steps: 9pm-12am","steps: 2pm-5pm","steps: 5pm-9pm","steps: 9pm-12am","steps: 5pm-9pm","steps: 9pm-12am","steps: 9pm-12am","steps: 12am-5am","steps: 12am-5am","steps: 5am-8am","steps: 12am-5am","steps: 5am-8am","steps: 8am-11am","steps: 12am-5am","steps: 5am-8am","steps: 8am-11am","steps: 11am-2pm","steps: 12am-5am","steps: 5am-8am","steps: 8am-11am","steps: 11am-2pm","steps: 2pm-5pm","steps: 12am-5am","steps: 5am-8am","steps: 8am-11am","steps: 11am-2pm","steps: 2pm-5pm", "steps: 5pm-9pm"))

combined_predictions_filtered <- combined_predictions %>%
  # assuming that j_label is a character representation of j_index
  # convert j_label to numeric for comparison if necessary
  # mutate(j_index = as.numeric(sub("steps: ", "", j_label))) %>%
  anti_join(excluded_indices, by = c("i_label" = "i_label", "j_label" = "j_label"))

custom_colors <- c("#FB8C00", "#558B2F")
plot_interaction <- ggplot(combined_predictions_filtered, aes(x = steps, y = CGM_predicted, color = sspg_status)) +
  geom_line(size = 1.00) +
  geom_text(data = combined_predictions_filtered, aes(label = perm.p.value.y, x = 0.4, y = Inf),  # Adjust y for better visibility
            color = "black", size = 10, vjust = 1.1) +
  scale_color_manual(values = custom_colors) + # set custom colors
  facet_grid(j_label~i_label, switch = "y", scales = "free") +#facet_grid(j_label~i_label,switch = "x", scales = "free") + # Using custom labels for facets
  scale_y_continuous("CGM (mean; mg/dL)", position="right") +   # put the y-axis labels on the right
  labs(title = "Interaction Effect of Steps and SSPG Status on CGM",
       #x = "Step Count",
       #y = "Predicted CGM (Mean)",
       color = "SSPG Status") +
  theme_minimal() +
  theme(panel.border = element_rect(colour = "black", fill=NA, size=0.5), 
        panel.background = element_rect(fill="white", colour="white", size=0.1, linetype="solid", color="black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        axis.ticks.x = element_line(colour="black"),
        axis.ticks.y = element_line(colour="black"),
        axis.line = element_line(colour = 'black', size = 0.1), 
        axis.title=element_text(size=20,face="bold"), 
        plot.title = element_text(size=21, face="bold"),
        axis.text.x = element_text(colour = "black", size=8, angle = 90, vjust = 0.5, hjust = 1, face = "bold"), 
        axis.text.y = element_text(colour = "black",size=8, face = "bold"),
        strip.text.x = element_text(size=16, face = "bold", angle = 90, hjust = 0),  
        strip.text.y = element_text(size=16, face = "bold", angle = 0, hjust = 0),
        strip.placement = "outside",  # Place strips outside the panels
        strip.background = element_blank(), 
        plot.margin = unit(c(3, 1, 1, 1), "lines"), # increase top margin
        legend.title = element_text(size=16,face = "bold"),
        legend.text = element_text(size=16,face = "bold"),
        legend.key.size = unit(1.5, 'lines'),
        legend.position = "top", # coordinates (0, 1) place the legend at the top left
        legend.justification = c(0, 1)) # increase legend key size)

# save the plot
ggsave("interaction_plot_grid_all_finalized.pdf", plot = plot_interaction, width = 15, height = 8, dpi = 300)

j_labels_to_include <- c("steps: 12am-5am", "steps: 8am-11am", "steps: 11am-2pm", "steps: 2pm-5pm")

# filter the DataFrame based on j_label values
filtered_df <- combined_predictions_filtered[combined_predictions_filtered$j_label %in% j_labels_to_include, ]
plot_interaction <- ggplot(filtered_df, aes(x = steps, y = CGM_predicted, color = sspg_status)) +
  geom_line(size = 1.00) +
  scale_color_manual(values = custom_colors) + # set custom colors
  facet_grid(j_label~i_label,switch = "y", scales = "free") + # using custom labels for facets
  scale_y_continuous("CGM (mean; mg/dL)", position="right") +   # put the y-axis labels on the right
  geom_text(data = filtered_df, aes(label = perm.p.value.y, x = 0.4, y = Inf),#x = 0.4, y = Inf),  # adjust y for better visibility
            color = "black", size = 10, vjust = 1.1) +
  labs(title = "Interaction Effect of Steps and SSPG Status on CGM",
       x = "Step Count",
       y = "Predicted CGM (Mean)",
       color = "SSPG Status") +
  theme_minimal() +
  theme(panel.border = element_rect(colour = "black", fill=NA, size=0.5), 
        panel.background = element_rect(fill="white", colour="white", size=0.1, linetype="solid", color="black"),
        panel.grid.major = element_blank(),#element_line(size = 0.1, linetype = 'solid', colour = "grey"),
        panel.grid.minor = element_blank(),#element_line(size = 0.1, linetype = 'solid', colour = "grey"), 
        axis.ticks.x = element_line(colour="black"),
        axis.ticks.y = element_line(colour="black"),
        axis.line = element_line(colour = 'black', size = 0.1), 
        axis.title=element_text(size=20,face="bold"), 
        plot.title = element_text(size=21, face="bold"),
        axis.text.x = element_text(colour = "black", size=8, angle = 90, vjust = 0.5, hjust = 1, face = "bold"), 
        axis.text.y = element_text(colour = "black",size=8, face = "bold"),
        strip.text.x = element_text(size=16, face = "bold", angle = 90, hjust = 0),  # horizontal for x
        strip.text.y = element_text(size=16, face = "bold", angle = 0, hjust = 0), # vertical for y
        #strip.text = element_text(size=16, face = "bold", angle = 0),
        strip.text = element_text(size=16, face = "bold"),
        strip.placement = "outside",  # place strips outside the panels
        strip.background = element_blank(), 
        legend.title = element_text(size=16,face = "bold"),
        legend.text = element_text(size=16,face = "bold"),
        legend.key.size = unit(1.5, 'lines')) # increase legend key size)

# save the plot
ggsave("interaction_plot_grid_filtered_finalized.pdf", plot = plot_interaction, width = 15, height = 8, dpi = 300)

#export plots
plots_arrange <- marrangeGrob(plots, nrow=1, ncol=1)

