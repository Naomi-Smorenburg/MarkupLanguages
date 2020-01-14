###################################
# Prepare R 
###################################

# set working directory
# Set this to your own device 

# set seed for replication
set.seed(1234)

# load dataset from excel 
data <- read_excel()

# inspect dataset
View(data)

# install & load required packages 
install.packages('ranger')
install.packages("ggfortify")
install.packages("survminer")
install.packages('eha')
library(survival)
library(ranger)
library(ggplot2)
library(dplyr)
library(ggfortify)
library(tidyverse)
library(survminer)
library(eha)

###################################
# Data manipulation  
###################################

# make sure dates are in the right format
data <- data %>% mutate(DATEONSET = as.Date(DATEONSET, format = '%Y-%m-%d'))
data <- data %>% mutate(DATEDX = as.Date(DATEDX, format = '%Y-%m-%d'))
data <- data %>% mutate(DATEEVENT = as.Date(DATEEVENT, format = '%Y-%m-%d'))

# add a variable indicating days from disease onset (first symptoms) to event/last follow up
data <- data %>% mutate(time_onset = DATEEVENT-DATEONSET)

# add a variable indicating days from disease diagnosis to event/last follow up
data <- data %>% mutate(time_diagnosis = DATEEVENT-DATEDX)

# check data
data

###################################
# Use the time of onset, indicating time which patients got first symptoms 
###################################

# Kaplan Meier Survival Curve
km <- with(data, Surv(time_onset, STATUS))
head(km,80)

km_fit <- survfit(Surv(time_onset, STATUS) ~ 1, data=data)
summary(km_fit, times = c(1,30,60,90*(1:30)))

# Plot KM curve 
autoplot(km_fit)

# Plot curve using ggplot 
ggsurvplot(km_fit, data = data, pval = TRUE)

###################################
# Estimate hazards
###################################

# Create survival object
surv <- Surv(time=data$time_onset, event = data$STATUS)
surv

reg_fit <- coxph(formula=surv~STATUS,data = data)
summary(reg_fit)

# Now get baseline curve
baseline <- basehaz(reg_fit)
mean(baseline$hazard)

# hazard = 0.16

###################################
# Use the diagnosis time, indicating date of diagnosis, which has more data 
###################################
# Kaplan Meier Survival Curve
km <- with(data, Surv(time_diagnosis, STATUS))
head(km,80)

km_fit <- survfit(Surv(time_diagnosis, STATUS) ~ 1, data=data)
summary(km_fit, times = c(1,30,60,90*(1:30)))

# Plot KM curve 
autoplot(km_fit)

# Plot curve using ggplot 
ggsurvplot(km_fit, data = data, pval = TRUE)


###################################
# Estimate hazards
###################################

# Create survival object
surv <- Surv(time=data$time_diagnosis, event = data$STATUS)
surv

reg_fit <- coxph(formula=surv~STATUS,data = data)
summary(reg_fit)

# Now get baseline curve
baseline <- basehaz(reg_fit)
mean(baseline$hazard)

