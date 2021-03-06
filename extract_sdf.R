# This script extracts and plots information theoretic SDFs
# which are extracted using the methodology outlined in:
# Information Theoretic Approach to High Dimensional
# "Multiplicative Models: Stochastic Discount Factor
# and Treatment Effect" by Chen Qui and Taisuke AND
# An Information-Theoretic Asset Pricing Model by
# Anisha Ghosh, Christian Julliard and Alex P. Taylor

# Clear workspace
rm(list = ls())

# Check for required packages and install them #
list.of.packages <- c("xts", "tidyverse","lubridate","lbfgs","PerformanceAnalytics","devtools",
                      "pracma")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)

# Date/Time and data wrangling #
library(xts) ; library(tidyverse) ; library(lubridate) ; library(PerformanceAnalytics)
library(pracma)
# Package for solving for the SDF #
library(lbfgs)
# devtools
library(devtools)
# ============================================================ #
# Fetch sample data from Kenneth French's website #
# In this example we will use the FF 25 Size-BTM sorted portfolios #

ff_url <- paste("https://mba.tuck.dartmouth.edu/pages/faculty/ken.french/ftp/",
                "25_Portfolios_5x5_CSV.zip",sep = "")

# Create temp_file to store the file
temp_file <- tempfile()
download.file(ff_url, temp_file)
ff_factors_raw_data <- unzip(temp_file)

# Skipping the first 15 rows
ff_factors_raw_data <- read_csv(ff_factors_raw_data, skip = 15)
df = ff_factors_raw_data[1:1132,]
colnames(df)[1] <- "date"
df$date = paste(df$date,"01",sep="")
df$date = ymd(df$date)

# Convert the datatype to xts 
df = xts(df[,2:ncol(df)], order.by = df$date)
# Data in % form - convert to decimal
df = df/100
rm(ff_factors_raw_data)

# Download the rf-free rate #
ff_url <- paste("https://mba.tuck.dartmouth.edu/pages/faculty/ken.french/ftp/",
                "F-F_Research_Data_Factors_CSV.zip",sep = "")

# Create temp_file to store the file
temp_file <- tempfile()
download.file(ff_url, temp_file)
ff_factors_raw_data <- unzip(temp_file)

# Skipping the first 15 rows
rf <- read_csv(ff_factors_raw_data, skip = 2)
rf = rf[1:1132,c(1,5)]
colnames(rf)[1] <- "date"
rf$date = paste(rf$date,"01",sep="")
rf$date = ymd(rf$date)

# Convert the datatype to xts 
rf = xts(rf[,2:ncol(rf)], order.by = rf$date)
# Data in % form - convert to decimal
rf = rf/100
rm(ff_factors_raw_data)
# ============================================================ #
# Convert returns to excess returns #
df = df -  repmat(rf, 1, ncol(df))

# Extract the SDF #
# Load/source the functions direct from GitHub #

sdf_func = paste("https://github.com/Alexander-M-Dickerson/information-theoretic-sdfs/blob/main/",
                 "exponential_tilting_fn.R?raw=TRUE",sep = "")
sdf_grad = paste("https://github.com/Alexander-M-Dickerson/information-theoretic-sdfs/blob/main/",
                 "exponential_tilting_grad.R?raw=TRUE",sep = "")

devtools::source_url(sdf_func)
devtools::source_url(sdf_grad)

# Start the analysis from 1965 #
r = df["1965/"]

# Set different levels of l1 penalty for the optimization #
a = as.matrix( seq( 0.00, 0.15, 0.05) )
l1 = a[1]

no_firms_in_SDF = matrix(nrow = 1, ncol = nrow(a))
SDFs = matrix(nrow = nrow(r), ncol = nrow(a))
P = matrix(nrow = nrow(r), ncol = nrow(a))
LAMBDAs = matrix(nrow = ncol(r), ncol = nrow(a))

T1 = nrow(r)
rt = as.matrix(r)
j = 1

for (l1 in a){
  
  print(l1)
  init = matrix(0, nrow = 1, ncol = ncol(r))
  out <- lbfgs(    exponential_tilting_fn, exponential_tilting_grad, init, ffq=r, prec=0, 
                   invisible=1, orthantwise_c = l1,
                   linesearch_algorithm="LBFGS_LINESEARCH_BACKTRACKING",
                   gtol = 0.01)
  
  lambda1<-  as.matrix( out$par ) 
  LAMBDAs[,j] = lambda1
  SDFs[,j] <- (exp(r%*%lambda1)) / (sum(exp(r%*%lambda1)) * T1^-1 )
  no_firms_in_SDF[1,j] = colSums(lambda1 != 0)
  j = j + 1
  
  if (j == length(a) ){
    colnames(LAMBDAs)         = a
    colnames(SDFs)            = a
    colnames(no_firms_in_SDF) = a
  }
  
}


# Create an xts datatype for the sdfs #
r = as.xts(r)
startDate = as.POSIXct(index(r)[1], tz="UTC")
endDate = as.POSIXct(index(r)[nrow(r)], tz="UTC")
bin = seq(from = startDate,to = endDate,by = "months")
sdf = xts(x  = SDFs, order.by = index(r))
plot(sdf)
# ============================================================ #
# Check for yourself the pricing errors are 0 as per the Euler equation #
# For l1 > 0 , the pricing error will equal the l1 penalty 
# At the asset level
round( t(as.matrix(r))%*%sdf , 2 )
# ============================================================ #
