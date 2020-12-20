# Clear workspace
rm(list = ls())
# Set directory
directory = ""
setwd(directory)

# Check for required packages and install them #
list.of.packages <- c("xts", "tidyverse","lubridate","lbfgs","PerformanceAnalytics")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)

# Date/Time and data wrangling #
library(xts) ; library(tidyverse) ; library(lubridate) ; library(PerformanceAnalytics)
# Package for solving for the SDF #
library(lbfgs)


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
# XTS
df = xts(df[,2:ncol(df)], order.by = df$date)
# Data in % form - convert to decimal
df = df/100
rm(ff_factors_raw_data)

# Extract the SDF #
source("exponential_tilting_fn.R")
source("exponential_tilting_grad.R")

r = df["1965/"]

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
  
  lambda1<-  as.matrix( out$par ) # this yields the lambda's that have been solved #
  #lambda1[ lambda1 > 0 ] = 0
  LAMBDAs[,j] = lambda1
  SDFs[,j] <- (exp(r%*%lambda1)) / (sum(exp(r%*%lambda1)) * T1^-1 )
  #(exp(r%*%lambda1)) / (sum(exp(r%*%lambda1)) * T1^-1 )
  #
  #1 / (     1 +  r%*%lambda1    )
  #
  #
  #P[,j] <-  (exp(r%*%lambda1)) / (sum(exp(r%*%lambda1)) * T1^-1 )
  #r%*%lambda1    
  no_firms_in_SDF[1,j] = colSums(lambda1 != 0)
  j = j + 1
  
  
}

r = as.xts(r)
startDate = as.POSIXct(index(r)[1], tz="UTC")
endDate = as.POSIXct(index(r)[nrow(r)], tz="UTC")
bin = seq(from = startDate,to = endDate,by = "months")
sdf = xts(x  = SDFs, order.by = index(r))
plot(sdf)
