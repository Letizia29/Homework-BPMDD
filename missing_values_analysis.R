library(VIM)
missing_values_pd <- read.csv("missing_values_pd.csv", header = FALSE)
missing_values_hc <- read.csv("missing_values_hc.csv", header = FALSE)

missing_values_hc = t(missing_values_hc)
missing_values_pd = t(missing_values_pd)

matrixplot(missing_values_hc)
matrixplot(missing_values_pd)

barMiss(missing_values_hc)
barMiss(missing_values_pd)

aggr(missing_values_hc, numbers=TRUE, prop=FALSE)


## load data
data = read.csv("Patient_Master.csv")

aggr(data)


new_data_hc = read.csv("new_data_hc.csv")
