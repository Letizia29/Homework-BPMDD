library(VIM)

# Load data
new_data_hc = read.csv("new_data_hc.csv")
new_data_pd = read.csv("new_data_pd.csv")

# labels
names = c("NP1PTOT", "NP1RTOT", "NP2TOT", "NP3TOT", "NP4TOT", "GENETICS", "FAMILIARITY", "ETHNICITY", "SEX", "AGE", "HEIGHT", "WEIGHT", "HAND", "PRIM_DIAG")

# missing values patterns
matrixplot(new_data_hc, cex.axis=0.7, ylab='Subjects')
title('HC - missing values patterns')
matrixplot(new_data_pd, cex.axis=0.7, ylab='Subjects')
title('PD - missing values patterns')

# proportions of missings
colnames(new_data_hc) = names

barMiss(new_data_hc, pos=5)
barMiss(new_data_pd, pos=14)

# combinazioni di missing values
aggr(new_data_hc, labels=names, cex.axis=0.8)
aggr(new_data_pd, labels=names, cex.axis=0.8)

