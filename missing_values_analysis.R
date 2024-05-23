library(VIM)

# Load data
new_data_hc = read.csv("new_data_hc.csv")
new_data_pd = read.csv("new_data_pd.csv")

# immagine di
matrixplot(new_data_hc)
matrixplot(new_data_pd)

# immagine di
barMiss(new_data_hc)
barMiss(new_data_pd)

# labels
names = list("NP1PTOT", "NP1RTOT", "NP2TOT", "NP3TOT", "NP4TOT", "GENETICS", "FAMILIARITY", "ETHNICITY", "SEX", "AGE", "HEIGHT", "WEIGHT", "HAND", "PRIM_DIAG")


# combinazioni di missing values
aggr(new_data_hc, labels=names, cex.axis=0.8)
aggr(new_data_pd, labels=names, cex.axis=0.8)

