# read in data. Convert to matrix form. Keep the names of the characters
# in the right order in char.names
#
# M - original data
# M.smoothed - data used for fitting model
# M.1 - first 2 books, used for validation
# M.2 - first 3 books, used for validation
# char.names -- vector of names of characters for labelling plots.

dat <- read.table("got.txt", header=T, stringsAsFactors=F)

char.names <- dat[,1]
char.names[6] <- "Jon Snow"
char.names[22] <- "Jon C"
char.names[21] <- "Quentyn"
char.names[23] <- "Melisandre"
char.names[24] <- "Barristan"

dat <- as.matrix(dat[,-c(1,7)])
dat <- matrix(as.numeric(dat), nr=nrow(dat))
M <- dat # data matrix

# split books 4 and 5
M2 <- M
M2[,4] <- (M[,4]+M[,5])*sum(M[,4])/sum(M[,4:5])
M2[,5] <- (M[,4]+M[,5])*sum(M[,5])/sum(M[,4:5])

# combine books 4 and 5
M3 <- M[,1:4]
M3[,4] <- M[,4]+M[,5]

M.smoothed <- M2
M.combined <- M3

# for validation

M.1 <- M[1:9, 1:2]
M.2 <- M[1:12, 1:3]
