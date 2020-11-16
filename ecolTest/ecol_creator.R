setwd("C:/Users/user/Desktop/R/ecolTest/")
library(roxygen2)
library("devtools")
create("ecolTest")
setwd("./ecolTest/")
document()
setwd("..")
install("ecolTest")

library("ecolTest")

help(Hutcheson.T.test)

