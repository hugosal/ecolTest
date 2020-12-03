setwd("C:/Users/user/Desktop/R/THutcheson/Hutcheson-T-test/trunk/ecolTest")
library("roxygen2")
library("devtools")
#create("ecolTest")
setwd("./ecolTest/")
#document()
setwd("..")

install("ecolTest")

library("ecolTest")

help(Hutcheson.T.test)

Hutcheson.T.test(x=c(1,2,3,4,5),y = c(4,5,7,8,9))
