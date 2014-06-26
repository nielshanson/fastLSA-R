

n <- 1250
m <- matrix(sample(0:1,n^2,replace=T),nrow=n)
system.time(mutlinksC(m,num_threads=1))
m <- matrix(sample(0:1,n^2,replace=T),nrow=n)
system.time(mutlinksC(m,num_threads=2))
m <- matrix(sample(0:1,n^2,replace=T),nrow=n)
system.time(mutlinksC(m,num_threads=4))
m <- matrix(sample(0:1,n^2,replace=T),nrow=n)
system.time(mutlinksC(m,num_threads=8))
m <- matrix(sample(0:1,n^2,replace=T),nrow=n)
system.time(mutlinksC(m,num_threads=16))

# fastLSA code
m=4
n=3
mat <- matrix(1:(m*n),nrow=m, ncol=n)
mat
fastLSA(mat)

test_mat <-read.table("time_series.txt")
fastLSA(test_mat[1:30,],d=5, minLSA=0.2, alpha=0.0001, rez=1000000)

library('Matrix')
library('igraph')
?igraph
