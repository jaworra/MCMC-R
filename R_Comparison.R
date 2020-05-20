#R Matrix Testing

#Initialise matrix
matris1=matrix(rnorm(10),1000,1000)
matris2=matrix(rnorm(10),1000,1000)


# RNative
ptm <- proc.time() # Start the clock!
matris3 = matris1 * matris2
R_matrix_multiplication  = proc.time() - ptm
R_matrix_multiplication  = R_matrix_multiplication[3]
R_matrix_multiplication
#-----------------------------------

# RLoop Through list
# Loop through the vector, adding one
ptm <- proc.time()
for (i in 1:100000){
  h[i] <- g[i] + 1
}
R_matrix_multiplication = proc.time() - ptm
R_matrix_multiplication = R_matrix_multiplication[3]
R_matrix_multiplication
#-----------------------------------



# lapply Multiply all values by 10:
#with lapply
ptm <- proc.time()
lapply_matrix = lapply(matris1,function(x) 10 * x)
R_matrix_lapply = proc.time() - ptm
R_matrix_lapply = R_matrix_lapply[3]
#-----------------------------------


# loop through each value and add 10:
#with for loop
mat = matris1
ptm <- proc.time()
# Create the loop with r and c to iterate over the matrix
for (r in 1:nrow(mat))   
  for (c in 1:ncol(mat))  
    mat[r,c] = mat[r,c] * 10
R_matrix_for = proc.time() - ptm
R_matrix_for
#-----------------------------------


# sapply Multiply all values by 10:
#with sapply
ptm <- proc.time()
sapply_matrix = sapply(matris1,function(x) 10 * x)
R_matrix_sapply = proc.time() - ptm
R_matrix_sapply
#-----------------------------------


# parallelize computations using a cluster.
# equivalet with parrellel sapply, Multiply all values by 10:
library("parallel")

#------------2cores-----------------------
cl = makeCluster(2) #set number of cluster objects
ptm <- proc.time()
parSapply_parrallel_matrix = parSapply(cl,matris1,function(x) 10 * x)

stopCluster(cl) #always stop cluter, avoid memory leak
parSapply2_parrallel_matrix_time = proc.time() - ptm
parSapply2_parrallel_matrix_time


#------------3cores-----------------------
cl = makeCluster(3) #set number of cluster objects
ptm <- proc.time()
parSapply_parrallel_matrix = parSapply(cl,matris1,function(x) 10 * x)

stopCluster(cl) #always stop cluter, avoid memory leak
parSapply3_parrallel_matrix_time = proc.time() - ptm
parSapply3_parrallel_matrix_time


#------------4cores-----------------------
cl = makeCluster(4) #set number of cluster objects
ptm <- proc.time()
parSapply_parrallel_matrix = parSapply(cl,matris1,function(x) 10 * x)

stopCluster(cl) #always stop cluter, avoid memory leak
parSapply4_parrallel_matrix_time = proc.time() - ptm
parSapply4_parrallel_matrix_time



#Running jobs in parallel incurs overhead. Only if the jobs you fire at the worker nodes take a significant amount of time does parallelization improve overall performance. 
#-----------------------------------

results_survey <- data.frame("for_loop" = c(R_matrix_for),
                     "sapply" = c(R_matrix_sapply),
                     "parallel_2" = c(parSapply2_parrallel_matrix_time),
                     "parallel_3" = c(parSapply3_parrallel_matrix_time),
                     "parallel_4" = c(parSapply4_parrallel_matrix_time)
                     )

results_survey<- head(results_survey,-2)#remove user.child and sys.child recored
results_survey


#three decimal places
is.num <- sapply(results_survey, is.numeric)
results_survey[is.num] <- lapply(results_survey[is.num], round, 3)
library(gridExtra)
#table output
grid.table(results_survey)



x <- c(1,10,100,1000,10000,20000)
#sapply
x <- c(1,10,100,1000,2000)
i=0
l = 0
remove(l)
l <- list();
for (val in x) {
  i = i+1
  matris_temp =  matrix(rnorm(10),val,val)
  
  ptm <- proc.time()
  sapply_matrix = sapply(matris_temp,function(x) 10 * x)
  R_matrix_sapply = proc.time() - ptm
  l[i] <- R_matrix_sapply[3]
}
remove(matris_temp)
remove(sapply_matrix)
sapply_matrix_l <- l



x <- c(1,10,100,1000,2000)
i=0
l = 0
remove(l)
l <- list();

for (val in x) {
  i = i+1
  matris_temp =  matrix(rnorm(10),val,val)
  
  cl = makeCluster(2) #set number of cluster objects
  ptm <- proc.time()
  parSapply_parrallel_matrix = parSapply(cl,matris_temp,function(x) 10 * x)
  
  stopCluster(cl) #always stop cluter, avoid memory leak
  parSapply2_parrallel_matrix_time = proc.time() - ptm

  l[i] <- parSapply2_parrallel_matrix_time[3]
}
remove(matris_temp)
remove(sapply_matrix)
parSapply <- l


results_survey<- head(results_survey,-2)



x <- c(1,10,100,1000,2000)
i=0
l = 0
remove(l)
l <- list();
for (val in x) {
  i = i+1
  matris_temp =  matrix(rnorm(10),val,val)
  
  cl = makeCluster(4) #set number of cluster objects
  ptm <- proc.time()
  parSapply_parrallel_matrix = parSapply(cl,matris_temp,function(x) 10 * x)
  
  stopCluster(cl) #always stop cluter, avoid memory leak
  parSapply4_parrallel_matrix_time = proc.time() - ptm
  
  l[i] <- parSapply4_parrallel_matrix_time[3]
}
remove(matris_temp)
remove(sapply_matrix)
parSapply4 <- l















a
is.num <- sapply(R_matrix_sapply, is.numeric)
R_matrix_sapply[is.num] <- lapply(R_matrix_sapply[is.num], round, 3)

R_matrix_sapply



#R_matrix_sapply[3][1]

#matris1=matrix(rnorm(10),1000,1000)
#.atris2=matrix(rnorm(10),1000,1000)

a.list <- lapply(R_matrix_sapply[3], unname)

R_matrix_sapply[3]


matris_temp

R_matrix_sapply



# Run x number of array 

# GPU 
#library(tensorflow)

# Rewrite in C/C++


