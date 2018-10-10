library(arules)

#for each cluster, loop through support and minlen, interested in the most popular combination
most.pop.combs <- function(clust_samps){
  most_popular <- c()
  for (j in 3:10){
    clust_common <- eclat(clust_samps, parameter = list(support=0.2, minlen = j, maxlen = j, tidLists=TRUE))
    if (dim(tidLists(clust_common))[1] == 0){break}
    else {
      comb <- as(tidLists(clust_common), "list")
      comb_max <- lengths(comb)
      max_ind <- which(comb_max==max(comb_max))
      if (comb_max[max_ind[1]] > 100){
        for (i in 1:length(max_ind)){
          most_popular <- c(most_popular,comb_max[max_ind[i]])
        }
      }
      else{break}
    }
  }
  return(most_popular)
}

# for cluster i, will output cli_most_pop listing most popular combinations of OTUs that belong to cluster i, numbers represent the number of samples that the combination appears in
for (i in c(1:13)){
  clust.samples.py <- scan(paste("tup_cluster",i,".txt",sep=""), what="", sep="\n") # read in tuples generated from python that parsed through all the samples and identified members of the cluster in each sample 
  clust.samples.edited <- strsplit(clust.samples.py, "[[:space:]]+")
  clust_pop <- most.pop.combs(clust.samples.edited)
  assign(paste("cl",i,"_most_pop",sep=""),clust_pop)
}