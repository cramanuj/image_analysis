## Citation:
## https://bmcbioinformatics.biomedcentral.com/articles/10.1186/s12859-016-1083-9
#######################################################
## Compute Bartlett test statistic
#######################################################
bartlettTestStat = function(dataMatrix)
{
  means = rowMeans(dataMatrix)
  k = nrow(dataMatrix)
  n = ncol(dataMatrix)
  N = n*k
  vars = apply(dataMatrix,1,var)
  var.pooled = sum(vars) / (k-1)
  numerator = (N-k) * log(var.pooled) - (n-1) * sum(log(vars))
  denom = 1 + ( ( k/(n-1) - 1/(N-k) ) / (3* (k-1) ) )
  bt = numerator/denom
  return(bt)
}

##################################################################################################
## Estimate optimum parameters for variance stabilization in microarray data.
## Input:
#     data: intesity data matrix
#     cfLow, cfHigh: lowest and highest possible values for cofactor (log scale)
#     frac: fraction of differentially expressed genes used in variance stabilization (<0 & >=1)
# output: optimum cofactor
##################################################################################################
transformData = function(data, cfLow=0, cfHigh=300, frac=1){
  if(frac>1 || frac <= 0)
    stop(" 0< frac<=1 ")
  if(cfLow>=cfHigh)
  {
    print("Warning: cfLow>=cfHigh, using default values")
    cfLow=0
    cfHigh=10
  }

  cat("====================================================================\n")
  cat("Finding optimum cofactor for asinh transformation\n")
  cat("====================================================================\n")
  cat(sprintf("%15s     %15s \n", "cofactor(log scale)", "Bartlett\'s stat"))
  cat("====================================================================\n")

  cofactors = seq(cfLow,cfHigh,1)
  bartlett = NULL
  for(cf in cofactors)
  {
    data.t = asinh(data/exp(cf))

    #find the non-differentially exprssed genes
    if(frac<1)
    {
      diff = quantile(abs(data.t[,1] - data.t[,2]), frac)
      keep.idx = which(abs(data.t[,1] - data.t[,2]) <= diff)
      data.t = data.t[keep.idx,]
    }
    bt=bartlettTestStat(data.t)
    bartlett = c(bartlett, bt)
    cat(sprintf("%10d %25.2f \n", cf, bt))
  }

  minIdx = which.min(bartlett)
  cat("\n Optimum cofactor :", sprintf("exp(%d)",cofactors[minIdx]), "\n")
  cat("====================================================================\n\n")

  plot(cofactors, bartlett, type='o', pch=16,
       xlab="Cofactors (log scale)", ylab="Bartlett's statistics",
       main = paste("Optimum cofactor: ",
                    sprintf("exp(%d)",cofactors[minIdx]), sep=""))
  points(cofactors[minIdx], bartlett[minIdx], pch=16, col='red')

  data.t = asinh(data/exp(cofactors[minIdx]))
  return(data.t)
}
