###### Parameters for each model configuration, to be used as input files for fitting models on HPC server --------
para.forhpc <- function(){
  count = 0
  for(n in c(40, 200, 1000)){
    # for(type in c(paste0("Sparse-Ex", seq(1,4)), paste0("Dense-Ex", seq(1,2)))){
    for(type in c(paste0("VS-Ex", seq(1,3)), paste0("GR-Ex", seq(4,6))) ){
      for(p in c(12, n/2, n-1)){
        for(snr in c(0.2, 1, 8.5)){
          # if(type == "Exponential"){
          #   count = count + 1
          #   write.table(c(type, n, p, snr),
          #               file = paste(base, "/code/run_model/para_forhpc/",count,".txt",sep=''), sep = "", col.names = FALSE, row.names = FALSE)
          # }else{
          for(rho in c(0, 0.5, 0.9)){
            count = count + 1
            write.table(c(type, n, p, snr, rho),
                        file = paste(base, "/code/run_model/para_forhpc/",count,".txt",sep=''), sep = "", col.names = FALSE, row.names = FALSE)
          }
          #}
        }
      }
    }
  }
  
  for(n in c(10, 40)){
    for(type in c(paste0("GR-Ex", seq(1,3)))){
      for(p in c(6)){
        for(snr in c(0.2, 1, 8.5)){
          # if(type == "Exponential"){
          #   count = count + 1
          #   write.table(c(type, n, p, snr),
          #               file = paste(base, "/code/run_model/para_forhpc/",count,".txt",sep=''), sep = "", col.names = FALSE, row.names = FALSE)
          # }else{
          for(rho in c(0, 0.5, 0.9)){
            count = count + 1
            write.table(c(type, n, p, snr, rho),
                        file = paste(base, "/code/run_model/para_forhpc/",count,".txt",sep=''), sep = "", col.names = FALSE, row.names = FALSE)
          }
          #}
        }
      }
    }
  }
}
# base = "/Volumes/HDD/Dropbox/Sen/Research/Model_selection/RAICc"
# para.forhpc()
###### Fit the models on HPC server --------
base = "/scratch/st1864/RAICc" # the basis of the directory to store all the outputs
## Environment parameters 
args = (commandArgs(TRUE))
randomx = as.logical( as.character(args[1]) )
# restriction_general = as.logical( as.character(args[2]) )
num = as.numeric(Sys.getenv("SLURM_ARRAY_TASK_ID"))

## create folders to store the results (both temporary and final results)
for(x_type in c("fixedx", "randomx")){
  # for(problem_type in c("subset_selection", "general_restriction")){
  #   dir.create(paste(base, "results", x_type, problem_type, sep="/"), showWarnings = FALSE, recursive = TRUE)
  #   for(tmp_folder in c("betahat", "result_intermediate")){
  #     dir.create(paste(base, "tmp", x_type, problem_type, tmp_folder, sep="/"), showWarnings = FALSE, recursive = TRUE)
  #   }
  # }
  dir.create(paste(base, "results", x_type, sep="/"), showWarnings = FALSE, recursive = TRUE)
  for(tmp_folder in c("betahat", "result_intermediate")){
    dir.create(paste(base, "tmp", x_type, tmp_folder, sep="/"), showWarnings = FALSE, recursive = TRUE)
  }
}

## Read the parametesrs
# if(restriction_general){
#   para = c(t(read.table(file = paste(base, "/para_forhpc/general_restriction/", num, ".txt",sep=''), sep = "", header=FALSE)))
# }else{
#   para = c(t(read.table(file = paste(base, "/para_forhpc/subset_selection/", num, ".txt",sep=''), sep = "", header=FALSE)))
# }
para = c(t(read.table(file = paste(base, "/para_forhpc/", num, ".txt",sep=''), sep = "", header=FALSE)))
type = para[1]
n = as.numeric(para[2])
p = as.numeric(para[3])
snr = as.numeric(para[4])
# if(type != "Exponential"){
rho = as.numeric(para[5])
# }else{
#   rho = NULL
# }

## Source all the required libraries and functions
## The code and results are stored at different file systems on the server
source(paste0(getwd(), "/utils.R"))

# # Existing results
# allfiles = list.files(paste0(base, "/tmp"), pattern = paste0(filename, "_rep") )
# rep = sub(paste0(filename, "_rep"), "", allfiles)
# rep = max( as.numeric( sub(".rds", "", rep) ) )
# result_existing = readRDS(paste0(base, "/tmp/", filename, "_rep", rep, ".rds"))

ptm = proc.time()
#result = run.all(n, p, rho, snr, type, nrep=1000, random.x=randomx, restriction.general=restriction_general, write.to.file=TRUE)
result = run.all(n, p, rho, snr, type, nrep=1000, random.x=randomx, write.to.file=TRUE)
proc.time() - ptm



# type = "Sparse-Ex2"
# n = 1000
# p = 1000
# snr = 8.5
# rho = 0.5
# 
# calc.ic.all <- function(fit, y, df, sigma.sq = NULL) {
#   y = matrix(y, ncol = 1)
#   df = matrix(df, nrow = 1)
#   if (ncol(fit) != ncol(df)) {
#     stop("the number of coef vectors does not match the number of df")
#   }
#   n = nrow(y)
#   rss = colSums(sweep(fit, 1, y, "-")^2)
#   
#   ic = list()
#   #ic$aic = n*log(rss/n) + 2*(df + 1) + n
#   #ic$bic = n*log(rss/n) + log(n)*(df + 1) + n
#   
#   #ic$cp = rss/n + 2*sigma.sq*df/n
#   #ic$rcp = rss/n + sigma.sq*df*(2 + (df + 1)/(n - df - 1))/n
#   #ic$cphat = rss/n + 2*rss*df/((n-df)*n)
#   ic$sp = (n - 1)*rss/((n - df) * (n - df - 1))
#   ic$gcv = n*rss/(n - df)^2
#   
#   # df[which(df >= n - 2)] = n - 3
#   #ic$aicc = n*log(rss/n) + n*(n + df)/(n - df - 2)
#   #ic$raicc = n*log(rss/n) + n^2*(n - 1)/((n - df - 1)*(n - df - 2))
#   
#   return(ic)
# }
# run.all <- function(n, p, rho, snr, type, nrep, random.x=TRUE, restriction.general=FALSE, write.to.file=FALSE){
#   
#   ## Specify filenames, filepaths and etc
#   if(random.x){
#     seed_x = (nrep+1):(2*nrep)
#     filepath = "randomx/"
#   }else{
#     seed_x = rep(nrep+1, nrep)
#     filepath = "fixedx/"
#   }
#   snr_all = c(0.2, 1, 8.5)
#   names(snr_all) = c("lsnr", "msnr", "hsnr")
#   snr_name = names(snr_all)[which(snr_all == snr)]
#   
#   if(is.null(rho)){
#     filename = paste0(type, "_n", n, "_p", p, "_", snr_name)
#   }else{
#     filename = paste0(type, "_n", n, "_p", p, "_", snr_name, "_rho", gsub("[.]","",as.character(rho)))
#   }
#   
#   if(!restriction.general){
#     filepath = paste0(filepath, "subset_selection")
#     ## Generate beta and Sigma
#     beta_Sigma = gen.beta.Sigma(p, rho, type)
#     m = p:0
#   }else{
#     filepath = paste0(filepath, "general_restriction")
#     ## Generate beta and Sigma
#     beta_Sigma = gen.beta.Sigma.restriction(p, rho, type)
#     m = unlist( lapply(beta_Sigma$restriction_candidate$R, nrow) )
#   }
#   
#   ## Check if some intermediate results already exist
#   result = betahat_all = list()
#   # if(paste0(filename, ".rds") %in% list.files(paste0(base, "/tmp/", filepath, "/result_intermediate"))){
#   #   result = readRDS(paste0(base, "/tmp/", filepath, "/result_intermediate/", filename, ".rds"))
#   # }
#   rep_start = length(result) + 1
#   
#   for(rep in rep_start:nrep){
#     ## Generate data
#     data = gen.data(n, p, snr, type, beta_Sigma$beta, beta_Sigma$Sigma.half, seed.y=rep, seed.x=seed_x[rep])
#     
#     i_allmethods = list()
#     
#     betahat = coef.nest(data$x, data$y, intercept=FALSE)
#     ## Coefficient vectors for all candidates and corresponding fitted values
#     muhat = data$x %*% betahat
#     
#     ## Information criteria
#     ic_all = calc.ic.all(muhat, data$y, df = p-m , sigma.sq = NULL) 
#     i_ic_all = lapply(ic_all, which.min)
#     names(i_ic_all) = paste("ic", names(i_ic_all), sep="_")
#     i_allmethods = c(i_allmethods, i_ic_all)
#     
#     ## Evaluate the selection rules
#     result[[rep]] = eval.metrics(unlist(i_allmethods), muhat, betahat, data$x, data$y, data$sigma, data$mu, beta_Sigma$beta, beta_Sigma$Sigma, m)
#     betahat_all[[rep]] = betahat
#     
#     if(rep %% 100 == 0){
#       print(rep)
#       ## Save temporary results
#       if(write.to.file){
#         saveRDS(result, paste0(base, "/tmp/", filepath, "/result_intermediate/", filename, ".rds"))
#         saveRDS(betahat_all[(rep-100+1):rep], paste0(base, "/tmp/", filepath, "/betahat/", filename, "_rep", rep, ".rds"))
#         betahat_all = list()
#       }
#     }
#   }
#   
#   ## Adjust the layout of results for the final output
#   allmethods = strsplit(names(i_allmethods), split="_")
#   result_final = list()
#   for(i in 1:length(allmethods)){
#     tmp = do.call(cbind, lapply(1:nrep, function(rep){unlist(lapply(result[[rep]], "[[", i))  }))
#     if(length(allmethods[[i]]) == 1){
#       result_final[[ allmethods[[i]] ]] = split(tmp, row(tmp, as.factor=TRUE))
#     }else{
#       result_final[[allmethods[[i]][1]]][[allmethods[[i]][2]]] = split(tmp, row(tmp, as.factor=TRUE))
#     }
#   }
#   if(write.to.file){
#     saveRDS(result_final, paste0(base, "/results/", filepath, "/", filename, ".rds"))
#   }
#   
#   return(result_final)
# }
# 
# ptm = proc.time()
# result = run.all(n, p, rho, snr, type, nrep=1000, random.x=FALSE, restriction.general=FALSE, write.to.file=FALSE)
# proc.time() - ptm
# 
# boxplot(result$ic$gcv$nvar, result$ic$sp$nvar)
# 
# filename = paste0(type, "_n", n, "_p", 996, "_", "hsnr", "_rho", gsub("[.]","",as.character(rho)))
# result = readRDS(paste0(base_results, "/", type.x, "/", type.problem, "/", filename, ".rds"))


# boxplot( result$ic$raicc$lossF, result$ic$aicc$lossF )
# boxplot( sqrt(result$ic$raicc$lossF), sqrt(result$ic$aicc$lossF) )
# boxplot( log(result$ic$raicc$lossF), log(result$ic$aicc$lossF) )
# mean (log(result$ic$raicc$lossF))
# mean (log(result$ic$aicc$lossF))

# snr_all = c(0.2, 1, 8.5)
# names(snr_all) = c("lsnr", "msnr", "hsnr")
# 
# type.x = "randomx"
# type.x = "fixedx"
# allfiles = list.files(paste0("/scratch/st1864/RAICc/results/", type.x, "/general_restriction"))
# 
# missing = c()
# for(num in 1:378){
#   para = c(t(read.table(file = paste("/scratch/st1864/RAICc/para_forhpc/general_restriction/", num, ".txt",sep=''), sep = "", header=FALSE)))
#   type = para[1]
#   n = as.numeric(para[2])
#   p = as.numeric(para[3])
#   snr = as.numeric(para[4])
#   rho = as.numeric(para[5])
#   snr_name = names(snr_all)[which(snr_all == snr)]
# 
#   filename = paste0(type, "_n", n, "_p", p, "_", snr_name, "_rho", gsub("[.]","",as.character(rho)))
# 
#   if(!paste0(filename, ".rds") %in% allfiles){
#     missing = c(missing, num)
#   }else{
#     test = readRDS(paste0("/scratch/st1864/RAICc/results/", type.x, "/general_restriction/", filename, ".rds"))
#   }
# }
# paste(missing, collapse = ',')


#
# tmp = Diagonal(p)
# R = lapply(1:nrow(tmp), function(i){tmp[i:p,,drop=FALSE]})
# r = lapply(R, function(RR){Matrix(rep(0, nrow(RR)), ncol=1)})



# test = run.all(n, p, rho, snr, type, nrep=20)
# ptm = proc.time()
# result = run.all(n, p, rho, snr, type, nrep=100, random.x=TRUE, restriction.general = TRUE, write.to.file=FALSE)
# proc.time() - ptm




