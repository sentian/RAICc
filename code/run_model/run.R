###### Parameters for each model configuration, to be used as input files for fitting models on HPC server --------
para.forhpc <- function(){
  count = 0
  for(n in c(40, 200, 1000)){
    for(type in c(paste0("Sparse-Ex", seq(1,4)), paste0("Dense-Ex", seq(1,2)))){
      for(p in c(12, n/2, n-4)){
        for(snr in c(0.2, 1, 8.5)){
          # if(type == "Exponential"){
          #   count = count + 1
          #   write.table(c(type, n, p, snr),
          #               file = paste(base, "/code/run_model/para_forhpc/",count,".txt",sep=''), sep = "", col.names = FALSE, row.names = FALSE)
          # }else{
          for(rho in c(0, 0.5, 0.9)){
            count = count + 1
            write.table(c(type, n, p, snr, rho),
                        file = paste(base, "/code/run_model/para_forhpc/subset_selection/",count,".txt",sep=''), sep = "", col.names = FALSE, row.names = FALSE)
          }
          #}
        }
      }
    }
  }
}
para.forhpc.general <- function(){
  count = 0
  for(n in c(10, 100, 500)){
    for(type in c(paste0("Ex", seq(1,3)))){
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
                        file = paste(base, "/code/run_model/para_forhpc/general_restriction/",count,".txt",sep=''), sep = "", col.names = FALSE, row.names = FALSE)
          }
          #}
        }
      }
    }
  }
}
# base = "/Volumes/HDD/Dropbox/Sen/Research/Model_selection/RAICc"
# para.forhpc()
# para.forhpc.general()
###### Fit the models on HPC server --------
base = "/scratch/st1864/RAICc" # the basis of the directory to store all the outputs
## Environment parameters 
args = (commandArgs(TRUE))
randomx = as.logical( as.character(args[1]) )
restriction_general = as.logical( as.character(args[2]) )
num = as.numeric(Sys.getenv("SLURM_ARRAY_TASK_ID"))

## create folders to store the results (both temporary and final results)
for(x_type in c("fixedx", "randomx")){
  for(problem_type in c("subset_selection", "general_restriction")){
    dir.create(paste(base, "results", x_type, problem_type, sep="/"), showWarnings = FALSE, recursive = TRUE)
    for(tmp_folder in c("betahat", "result_intermediate")){
      dir.create(paste(base, "tmp", x_type, problem_type, tmp_folder, sep="/"), showWarnings = FALSE, recursive = TRUE)
    }
  }
}

## Read the parametesrs
if(restriction_general){
  para = c(t(read.table(file = paste(base, "/para_forhpc/general_restriction/", num, ".txt",sep=''), sep = "", header=FALSE)))
}else{
  para = c(t(read.table(file = paste(base, "/para_forhpc/subset_selection/", num, ".txt",sep=''), sep = "", header=FALSE)))
}
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
result = run.all(n, p, rho, snr, type, nrep=1000, random.x=randomx, restriction.general=restriction_general, write.to.file=TRUE)
proc.time() - ptm

# boxplot( result$ic$raicc$lossF, result$ic$aicc$lossF )
# boxplot( sqrt(result$ic$raicc$lossF), sqrt(result$ic$aicc$lossF) )
# boxplot( log(result$ic$raicc$lossF), log(result$ic$aicc$lossF) )
# mean (log(result$ic$raicc$lossF))
# mean (log(result$ic$aicc$lossF))

# snr_all = c(0.2, 1, 8.5)
# names(snr_all) = c("lsnr", "msnr", "hsnr")
# 
# allfiles = list.files("/scratch/st1864/RAICc/results/randomx/subset_selection")
# missing = c()
# for(num in 1:486){
#   para = c(t(read.table(file = paste("/scratch/st1864/RAICc/para_forhpc/subset_selection/", num, ".txt",sep=''), sep = "", header=FALSE)))
#   type = para[1]
#   n = as.numeric(para[2])
#   p = as.numeric(para[3])
#   snr = as.numeric(para[4])
#   rho = as.numeric(para[5])
#   snr_name = names(snr_all)[which(snr_all == snr)]
# 
#   filename = paste0(type, "_n", n, "_p", p, "_", snr_name, "_rho", gsub("[.]","",as.character(rho)))
#   # 
#   # if(!paste0(filename, ".rds") %in% allfiles){
#   #   missing = c(missing, num)
#   # }
#   test = readRDS(paste0("/scratch/st1864/RAICc/results/randomx/subset_selection/", filename, ".rds"))
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


