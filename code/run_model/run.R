###### Parameters for each model configuration, to be used as input files to fit models on the HPC server --------
para.forhpc <- function(){
  count = 0
  for(n in c(40, 200, 1000)){
    for(type in c(paste0("VS-Ex", seq(1,3)), paste0("GR-Ex", seq(4,6))) ){
      for(p in c(12, n/2, n-1)){
        for(snr in c(0.2, 1, 8.5)){
          for(rho in c(0, 0.5, 0.9)){
            count = count + 1
            write.table(c(type, n, p, snr, rho),
                        file = paste(base, "/code/run_model/para_forhpc/",count,".txt",sep=''), sep = "", col.names = FALSE, row.names = FALSE)
          }
        }
      }
    }
  }
  
  for(n in c(10, 40)){
    for(type in c(paste0("GR-Ex", seq(1,3)))){
      for(p in c(6)){
        for(snr in c(0.2, 1, 8.5)){
          for(rho in c(0, 0.5, 0.9)){
            count = count + 1
            write.table(c(type, n, p, snr, rho),
                        file = paste(base, "/code/run_model/para_forhpc/",count,".txt",sep=''), sep = "", col.names = FALSE, row.names = FALSE)
          }
        }
      }
    }
  }
}
# base = "/Volumes/HDD/Dropbox/Sen/Research/Model_selection/RAICc"
# para.forhpc()

###### Fit the models on HPC server --------
base = "/scratch/st1864/RAICc" # the basis directory to store all the outputs
## Environment parameters specified in the .sh file for HPC
args = (commandArgs(TRUE))
randomx = as.logical( as.character(args[1]) )
num = as.numeric(Sys.getenv("SLURM_ARRAY_TASK_ID"))

## create folders to store the results (both temporary and final results)
for(x_type in c("fixedx", "randomx")){
  dir.create(paste(base, "results", x_type, sep="/"), showWarnings = FALSE, recursive = TRUE)
  for(tmp_folder in c("betahat", "result_intermediate")){
    dir.create(paste(base, "tmp", x_type, tmp_folder, sep="/"), showWarnings = FALSE, recursive = TRUE)
  }
}

## Read the parametesrs
para = c(t(read.table(file = paste(base, "/para_forhpc/", num, ".txt",sep=''), sep = "", header=FALSE)))
type = para[1]
n = as.numeric(para[2])
p = as.numeric(para[3])
snr = as.numeric(para[4])
rho = as.numeric(para[5])

## Source all the required libraries and functions
## The code and results are stored at different file systems on the HPC server
source(paste0(getwd(), "/utils.R"))

## Fit the model
ptm = proc.time()
result = run.all(n, p, rho, snr, type, nrep=1000, random.x=randomx, write.to.file=TRUE)
proc.time() - ptm
