library(xtable)

## Specify the directory to save the tables
# The default directory on my machine is 'code/', change it to its parent directory 'BOSSreg/'
# base = getwd()

# base = "/Users/sentian/Dropbox/Sen/Research/Model_selection/RAICc"
base = "/Volumes/HDD/Dropbox/Sen/Research/Model_selection/RAICc"
base_tables = paste0(base, "/paper/tables")
## Specify the directory to read the results
base_results = paste0(base, "/code/run_model/results")
## Create the directory to save the tables
dir.create(paste0(base_tables, "/supplement"), recursive = TRUE, showWarnings = FALSE)

source(paste0(base, "/code/utils.R"))
## Calculate the evaluation metrics: % worse than oracle method, relative efficiency, sparsistency (extra variables)
summary.result <- function(result, rmse.oracle){
  # n = dim(data$x[[1]])[1]
  # % worse than the oracle performance (best possible BS or BOSS)
  percentloss = unlist(lapply(result[["rmse"]], function(xx){round(100*(xx / rmse.oracle - 1), 0)}))
  
  # Sparsistency (extra variables)
  nvar.sparsistency <- function(nvar_method, sparsistency_method){
    if(is.na(sparsistency_method)){
      as.character(round(nvar_method, 1))
    }else{
      # True model
      paste(round(nvar_method, 1), "(", round(sparsistency_method, 1), ")", sep="")
    }
  }
  numvar = unlist(Map(nvar.sparsistency, result$nvar, result$sparsistency))
  
  # # Relative efficiency
  # mu = data$x %*% data$beta
  # rmse_null = sqrt(sum(mu^2) / n) # null model
  # betahat_full = ginv(data$x) %*% data$y # full OLS
  # rmse_full = mean( sqrt(colSums(sweep(data$x%*%betahat_full, 1, mu,'-')^2)/n) )
  # rmse_method = unlist(result$rmse)
  # rmse_min = min(c(rmse_method, rmse_null, rmse_full))
  # efficiency = round(rmse_min/rmse_method, 2)
  
  # output = list(percentloss = percentloss, efficiency = efficiency, numvar = numvar)
  output = list(percentloss = percentloss, numvar = numvar)
  
  # Different styles of output for different types of tables
  return(output)
}
## Evaluation metrics for all scenarios considered in a table
df.summary.result <- function(n, p, snr, type, rho, methods){
  count = 0
  output = list()
  if(!is.null(rho)){
    # X is general
    for(i in 1:length(n)){
      for(j in 1:length(snr)){
        for(m in 1:length(rho)){
          for(k in 1:length(p)){
            metrics = list()
            count = count + 1
            for(l in 1:length(type)){
              # Read the simulation results
              filename = paste0(type[l], "_n", n[i], "_p", p[k], "_", names(snr)[j], "_rho", gsub("[.]","",as.character(rho[m])))
              result = readRDS(paste0(base_results, "/", filename, ".rds"))
              result_target = lapply(methods, function(xx){result[[xx]]})
         
              # Best possible candidate
              rmse_oracle = mean(result$best$rmse)
              # The evaluation metrics for the specified methods
              tmp_function <- function(xx){
                if(!is.null(xx[[which.metric]])){
                  mean(xx[[which.metric]])
                }else{
                  NA
                }
              }
              result_tosummary = list()
              for(which.metric in c("rmse", "sparsistency", "nvar")){
                result_tosummary[[which.metric]] = lapply(result_target, tmp_function)
              }
              # Calculate the three metrics to be presented in the table
              metrics[[l]] = summary.result(result_tosummary, rmse_oracle)
            }
            output[[count]] = lapply(1:2, function(ii){do.call(c, lapply(metrics, "[[", ii))})
          }
        }
      }
    }
  }else{
    # X is orthogonal
    for(i in 1:length(n)){
      for(j in 1:length(snr)){
        for(k in 1:length(p)){
          metrics = list()
          count = count + 1
          for(l in 1:length(type)){
            # Read the simulation results
            filename = paste0(type[l], "_n", n[i], "_p", p[k], "_", names(snr)[j])
            result = readRDS(paste0(base_results, "/", filename, ".rds"))
            result_target = lapply(methods, function(xx){result[[xx]]})
            
            # Best possible candidate
            rmse_oracle = mean(result$best$rmse)
            # The evaluation metrics for the specified methods
            tmp_function <- function(xx){
              if(!is.null(xx[[which.metric]])){
                mean(xx[[which.metric]])
              }else{
                NA
              }
            }
            result_tosummary = list()
            for(which.metric in c("rmse", "sparsistency", "nvar")){
              result_tosummary[[which.metric]] = lapply(result_target, tmp_function)
            }
            # Calculate the three metrics to be presented in the table
            metrics[[l]] = summary.result(result_tosummary, rmse_oracle)
          }
          output[[count]] = lapply(1:2, function(ii){do.call(c, lapply(metrics, "[[", ii))})
        }
      }
    }
  }
  output = do.call(rbind, lapply(1:2, function(ii){do.call(rbind, lapply(output, "[[", ii))}))
  return(output) 
}

## Supplement material
table.supplement <- function(scale.box = 0.52){
  # Parameters
  n = c(200, 2000)
  pnratio = c(0.1, 0.5, 0.98)
  rho = c(0, 0.5, 0.9)
  snr = c(8.5, 1, 0.2)
  names(snr) = c("hsnr", "msnr", "lsnr")
  nrep = 1000
  
  for(type in c(paste0("Sparse-Ex", 1:4), "Dense", "Omit", "Exponential")){
    title = paste0("The performance of various selection rules, ", type)
    filename = type
    
    if(type == "Exponential"){
      rho = NULL
    }
    
    output = list()
    for(i in 1:length(n)){
      output[[i]] = df.summary.result(n[i], round(pnratio*n[i]), snr, type, rho, 
                                      methods = do.call(c, list(lapply(c("raicc", "aicc", "rcp", "cp"), function(xx){c("ic", xx)}),
                                                                lapply(c("tenfold", "loo"), function(xx){c("cv", xx)}),
                                                                lapply(c("aic", "bic", "gcv"), function(xx){c("ic", xx)})))
      )
    }
    output = do.call(cbind, output)
    
    if(type == "Exponential"){
      output = cbind(rep(c("\\midrule\\multirow{3}[2]{*}{hsnr}", "", "", 
                           "\\midrule\\multirow{3}[2]{*}{msnr}", "", "",
                           "\\midrule\\multirow{3}[2]{*}{lsnr}", "", ""), 2),
                     rep(paste0("p/n=", c(0.1, 0.5, 0.98)), 6),
                     output)
      
      command = c(paste("\\toprule \n",
                        "\\multicolumn{1}{|r}{} &   & \\multicolumn{9}{c||}{$n=200$} & \\multicolumn{9}{c|}{$n=2000$} \\\\\n",
                        "\\cmidrule{3-20}\\multicolumn{1}{|c}{} &       & RAICc & AICc  & RC$_p$ & C$_p$ & 10F CV & LOO CV & AIC   & BIC   & GCV   & RAICc & AICc  & RC$_p$ & C$_p$ & 10F CV & LOO CV & AIC   & BIC   & GCV       \\\\\n",
                        "\\cmidrule{3-20}\\multicolumn{1}{|c}{} &       & \\multicolumn{18}{c|}{\\% worse than the best possible candidate} \\\\\n"),
                  paste("\\midrule \n",
                        "\\multicolumn{1}{|c}{} &       & \\multicolumn{18}{c|}{Number of variables} \\\\\n"),
                  paste("\\bottomrule \n")
                  )
      print(xtable(output,
                   align = "l|c|c|cc|cc|cc|ccc||cc|cc|cc|ccc|",  # align and put a vertical line (first "l" again represents column of row numbers)
                   label = paste0("tab:", type),
                   caption = title),
            #size = size, #Change size; useful for bigger tables "normalsize" "footnotesize"
            scalebox = scale.box,
            caption.placement = "top",
            include.rownames = FALSE, 
            include.colnames = FALSE, 
            hline.after = NULL, 
            floating = TRUE, # whether \begin{Table} should be created (TRUE) or not (FALSE)
            sanitize.text.function = force, # Important to treat content of first column as latex function
            add.to.row = list(pos = list(-1,
                                         9,
                                         nrow(output)),
                              command = command
            ),
            file = paste0(base_tables, "/supplement/", filename, ".tex"), compress = FALSE
      )
    }else{
      output = cbind(rep(c("\\midrule\\multirow{9}[6]{*}{hsnr}", "", "", "\\cmidrule{2-21}", "", "", "\\cmidrule{2-21}", "", "",
                           "\\midrule\\multirow{9}[6]{*}{msnr}", "", "", "\\cmidrule{2-21}", "", "", "\\cmidrule{2-21}", "", "",
                           "\\midrule\\multirow{9}[6]{*}{lsnr}", "", "", "\\cmidrule{2-21}", "", "", "\\cmidrule{2-21}", "", ""), 2),
                     rep(c("\\multirow{3}[2]{*}{$\\rho=0$}", "", "", "\\multirow{3}[2]{*}{$\\rho=0.5$}", "", "", "\\multirow{3}[2]{*}{$\\rho=0.9$}", "", ""), 6),
                     rep(paste0("p/n=", c(0.1, 0.5, 0.98)), 18),
                     output)
      
      command = c(paste("\\toprule \n",
                        "\\multicolumn{1}{|r}{} & \\multicolumn{1}{r}{} &       & \\multicolumn{9}{c||}{$n=200$}                                         & \\multicolumn{9}{c|}{$n=2000$} \\\\\n",
                        "\\cmidrule{4-21}\\multicolumn{1}{|c}{} & \\multicolumn{1}{c}{} &       & RAICc & AICc  & RC$_p$ & C$_p$ & 10F CV & 200F CV & AIC   & BIC   & GCV   & RAICc & AICc  & RC$_p$ & C$_p$ & 10F CV & 200F CV & AIC   & BIC   & GCV       \\\\\n",
                        "\\cmidrule{4-21}\\multicolumn{1}{|c}{} & \\multicolumn{1}{c}{} &       & \\multicolumn{18}{c|}{\\% worse than the best possible candidate} \\\\\n"),
                  paste("\\midrule \n",
                        "\\multicolumn{1}{|c}{} & \\multicolumn{1}{c}{} &       & \\multicolumn{18}{c|}{Number of variables (Sparsistency)} \\\\\n"),
                  paste("\\bottomrule \n")
                  )
      print(xtable(output,
                   align = "l|c|c|c|cc|cc|cc|ccc||cc|cc|cc|ccc|",  # align and put a vertical line (first "l" again represents column of row numbers)
                   label = paste0("tab:", type),
                   caption = title),
            #size = size, #Change size; useful for bigger tables "normalsize" "footnotesize"
            scalebox = scale.box,
            caption.placement = "top",
            include.rownames = FALSE, 
            include.colnames = FALSE, 
            hline.after = NULL, 
            floating = TRUE, # whether \begin{Table} should be created (TRUE) or not (FALSE)
            sanitize.text.function = force, # Important to treat content of first column as latex function
            add.to.row = list(pos = list(-1,
                                         27,
                                         nrow(output)),
                              command = command
            ),
            file = paste0(base_tables, "/supplement/", filename, ".tex"), compress = FALSE
      )
    }
    
  }
}
table.supplement()



## Evaluation metrics for all scenarios considered in a table
df.summary.result <- function(n, p, snr, type, rho, type.x, type.problem){
  methods = do.call(c, list(lapply(c("raicc", "aicc", "rcp", "cp"), function(xx){c("ic", xx)}),
                            lapply(c("tenfold", "loo"), function(xx){c("cv", xx)}),
                            lapply(c("sp", "cphat", "gcv", "bic"), function(xx){c("ic", xx)})))
  allmetrics = c("rmse", "rr", "rte", "pve")
  if(type.x == 'fixedx'){
    allmetrics = c(allmetrics, "fkl")
  }else{
    allmetrics = c(allmetrics, "rkl")
  }
  if(type.problem == "subset_selection"){
    allmetrics = c(allmetrics, "nvar")
  }else{
    allmetrics = c(allmetrics, "i_allmethods")
  }
  count = 0
  output = list()
  
  for(i in 1:length(n)){
    for(j in 1:length(snr)){
      for(m in 1:length(rho)){
        for(k in 1:length(p)){
          metrics = list()
          count = count + 1
          
          # Read the simulation results
          filename = paste0(type, "_n", n[i], "_p", p[k], "_", names(snr)[j], "_rho", gsub("[.]","",as.character(rho[m])))
          result = readRDS(paste0(base_results, "/", type.x, "/", type.problem, "/", filename, ".rds"))
          result_target = lapply(methods, function(xx){result[[xx]]})
          
          # The evaluation metrics for the specified methods
          tmp_function <- function(xx){
            if(which.metric == "fkl" | which.metric == "rkl"){
              round(mean(xx[[which.metric]]) / n[i], 3)
            }else if(which.metric == "i_allmethods" & type.problem == "general_restriction"){
              round(100 * sum(xx[[which.metric]] %in% seq(4,6)) / length(xx[[which.metric]]), 1)
            }else{
              round(mean(xx[[which.metric]]), 3)
            }
          }
          result_tosummary = list()
          for(which.metric in allmetrics){
            result_tosummary[[which.metric]] = unlist( lapply(result_target, tmp_function) )
          }
          output[[count]] = result_tosummary
        }
      }
    }
  }
  
  output = lapply(1:length(allmetrics), function(ii){do.call(rbind,lapply(output, "[[", ii))})
  output = cbind( do.call(rbind, output[c(1,3,5)]), do.call(rbind, output[c(2, 4, 6)]) )
  return(output) 
}
table.supplement <- function(scale.box = 0.52, type.problem=c("subset_selection", "general_restriction")){
  type.problem = match.arg(type.problem)
  
  if(type.problem == "subset_selection"){
    alltypes = c(paste0("Sparse-Ex", 1:4), paste0("Dense-Ex", 1:2))
  }else{
    alltypes = paste0("Ex", 1:3)
  }
  
  # Parameters
  rho = c(0, 0.5, 0.9)
  snr = c(8.5, 1, 0.2)
  names(snr) = c("hsnr", "msnr", "lsnr")
  nrep = 1000
  
  for(type_x in c("fixedx", "randomx")){
    for(type in alltypes){
      if(type.problem == "subset_selection"){
        for(n in c(40, 200, 1000)){
            p = c(12, n/2, n-4)
          
          if(type_x == "fixedx"){
            title = paste0("Fixed X, ", type, ", n=", n)
          }else{
            title = paste0("Random X, ", type, ", n=", n)
          }
          filename = paste0(type, "_n", n)

          output = df.summary.result(n, p, snr, type, rho, type_x, type.problem)
          output = cbind(rep(c("\\midrule\\multirow{9}[6]{*}{hsnr}", "", "", "\\cmidrule{2-23}", "", "", "\\cmidrule{2-23}", "", "",
                               "\\midrule\\multirow{9}[6]{*}{msnr}", "", "", "\\cmidrule{2-23}", "", "", "\\cmidrule{2-23}", "", "",
                               "\\midrule\\multirow{9}[6]{*}{lsnr}", "", "", "\\cmidrule{2-23}", "", "", "\\cmidrule{2-23}", "", ""), 3),
                         rep(c("\\multirow{3}[2]{*}{$\\rho=0$}", "", "", "\\multirow{3}[2]{*}{$\\rho=0.5$}", "", "", "\\multirow{3}[2]{*}{$\\rho=0.9$}", "", ""), 9),
                         rep(paste0("p=", p), 27),
                         output)
          
          command = c(paste("\\toprule \n",
                            "\\multicolumn{1}{|c}{} & \\multicolumn{1}{c}{} &       & RAICc & AICc  & RC$_p$ & C$_p$ & 10F CV & LOO CV & S$_p$ & $\\widehat{\\text{C}}_p$ & GCV   & BIC   & RAICc & AICc  & RC$_p$ & C$_p$ & 10F CV & LOO CV & S$_p$ & $\\widehat{\\text{C}}_p$ & GCV   & BIC  \\\\\n",
                            "\\cmidrule{4-23}\\multicolumn{1}{|c}{} & \\multicolumn{1}{c}{} &       & \\multicolumn{10}{c||}{RMSE}                                                   & \\multicolumn{10}{c|}{Relative Risk}       \\\\\n"),
                      paste("\\midrule \n",
                            "\\multicolumn{1}{|c}{} & \\multicolumn{1}{c}{} &       & \\multicolumn{10}{c||}{Relative Test Error}                                    & \\multicolumn{10}{c|}{Proportion of Variance Explained} \\\\\n"),
                      paste("\\midrule \n",
                            "\\multicolumn{1}{|c}{} & \\multicolumn{1}{c}{} &       & \\multicolumn{10}{c||}{Kullback-Leibler Discrepancy}                                    & \\multicolumn{10}{c|}{Number of Variables} \\\\\n"),
                      paste("\\bottomrule \n")
          )
          print(xtable(output,
                       align = "l|c|c|c|cc|cc|cc|ccc|c||cc|cc|cc|ccc|c|",  # align and put a vertical line (first "l" again represents column of row numbers)
                       # label = paste0("tab:", filename),
                       caption = title),
                #size = size, #Change size; useful for bigger tables "normalsize" "footnotesize"
                scalebox = scale.box,
                caption.placement = "top",
                include.rownames = FALSE, 
                include.colnames = FALSE, 
                hline.after = NULL, 
                floating = TRUE, # whether \begin{Table} should be created (TRUE) or not (FALSE)
                sanitize.text.function = force, # Important to treat content of first column as latex function
                add.to.row = list(pos = list(-1,
                                             27,
                                             54,
                                             nrow(output)),
                                  command = command
                ),
                file = paste0(base_tables, "/supplement/", type_x, "/", type.problem, "/", filename, ".tex"), compress = FALSE
          )
        }
      }else{
        n = c(10, 100, 500)
        p = 6
        
        if(type_x == "fixedx"){
          title = paste0("Fixed X, ", type, ", p=", p)
        }else{
          title = paste0("Random X, ", type, ", p=", p)
        }
        filename = paste0(type)
        
        output = df.summary.result(n, p, snr, type, rho, type_x, type.problem)
        output = cbind(rep(c("\\midrule\\multirow{9}[6]{*}{$n=10$}", "", "", "\\cmidrule{2-23}", "", "", "\\cmidrule{2-23}", "", "",
                             "\\midrule\\multirow{9}[6]{*}{$n=100$}", "", "", "\\cmidrule{2-23}", "", "", "\\cmidrule{2-23}", "", "",
                             "\\midrule\\multirow{9}[6]{*}{$n=500$}", "", "", "\\cmidrule{2-23}", "", "", "\\cmidrule{2-23}", "", ""), 3),
                       rep(c("\\multirow{3}[2]{*}{hsnr}", "", "", "\\multirow{3}[2]{*}{msnr}", "", "", "\\multirow{3}[2]{*}{lsnr}", "", ""), 9),
                       rep(c("$\\rho=0$", "$\\rho=0.5$", "$\\rho=0.9$"), 27),
                       output)
        
        command = c(paste("\\toprule \n",
                          "\\multicolumn{1}{|c}{} & \\multicolumn{1}{c}{} &       & RAICc & AICc  & RC$_p$ & C$_p$ & 10F CV & LOO CV & S$_p$ & $\\widehat{\\text{C}}_p$ & GCV   & BIC   & RAICc & AICc  & RC$_p$ & C$_p$ & 10F CV & LOO CV & S$_p$ & $\\widehat{\\text{C}}_p$ & GCV   & BIC  \\\\\n",
                          "\\cmidrule{4-23}\\multicolumn{1}{|c}{} & \\multicolumn{1}{c}{} &       & \\multicolumn{10}{c||}{RMSE}                                                   & \\multicolumn{10}{c|}{Relative Risk}       \\\\\n"),
                    paste("\\midrule \n",
                          "\\multicolumn{1}{|c}{} & \\multicolumn{1}{c}{} &       & \\multicolumn{10}{c||}{Relative Test Error}                                    & \\multicolumn{10}{c|}{Proportion of Variance Explained} \\\\\n"),
                    paste("\\midrule \n",
                          "\\multicolumn{1}{|c}{} & \\multicolumn{1}{c}{} &       & \\multicolumn{10}{c||}{Kullback-Leibler Discrepancy}                                    & \\multicolumn{10}{c|}{\\% of Correct Restrictions Selected} \\\\\n"),
                    paste("\\bottomrule \n")
        )
        print(xtable(output,
                     align = "l|c|c|c|cc|cc|cc|ccc|c||cc|cc|cc|ccc|c|",  # align and put a vertical line (first "l" again represents column of row numbers)
                     # label = paste0("tab:", filename),
                     caption = title),
              #size = size, #Change size; useful for bigger tables "normalsize" "footnotesize"
              scalebox = scale.box,
              caption.placement = "top",
              include.rownames = FALSE, 
              include.colnames = FALSE, 
              hline.after = NULL, 
              floating = TRUE, # whether \begin{Table} should be created (TRUE) or not (FALSE)
              sanitize.text.function = force, # Important to treat content of first column as latex function
              add.to.row = list(pos = list(-1,
                                           27,
                                           54,
                                           nrow(output)),
                                command = command
              ),
              file = paste0(base_tables, "/supplement/", type_x, "/", type.problem, "/", filename, ".tex"), compress = FALSE
        )
        
        
      }
      
    }
  }
}
table.supplement(type.problem = "subset_selection")
table.supplement(type.problem = "general_restriction")
