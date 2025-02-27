library(ggplot2)
library(grid)
library(gridExtra)
library(ggsignif)

## Specify the directories 
base = getwd()
base_plots = paste0(base, '/paper/figures')
base_results = paste0(base, "/code/run_model/results") # the directory where simulation results are stored
## Create the directory 
dir.create(base_plots, recursive = TRUE)

## Functions to be used throughout this file 
source(paste0(base, '/code/utils.R'))

## Extract the legend of a ggplot object
g_legend <- function(a.gplot){
  tmp <- ggplot_gtable(ggplot_build(a.gplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)
}

## Boxplot in the main text, for variable selection problem
plot.subsetselection <- function( snr_name, type.x ){
  ## parameters
  n = c(40, 1000)
  p = n - 1
  
  type = c("VS-Ex2", "VS-Ex3")
  type_name = c("Sparse", "Dense")
  rho = 0.5
  nrep = 1000
  
  ## evaluation metrics 
  if(type.x == "fixedx"){
    metrics = c("sqrtlossF", "logKLF", "nvar")
    metrics_names = c("RMSE for fixed-X (RMSEF)", "KL for fixed-X in log scale (logKLF)", "Size of the Selected Subset")
  }else{
    metrics = c("sqrtlossR", "logKLR", "nvar")
    metrics_names = c("RMSE for random-X (RMSER)", "KL for random-X in log scale (logKLR)", "Size of the Selected Subset")
  }
  metrics_brief = c("RMSE", "log(KL)", "# Variables")
  
  ## selection rules
  methods = do.call(c, list(lapply(c("raicc", "aicc", "rcp", "cp"), function(xx){c("ic", xx)}),
                            lapply(c("tenfold", "loo"), function(xx){c("cv", xx)}),
                            lapply(c("sp", "cphat", "gcv"), function(xx){c("ic", xx)}),
                            lapply(c("bic"), function(xx){c("ic", xx)})))
  methods_name = c('"RAICc"', '"AICc"', "RC[p]", "C[p]", '"10FCV"', '"LOOCV"', "S[p]", "FPE", '"GCV"', '"BIC"')
  groups = c(1, 1, 2, 2, 3, 3, 4, 4, 4, 5)
  
  ## create the data for the plot
  pp = list()
  for(i in 1:length(n)){
    for(j in 1:length(type)){
      ## simulation results
      filename = paste0(type[j], "_n", n[i], "_p", p[i], "_", snr_name, "_rho", gsub("[.]","",as.character(rho)))
      result = readRDS(paste0(base_results, "/", type.x, "/", filename, ".rds"))
      result_target = lapply(methods, function(xx){result[[xx]]})
      tmp_function <- function(xx){
        xx$sqrtlossF = sqrt( xx$lossF )
        xx$sqrtlossR = sqrt( xx$lossR )
        xx$logKLF = log( xx$KLF )
        xx$logKLR = log( xx$KLR )
        xx
      }
      result_target = lapply(result_target, tmp_function)
      ## create the data frame for plot
      df_toplot = data.frame( do.call(cbind, lapply(metrics, function(metric){unlist(lapply(result_target, function(xx){xx[[metric]]}))})) )
      colnames(df_toplot) = metrics
      df_toplot$ic.type = factor( rep(methods_name, each=nrep), levels=methods_name )
      df_toplot$group = factor( rep( groups, each=nrep) )
      
      ## make the plot
      for(k in 1:length(metrics)){
        ## calculate the mean of the evaluation metric for each criterion
        means = aggregate(formula(paste0(metrics[k], "~ic.type")), df_toplot, mean)
        means$group = factor( groups )
        
        p1 = ggplot(df_toplot, aes_string(x="ic.type", y=metrics[k], color="group")) 
        if(metrics[k] != "nvar"){
          means[,metrics[k]] = paste0(round(means[,metrics[k]], 2))
          
          p1 = p1 + geom_boxplot() 
          ## Signed rank test
          p1 = p1 + geom_signif( comparisons = list(c('"RAICc"', '"AICc"')), test = "wilcox.test", map_signif_level = FALSE,
                                 test.args=list(alternative = "two.sided", paired=TRUE, exact=FALSE), step_increase = .1,
                                 color = "black", textsize=5)
          if(metrics_brief[k] == "RMSE" & type.x == "randomx"){
            y_mean = min(df_toplot[,metrics[k]]) - 2 * sd(df_toplot[,metrics[k]])
          }else{
            y_mean = min(df_toplot[,metrics[k]]) - sd(df_toplot[,metrics[k]])
          }
        }else{
          means[,metrics[k]] = paste0(round(means[,metrics[k]], 1))
          
          p1 = p1 + geom_jitter(size = 1.2)
          if(p[i] == 999){
            y_mean = -50
          }else if(p[i] == 39){
            y_mean = -5
          }
          p1 = p1 + scale_y_continuous(
            breaks = trunc( seq(0, max(df_toplot$nvar), length.out=5) ) )
        }

        p1 = p1 + geom_text(data = means, aes_string(label = metrics[k], y = y_mean), color="black", size = 6)
        p1 = p1 + scale_x_discrete(name = "",
                                   breaks = methods_name,
                                   labels = parse(text = methods_name)
        )
        p1 = p1 + xlab("") + ylab(metrics_names[k])
        p1 = p1 + theme(legend.position = "none")
        p1 = p1 + ggtitle(paste0(type_name[j], ", n=", n[i], ", p=", p[i], ", ", metrics_brief[k]))
        p1 = p1 + theme(strip.text = element_text(size=20)) + theme(plot.title = element_text(size = 30, face = "bold", hjust=0.5), axis.text=element_text(size=15), axis.title=element_text(size=15))
        pp[[paste0(type_name[j], ", n=", n[i], ", ", metrics_brief[k])]] = p1
      }
    }
  }
  ## print to a file
  setEPS()
  postscript(file=paste0(base_plots, "/main/", type.x, "_VS_", snr_name, ".eps"), height=24, width=24)
  pp0 = grid.arrange(arrangeGrob(pp[[paste0("Sparse, n=40, RMSE")]],
                                 pp[[paste0("Sparse, n=40, log(KL)")]],
                                 pp[[paste0("Sparse, n=40, # Variables")]],
                                 pp[[paste0("Sparse, n=1000, RMSE")]],
                                 pp[[paste0("Sparse, n=1000, log(KL)")]],
                                 pp[[paste0("Sparse, n=1000, # Variables")]],
                                 pp[[paste0("Dense, n=40, RMSE")]],
                                 pp[[paste0("Dense, n=40, log(KL)")]],
                                 pp[[paste0("Dense, n=40, # Variables")]],
                                 pp[[paste0("Dense, n=1000, RMSE")]],
                                 pp[[paste0("Dense, n=1000, log(KL)")]],
                                 pp[[paste0("Dense, n=1000, # Variables")]],
                                 ncol=3))
  print(pp0)
  dev.off()
}
plot.subsetselection.maintext <- function(){
  plot.subsetselection("hsnr", "fixedx")
  plot.subsetselection("lsnr", "fixedx")
  plot.subsetselection("hsnr", "randomx")
  plot.subsetselection("lsnr", "randomx")
}
plot.subsetselection.maintext()

## Boxplot in the main text, for general restriction problem
plot.generalrestriction <- function( type.x ){
  ## parameters
  n = c(10, 40)
  p = 6
  type = c("GR-Ex1")
  snr_name = c("hsnr", "lsnr")
  names(snr_name) = c("High signal", "Low signal")
  rho = 0.5
  nrep = 1000
  
  ## evaluation metrics 
  if(type.x == "fixedx"){
    metrics = c("sqrtlossF", "logKLF", "nres")
    metrics_names = c("RMSE for fixed-X (sqrtlossF)", "KL for fixed-X in log scale (logKLF)", "Number of Restrictions for the Selected Model")
  }else{
    metrics = c("sqrtlossR", "logKLR", "nres")
    metrics_names = c("RMSE for random-X (sqrtlossR)", "KL for random-X in log scale (logKLR)", "Number of Restrictions for the Selected Model")
  }
  metrics_brief = c("RMSE", "log(KL)", "# Restrictions")
  
  ## selection rules
  methods = do.call(c, list(lapply(c("raicc", "aicc", "rcp", "cp"), function(xx){c("ic", xx)}),
                            lapply(c("tenfold", "loo"), function(xx){c("cv", xx)}),
                            lapply(c("sp", "cphat", "gcv"), function(xx){c("ic", xx)}),
                            lapply(c("bic"), function(xx){c("ic", xx)})))
  methods_name = c('"RAICc"', '"AICc"', "RC[p]", "C[p]", '"10FCV"', '"LOOCV"', "S[p]", "FPE", '"GCV"', '"BIC"')
  groups = c(1, 1, 2, 2, 3, 3, 4, 4, 4, 5)
  
  ## create the data for the plot
  pp = list()
  for(i in 1:length(n)){
    for(j in 1:length(snr_name)){
      ## simulation results
      filename = paste0(type, "_n", n[i], "_p", p, "_", snr_name[j], "_rho", gsub("[.]","",as.character(rho)))
      result = readRDS(paste0(base_results, "/", type.x, "/", filename, ".rds"))
      result_target = lapply(methods, function(xx){result[[xx]]})
      tmp_function <- function(xx){
        xx$sqrtlossF = sqrt( xx$lossF )
        xx$sqrtlossR = sqrt( xx$lossR )
        xx$logKLF = log( xx$KLF )
        xx$logKLR = log( xx$KLR )
        xx
      }
      result_target = lapply(result_target, tmp_function)
      ## create the data frame for plot
      df_toplot = data.frame( do.call(cbind, lapply(metrics, function(metric){unlist(lapply(result_target, function(xx){xx[[metric]]}))})) )
      colnames(df_toplot) = metrics
      df_toplot$ic.type = factor( rep(methods_name, each=nrep), levels=methods_name )
      df_toplot$group = factor( rep( groups, each=nrep) )
      
      ## make the plot
      for(k in 1:length(metrics)){
        ## calculate the mean of the evaluation metric for each criterion
        means = aggregate(formula(paste0(metrics[k], "~ic.type")), df_toplot, mean)
        means[,metrics[k]] = paste0(round(means[,metrics[k]], 2))
        means$group = factor( groups )
        
        p1 = ggplot(df_toplot, aes_string(x="ic.type", y=metrics[k], color="group")) 
        if(metrics[k] != "nres"){
          p1 = p1 + geom_boxplot() 
          ## Signed rank test
          p1 = p1 + geom_signif( comparisons = list(c('"RAICc"', '"AICc"')), test = "wilcox.test", map_signif_level = FALSE,
                                 test.args=list(alternative = "two.sided", paired=TRUE, exact=FALSE), step_increase = .1,
                                 color = "black", textsize = 5)
          y_mean = min(df_toplot[,metrics[k]]) - sd(df_toplot[,metrics[k]])

        }else{
          p1 = p1 + geom_jitter(size = 1.2)
          y_mean = -1
          p1 = p1 + scale_y_continuous(
            breaks = trunc( seq(0, max(df_toplot$nres)) ) )
        }

        p1 = p1 + geom_text(data = means, aes_string(label = metrics[k], y = y_mean), color="black", size = 6)
        p1 = p1 + scale_x_discrete(name = "",
                                   breaks = methods_name,
                                   labels = parse(text = methods_name)
                                   #labels = expression("RAICc", "AICc", "RC"[p], "C"[p], "10FCV", "LOOCV", "S"[p], tilde(C)[p], "GCV", "BIC")
        )
        p1 = p1 + xlab("") + ylab(metrics_names[k])
        p1 = p1 + theme(legend.position = "none")
        p1 = p1 + ggtitle(paste0("n=", n[i], ", ", "p=", p, ", ", names(snr_name)[j], ', ', metrics_brief[k]))
        p1 = p1 + theme(strip.text = element_text(size=20)) + theme(plot.title = element_text(size = 30, face = "bold", hjust=0.5), axis.text=element_text(size=15), axis.title=element_text(size=15))
        pp[[paste0("n=", n[i], ", ", names(snr_name)[j], ", ", metrics_brief[k])]] = p1
      }
    }
  }
  ## print to a file
  setEPS()
  postscript(file=paste0(base_plots, "/main/", type.x, "_", type, ".eps"), height=24, width=24)
  pp0 = grid.arrange(arrangeGrob(pp$`n=10, High signal, RMSE`,
                                 pp$`n=10, High signal, log(KL)`,
                                 pp$`n=10, High signal, # Restrictions`,
                                 pp$`n=40, High signal, RMSE`,
                                 pp$`n=40, High signal, log(KL)`,
                                 pp$`n=40, High signal, # Restrictions`,
                                 pp$`n=10, Low signal, RMSE`,
                                 pp$`n=10, Low signal, log(KL)`,
                                 pp$`n=10, Low signal, # Restrictions`,
                                 pp$`n=40, Low signal, RMSE`,
                                 pp$`n=40, Low signal, log(KL)`,
                                 pp$`n=40, Low signal, # Restrictions`,
                                 ncol=3))
  print(pp0)
  dev.off()
}
plot.generalrestriction.maintext <- function(){
  plot.generalrestriction("fixedx")
  plot.generalrestriction("randomx")
}
plot.generalrestriction.maintext()

## Boxplot in the main text, for the mixture of general restriction and variable selection
plot.generalsubset <- function( type.x ){
  ## parameters
  n = c(40, 1000)
  p = n - 1
  type = "GR-Ex4"
  snr_name = c("hsnr", "lsnr")
  names(snr_name) = c("High signal", "Low signal")
  rho = 0.5
  nrep = 1000
  
  ## evaluation metrics 
  if(type.x == "fixedx"){
    metrics = c("sqrtlossF", "logKLF", "nres")
    metrics_names = c("RMSE for fixed-X (sqrtlossF)", "KL for fixed-X in log scale (logKLF)", "Number of restrictions for the selected subset")
  }else{
    metrics = c("sqrtlossR", "logKLR", "nres")
    metrics_names = c("RMSE for random-X (sqrtlossR)", "KL for random-X in log scale (logKLR)", "Number of restrictions for the selected subset")
  }
  metrics_brief = c("RMSE", "log(KL)", "# Restrictions")
  
  ## selection rules
  methods = do.call(c, list(lapply(c("raicc", "aicc", "rcp", "cp"), function(xx){c("ic", xx)}),
                            lapply(c("tenfold", "loo"), function(xx){c("cv", xx)}),
                            lapply(c("sp", "cphat", "gcv"), function(xx){c("ic", xx)}),
                            lapply(c("bic"), function(xx){c("ic", xx)})))
  methods_name = c('"RAICc"', '"AICc"', "RC[p]", "C[p]", '"10FCV"', '"LOOCV"', "S[p]", "FPE", '"GCV"', '"BIC"')
  groups = c(1, 1, 2, 2, 3, 3, 4, 4, 4, 5)
  
  ## create the data for the plot
  pp = list()
  for(i in 1:length(n)){
    for(j in 1:length(snr_name)){
      ## simulation results
      filename = paste0(type, "_n", n[i], "_p", p[i], "_", snr_name[j], "_rho", gsub("[.]","",as.character(rho)))
      result = readRDS(paste0(base_results, "/", type.x, "/", "/", filename, ".rds"))
      result_target = lapply(methods, function(xx){result[[xx]]})
      tmp_function <- function(xx){
        xx$sqrtlossF = sqrt( xx$lossF )
        xx$sqrtlossR = sqrt( xx$lossR )
        xx$logKLF = log( xx$KLF )
        xx$logKLR = log( xx$KLR )
        xx
      }
      result_target = lapply(result_target, tmp_function)
      ## create the data frame for plot
      df_toplot = data.frame( do.call(cbind, lapply(metrics, function(metric){unlist(lapply(result_target, function(xx){xx[[metric]]}))})) )
      colnames(df_toplot) = metrics
      df_toplot$ic.type = factor( rep(methods_name, each=nrep), levels=methods_name )
      df_toplot$group = factor( rep( groups, each=nrep) )
      
      ## make the plot
      for(k in 1:length(metrics)){
        ## calculate the mean of the evaluation metric for each criterion
        means = aggregate(formula(paste0(metrics[k], "~ic.type")), df_toplot, mean)
        means$group = factor( groups )
        
        p1 = ggplot(df_toplot, aes_string(x="ic.type", y=metrics[k], color="group")) 
        if(metrics[k] != "nres"){
          means[,metrics[k]] = paste0(round(means[,metrics[k]], 2))
          
          p1 = p1 + geom_boxplot() 
          ## Signed rank test
          p1 = p1 + geom_signif( comparisons = list(c('"RAICc"', '"AICc"')), test = "wilcox.test", map_signif_level = FALSE,
                                 test.args=list(alternative = "two.sided", paired=TRUE, exact=FALSE), step_increase = .1,
                                 color = "black", textsize = 5)
          if(metrics_brief[k] == "RMSE" & type.x == "randomx"){
            y_mean = min(df_toplot[,metrics[k]]) - 2 * sd(df_toplot[,metrics[k]])
          }else{
            y_mean = min(df_toplot[,metrics[k]]) - sd(df_toplot[,metrics[k]])
          }
        }else{
          means[,metrics[k]] = paste0(round(means[,metrics[k]], 1))
          
          p1 = p1 + geom_jitter(size = 1.2)
          if(p[i] == 999){
            y_mean = -50
          }else if(p[i] == 39){
            y_mean = -5
          }else if(p[i] == 199){
            y_mean = -25
          }
          p1 = p1 + scale_y_continuous(
            breaks = trunc( seq(0, max(df_toplot$nres), length.out=5) ) )
        }

        p1 = p1 + geom_text(data = means, aes_string(label = metrics[k], y = y_mean), color="black", size = 6)
        p1 = p1 + scale_x_discrete(name = "",
                                   breaks = methods_name,
                                   labels = parse(text = methods_name)
                                   #labels = expression("RAICc", "AICc", "RC"[p], "C"[p], "10FCV", "LOOCV", "S"[p], tilde(C)[p], "GCV", "BIC")
        )
        p1 = p1 + xlab("") + ylab(metrics_names[k])
        p1 = p1 + theme(legend.position = "none")
        p1 = p1 + ggtitle(paste0("n=", n[i], ", p=", p[i], ", ", names(snr_name)[j], ", ", metrics_brief[k]))
        p1 = p1 + theme(strip.text = element_text(size=20)) + theme(plot.title = element_text(size = 25, face = "bold", hjust=0.5), axis.text=element_text(size=15), axis.title=element_text(size=15))
        pp[[paste0("n=", n[i], ", ", snr_name[j], ", ", metrics_brief[k])]] = p1
      }
    }
  }
  ## print to a file
  setEPS()
  postscript(file=paste0(base_plots, "/main/", type.x, "_", type, ".eps"), height=24, width=24)
  pp0 = grid.arrange(arrangeGrob(pp[[paste0("n=40, hsnr, RMSE")]],
                                 pp[[paste0("n=40, hsnr, log(KL)")]],
                                 pp[[paste0("n=40, hsnr, # Restrictions")]],
                                 pp[[paste0("n=1000, hsnr, RMSE")]],
                                 pp[[paste0("n=1000, hsnr, log(KL)")]],
                                 pp[[paste0("n=1000, hsnr, # Restrictions")]],
                                 pp[[paste0("n=40, lsnr, RMSE")]],
                                 pp[[paste0("n=40, lsnr, log(KL)")]],
                                 pp[[paste0("n=40, lsnr, # Restrictions")]],
                                 pp[[paste0("n=1000, lsnr, RMSE")]],
                                 pp[[paste0("n=1000, lsnr, log(KL)")]],
                                 pp[[paste0("n=1000, lsnr, # Restrictions")]],
                                 ncol=3))
  print(pp0)
  dev.off()
}
plot.generalsubset.maintext <- function(){
  plot.generalsubset("fixedx")
  plot.generalsubset("randomx")
}
plot.generalsubset.maintext()

## Boxplot in the supplemental materials
plot.subsetselection.generalsubset <- function( type.x, type, n, snr_name, rho ){
  ## parameters
  p = c(12, n/2, n-1)
  p_name = c("small", "medium", "large")
  nrep = 1000
  
  ## evaluation metrics
  if(type.x == "fixedx"){
    metrics = c("sqrtlossF", "logKLF")
    metrics_names = c("RMSE for fixed-X (RMSEF)", "KL for fixed-X in log scale (logKLF)")
  }else{
    metrics = c("sqrtlossR", "logKLR")
    metrics_names = c("RMSE for random-X (RMSER)", "KL for random-X in log scale (logKLR)")
  }
  metrics_brief = c("RMSE", "log(KL)")
  metrics_brief_unify = c(metrics_brief , "#")
  if(grepl("GR-", type)){
    metrics = c(metrics, "nres")
    metrics_names = c(metrics_names, "Number of Restrictions for the Selected Model")
    metrics_brief = c(metrics_brief, "# Restrictions")
  }else{
    metrics = c(metrics, "nvar")
    metrics_names = c(metrics_names, "Size of the Selected Subset")
    metrics_brief = c(metrics_brief, "# Variables")
  }
  
  ## selection rules
  methods = do.call(c, list(lapply(c("raicc", "aicc", "rcp", "cp"), function(xx){c("ic", xx)}),
                            lapply(c("tenfold", "loo"), function(xx){c("cv", xx)}),
                            lapply(c("sp", "cphat", "gcv"), function(xx){c("ic", xx)}),
                            lapply(c("bic"), function(xx){c("ic", xx)})))
  methods_name = c('"RAICc"', '"AICc"', "RC[p]", "C[p]", '"10FCV"', '"LOOCV"', "S[p]", "FPE", '"GCV"', '"BIC"')
  groups = c(1, 1, 2, 2, 3, 3, 4, 4, 4, 5)
  
  ## create the data for the plot
  pp = list()
  for(i in 1:length(p)){
    ## simulation results
    filename = paste0(type, "_n", n, "_p", p[i], "_", snr_name, "_rho", gsub("[.]","",as.character(rho)))
    result = readRDS(paste0(base_results, "/", type.x, "/", filename, ".rds"))
    result_target = lapply(methods, function(xx){result[[xx]]})
    tmp_function <- function(xx){
      xx$sqrtlossF = sqrt( xx$lossF )
      xx$sqrtlossR = sqrt( xx$lossR )
      xx$logKLF = log( xx$KLF )
      xx$logKLR = log( xx$KLR )
      xx
    }
    result_target = lapply(result_target, tmp_function)
    ## create the data frame for plot
    df_toplot = data.frame( 
      do.call(cbind, lapply(metrics, function(metric){unlist(lapply(result_target, function(xx){xx[[metric]]}))})) 
      )
    # df_toplot[,2] = log(df_toplot[,2])
    colnames(df_toplot) = metrics
    df_toplot$ic.type = factor( rep(methods_name, each=nrep), levels=methods_name )
    df_toplot$group = factor( rep( groups, each=nrep) )
    
    ## make the plot
    for(k in 1:length(metrics)){
      ## calculate the mean of the evaluation metric for each criterion
      means = aggregate(formula(paste0(metrics[k], "~ic.type")), df_toplot, mean)
      means[,metrics[k]] = paste0(round(means[,metrics[k]], 2))
      means$group = factor( groups )
      
      p1 = ggplot(df_toplot, aes_string(x="ic.type", y=metrics[k], color="group")) 
      if(metrics[k] %in% c("nvar", "nres")){
        p1 = p1 + geom_jitter(size = 1.2)
        y_mean = -round(p[i]/20)
      }else{
        p1 = p1 + geom_boxplot() 
        p1 = p1 + geom_signif( comparisons = list(c('"RAICc"', '"AICc"')), test = "wilcox.test", map_signif_level = FALSE,
                               test.args=list(alternative = "two.sided", paired=TRUE, exact=FALSE), step_increase = .1,
                               color = "black", textsize=5)
        y_mean = min(df_toplot[,metrics[k]]) - 2 * sd(df_toplot[,metrics[k]])
      }
      
      p1 = p1 + geom_text(data = means, aes_string(label = metrics[k], y = y_mean), color="black", size = 5)
      p1 = p1 + scale_x_discrete(name = "",
                                 breaks = methods_name,
                                 labels = parse(text = methods_name)
      )
      p1 = p1 + xlab("") + ylab(metrics_names[k])
      p1 = p1 + theme(legend.position = "none")
      p1 = p1 + ggtitle(paste0("p=", p[i], ", ", metrics_brief[k]))
      p1 = p1 + theme(strip.text = element_text(size=20)) + theme(plot.title = element_text(size = 30, face = "bold", hjust=0.5), axis.text=element_text(size=15), axis.title=element_text(size=15))
      pp[[paste0(p_name[i], "p, ", metrics_brief_unify[k])]] = p1
    }
  }
  
  ## print to a file
  setEPS()
  postscript(file=paste0(base_plots, "/supplement/", type.x, "_", type, "_n", n, "_", snr_name, "_rho", gsub("[.]","",as.character(rho)), ".eps"), height=14, width=24)
  pp0 = grid.arrange(arrangeGrob(pp[["smallp, RMSE"]],
                                 pp[["smallp, log(KL)"]],
                                 pp[["smallp, #"]],
                                 pp[["mediump, RMSE"]],
                                 pp[["mediump, log(KL)"]],
                                 pp[["mediump, #"]],
                                 pp[["largep, RMSE"]],
                                 pp[["largep, log(KL)"]],
                                 pp[["largep, #"]],
                                 ncol=3))
  print(pp0)
  dev.off()
}
plot.generalrestriction <- function( type.x, type, n, snr_name ){
  ## parameters
  p = 6
  rho = c(0, 0.5, 0.9)
  rho_name = c("small", "medium", "large")
  nrep = 1000
  
  ## evaluation metrics
  if(type.x == "fixedx"){
    metrics = c("sqrtlossF", "logKLF", "nres")
    metrics_names = c("RMSE for fixed-X (sqrtlossF)", "KL for fixed-X in log scale (logKLF)", "Number of Restrictions for the Selected Model")
  }else{
    metrics = c("sqrtlossR", "logKLR", "nres")
    metrics_names = c("RMSE for random-X (sqrtlossR)", "KL for random-X in log scale (logKLR)", "Number of Restrictions for the Selected Model")
  }
  metrics_brief = c("RMSE", "log(KL)", "# Restrictions")
  ## selection rules
  methods = do.call(c, list(lapply(c("raicc", "aicc", "rcp", "cp"), function(xx){c("ic", xx)}),
                            lapply(c("tenfold", "loo"), function(xx){c("cv", xx)}),
                            lapply(c("sp", "cphat", "gcv"), function(xx){c("ic", xx)}),
                            lapply(c("bic"), function(xx){c("ic", xx)})))
  methods_name = c('"RAICc"', '"AICc"', "RC[p]", "C[p]", '"10FCV"', '"LOOCV"', "S[p]", "FPE", '"GCV"', '"BIC"')
  groups = c(1, 1, 2, 2, 3, 3, 4, 4, 4, 5)
  
  ## create the data for the plot
  pp = list()
  for(i in 1:length(rho)){
    ## simulation results
    filename = paste0(type, "_n", n, "_p", p, "_", snr_name, "_rho", gsub("[.]","",as.character(rho[i])))
    result = readRDS(paste0(base_results, "/", type.x, "/", filename, ".rds"))
    result_target = lapply(methods, function(xx){result[[xx]]})
    tmp_function <- function(xx){
      xx$sqrtlossF = sqrt( xx$lossF )
      xx$sqrtlossR = sqrt( xx$lossR )
      xx$logKLF = log( xx$KLF )
      xx$logKLR = log( xx$KLR )
      xx
    }
    result_target = lapply(result_target, tmp_function)
    ## create the data frame for plot
    df_toplot = data.frame( 
      do.call(cbind, lapply(metrics, function(metric){unlist(lapply(result_target, function(xx){xx[[metric]]}))})) 
    )
    colnames(df_toplot) = metrics
    df_toplot$ic.type = factor( rep(methods_name, each=nrep), levels=methods_name )
    df_toplot$group = factor( rep( groups, each=nrep) )
    
    ## make the plot
    for(k in 1:length(metrics)){
      ## calculate the mean of the evaluation metric for each criterion
      means = aggregate(formula(paste0(metrics[k], "~ic.type")), df_toplot, mean)
      means[,metrics[k]] = paste0(round(means[,metrics[k]], 2))
      means$group = factor( groups )
      
      p1 = ggplot(df_toplot, aes_string(x="ic.type", y=metrics[k], color="group")) 
      if(metrics[k] != "nres"){
        p1 = p1 + geom_boxplot() 
        p1 = p1 + geom_signif( comparisons = list(c('"RAICc"', '"AICc"')), test = "wilcox.test", map_signif_level = FALSE,
                               test.args=list(alternative = "two.sided", paired=TRUE, exact=FALSE), step_increase = .1,
                               color = "black", textsize=5)
        y_mean = min(df_toplot[,metrics[k]]) - sd(df_toplot[,metrics[k]])
      }else{
        p1 = p1 + geom_jitter(size = 1.2)
        y_mean = -1
        p1 = p1 + scale_y_continuous(
          breaks = trunc( seq(0, max(df_toplot$nres)) ) )
      }
      
      p1 = p1 + geom_text(data = means, aes_string(label = metrics[k], y = y_mean), color="black", size = 5)
      p1 = p1 + scale_x_discrete(name = "",
                                 breaks = methods_name,
                                 labels = parse(text = methods_name)
                                 #labels = expression("RAICc", "AICc", "RC"[p], "C"[p], "10FCV", "LOOCV", "S"[p], tilde(C)[p], "GCV", "BIC")
      )
      p1 = p1 + xlab("") + ylab(metrics_names[k])
      p1 = p1 + theme(legend.position = "none")
      p1 = p1 + ggtitle(paste0("rho=", rho[i], ", ", metrics_brief[k]))
      p1 = p1 + theme(strip.text = element_text(size=20)) + theme(plot.title = element_text(size = 30, face = "bold", hjust=0.5), axis.text=element_text(size=15), axis.title=element_text(size=15))
      pp[[paste0(rho_name[i], "rho, ", metrics_brief[k])]] = p1
    }
  }
  ## print to a file
  setEPS()
  postscript(file=paste0(base_plots, "/supplement/", type.x, "_", type, "_n", n, "_", snr_name, ".eps"), height=14, width=24)
  pp0 = grid.arrange(arrangeGrob(pp[["smallrho, RMSE"]],
                                 pp[["smallrho, log(KL)"]],
                                 pp[["smallrho, # Restrictions"]],
                                 pp[["mediumrho, RMSE"]],
                                 pp[["mediumrho, log(KL)"]],
                                 pp[["mediumrho, # Restrictions"]],
                                 pp[["largerho, RMSE"]],
                                 pp[["largerho, log(KL)"]],
                                 pp[["largerho, # Restrictions"]],
                                 ncol=3))
  print(pp0)
  dev.off()
}
# Generate the figures, and write a .tex file, which has the LaTeX commands to include these figures
plot.supplement <- function(){
  ## Make plots in the supplement for the subset selection problem
  typex_name = c("Fixed-X", "Random-X")
  names(typex_name) = c("fixedx", "randomx")
  snr_fullname = c("high signal", "medium signal", "low signal")
  names(snr_fullname) = c("hsnr", "msnr", "lsnr")
  texfile = file(paste0(base, "/paper/plot_supplement.tex"), "w")
  count = 0
  
  for(type in c(paste0("VS-Ex", seq(1,3)), paste0("GR-Ex", seq(1,6)))){
    if(type %in% paste0("GR-Ex", seq(1,3))){
      n_all = c(10, 40)
    }else{
      n_all = c(40, 200, 1000)
    }
    for(n in n_all){
      for(snr_name in c("hsnr", "msnr", "lsnr")){
        count = count + 1
        if(count %% 9 == 0){
          print(paste0(count, "/72"))
        }
        if(type %in% paste0("GR-Ex", seq(1,3))){
          for(type.x in c("randomx", "fixedx")){
            plot.generalrestriction( type.x, type, n, snr_name )
            write(c("\\begin{figure}[!ht]",
                    "\\centering",
                    paste0("\\includegraphics[width=\\textwidth]{figures/supplement/", type.x, "_", type, "_n", n, "_", snr_name, ".eps}"),
                    paste0("\\caption{", type, ", ", "$n=", n, "$, ", snr_fullname[names(snr_fullname)==snr_name], ", and ", typex_name[names(typex_name)==type.x], ".}"),
                    "\\end{figure}"
            ), file = texfile, append=TRUE)
          }
          write("\\clearpage", file = texfile, append=TRUE)
        }else{
          for(rho in c(0, 0.5, 0.9)){
            for(type.x in c("randomx", "fixedx")){
              plot.subsetselection.generalsubset( type.x, type, n, snr_name, rho )
              write(c("\\begin{figure}[!ht]",
                      "\\centering",
                      paste0("\\includegraphics[width=\\textwidth]{figures/supplement/", type.x, "_", type, "_n", n, "_", snr_name, "_rho", gsub("[.]","",as.character(rho)), ".eps}"),
                      paste0("\\caption{", type, ", ", "$n=", n, "$, ", snr_fullname[names(snr_fullname)==snr_name], ", $\\rho=", rho, "$, and ", typex_name[names(typex_name)==type.x], ".}"),
                      "\\end{figure}"
              ), file = texfile, append=TRUE)
            }
            write("\\clearpage", file = texfile, append=TRUE)
          }
        }
      }
    }
  }
  
  close(texfile)
}
plot.supplement()
