library(ggplot2)
library(grid)
library(gridExtra)
library(ggsignif)

# base = "/Users/sentian/Dropbox/Sen/Research/Model_selection/RAICc"
base = "/Volumes/HDD/Dropbox/Sen/Research/Model_selection/RAICc"
base_plots = paste0(base, '/paper/figures')
base_results = paste0(base, "/code/run_model/results")
## Create the directory to save the plots
dir.create(base_plots, recursive = TRUE)

### Functions to be used throughout this file --------
source(paste0(base, '/code/utils.R'))

## Extract the legend of a ggplot object
g_legend <- function(a.gplot){
  tmp <- ggplot_gtable(ggplot_build(a.gplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)
}


## Boxplot in the main text, for subset selection problem
plot.subsetselection.alltypes <- function( snr_name, type.x, type.p=c("smallp","largep") ){
  ## parameters
  n = c(40, 1000)
  type.p = match.arg(type.p)
  if(type.p == "smallp"){
    p = rep(12, 2)
  }else{
    p = n - 4
  }
  type = c("Sparse-Ex2", "Dense-Ex1")
  type_name = c("Sparse", "Dense")
  # snr_name = "hsnr"
  # type.x = "fixedx"
  type.problem = "subset_selection"
  rho = 0.5
  nrep = 1000
  
  ## evaluation metrics 
  if(type.x == "fixedx"){
    metrics = c("loglossF", "logKLF", "nvar")
    metrics_names = c("Loss for fixed-X in log scale (loglossF)", "KL for fixed-X in log scale (logKLF)", "Size of the Selected Subset")
  }else{
    metrics = c("loglossR", "logKLR", "nvar")
    metrics_names = c("Loss for random-X in log scale (loglossR)", "KL for random-X in log scale (logKLR)", "Size of the Selected Subset")
  }
  metrics_brief = c("loss", "KL", "nvar")
  
  ## selection rules
  methods = do.call(c, list(lapply(c("raicc", "aicc", "rcp", "cp"), function(xx){c("ic", xx)}),
                            lapply(c("tenfold", "loo"), function(xx){c("cv", xx)}),
                            lapply(c("sp", "cphat", "gcv"), function(xx){c("ic", xx)}),
                            lapply(c("bic"), function(xx){c("ic", xx)})))
  # methods_name = c("RAICc", "AICc", "RCp", "Cp", "10FCV", "LOOCV", "Sp", "Cptilde", "GCV", "BIC")
  methods_name = c('"RAICc"', '"AICc"', "RC[p]", "C[p]", '"10FCV"', '"LOOCV"', "S[p]", "tilde(C)[p]", '"GCV"', '"BIC"')
  groups = c(1, 1, 2, 2, 3, 3, 4, 4, 4, 5)
  
  ## create the data for the plot
  pp = list()
  for(i in 1:length(n)){
    for(j in 1:length(type)){
      ## simulation results
      filename = paste0(type[j], "_n", n[i], "_p", p[i], "_", snr_name, "_rho", gsub("[.]","",as.character(rho)))
      result = readRDS(paste0(base_results, "/", type.x, "/", type.problem, "/", filename, ".rds"))
      result_target = lapply(methods, function(xx){result[[xx]]})
      tmp_function <- function(xx){
        xx$loglossF = log( xx$lossF )
        xx$loglossR = log( xx$lossR )
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
        if(metrics[k] != "nvar"){
          p1 = p1 + geom_boxplot() 
          p1 = p1 + geom_signif( comparisons = list(c('"RAICc"', '"AICc"')), test = "wilcox.test", map_signif_level = FALSE,
                                 test.args=list(alternative = "two.sided", paired=TRUE, exact=FALSE), step_increase = .1 )
          y_mean = min(df_toplot[,metrics[k]]) - sd(df_toplot[,metrics[k]])
        }else{
          p1 = p1 + geom_jitter(size = 1.2)
          if(p[i] == 996){
            y_mean = -50
          }else if(p[i] == 36){
            y_mean = -5
          }
          p1 = p1 + scale_y_continuous(
            breaks = trunc( seq(0, max(df_toplot$nvar), length.out=5) ) )
        }
        
        # p1 = p1 + geom_signif( comparisons = list(c("RAICc", "AICc"), c("RAICc", "RCp"), c("RAICc", "LOO CV"), c("RAICc", "GCV")), test = "wilcox.test", map_signif_level = FALSE,
        #                        test.args=list(alternative = "two.sided", paired=TRUE, exact=FALSE), step_increase = .1 )
        # p1 = p1 + coord_cartesian(ylim = ylim1)

        p1 = p1 + geom_text(data = means, aes_string(label = metrics[k], y = y_mean), size = 5)
        p1 = p1 + scale_x_discrete(name = "",
                                   breaks = methods_name,
                                   labels = parse(text = methods_name)
                                   #labels = expression("RAICc", "AICc", "RC"[p], "C"[p], "10FCV", "LOOCV", "S"[p], tilde(C)[p], "GCV", "BIC")
        )
        p1 = p1 + xlab("") + ylab(metrics_names[k])
        p1 = p1 + theme(legend.position = "none")
        p1 = p1 + ggtitle(paste0(type[j], ", n=", n[i], ", p=", p[i]))
        p1 = p1 + theme(strip.text = element_text(size=20)) + theme(plot.title = element_text(size = 30, face = "bold", hjust=0.5), axis.text=element_text(size=15), axis.title=element_text(size=15))
        pp[[paste0(type_name[j], ", n=", n[i], ", ", metrics_brief[k])]] = p1
      }
    }
  }
  ## print to a file
  setEPS()
  postscript(file=paste0(base_plots, "/main/", type.x, "/", type.problem, "/", type.p, "_", snr_name, ".eps"), height=24, width=24)
  pp0 = grid.arrange(arrangeGrob(pp[[paste0("Sparse, n=40, loss")]],
                                 pp[[paste0("Sparse, n=40, KL")]],
                                 pp[[paste0("Sparse, n=40, nvar")]],
                                 pp[[paste0("Sparse, n=1000, loss")]],
                                 pp[[paste0("Sparse, n=1000, KL")]],
                                 pp[[paste0("Sparse, n=1000, nvar")]],
                                 pp[[paste0("Dense, n=40, loss")]],
                                 pp[[paste0("Dense, n=40, KL")]],
                                 pp[[paste0("Dense, n=40, nvar")]],
                                 pp[[paste0("Dense, n=1000, loss")]],
                                 pp[[paste0("Dense, n=1000, KL")]],
                                 pp[[paste0("Dense, n=1000, nvar")]],
                                 ncol=3))
  print(pp0)
  dev.off()
}
plot.subsetselection.maintext <- function(){
  plot.subsetselection.alltypes("hsnr", "fixedx", type.p="largep")
  plot.subsetselection.alltypes("lsnr", "fixedx", type.p="largep")
  plot.subsetselection.alltypes("hsnr", "randomx", type.p="largep")
  plot.subsetselection.alltypes("lsnr", "randomx", type.p="largep")
  plot.subsetselection.alltypes("hsnr", "fixedx", type.p="smallp")
  plot.subsetselection.alltypes("lsnr", "fixedx", type.p="smallp")
  plot.subsetselection.alltypes("hsnr", "randomx", type.p="smallp")
  plot.subsetselection.alltypes("lsnr", "randomx", type.p="smallp")
}
plot.subsetselection.maintext()

## Boxplot in the main text, for general restriction problem
plot.generalrestriction.allsnr <- function( type.x ){
  ## parameters
  n = c(10, 500)
  p = 6
  type = c("Ex1")
  snr_name = c("hsnr", "lsnr")
  names(snr_name) = c("High signal", "Low signal")
  # type.x = "fixedx"
  type.problem = "general_restriction"
  rho = 0.5
  nrep = 1000
  
  ## evaluation metrics 
  if(type.x == "fixedx"){
    metrics = c("loglossF", "logKLF", "nres")
    metrics_names = c("Loss for fixed-X in log scale (loglossF)", "KL for fixed-X in log scale (logKLF)", "Number of Restrictions for the Selected Model")
  }else{
    metrics = c("loglossR", "logKLR", "nres")
    metrics_names = c("Loss for random-X in log scale (loglossR)", "KL for random-X in log scale (logKLR)", "Number of Restrictions for the Selected Model")
  }
  metrics_brief = c("loss", "KL", "nres")
  
  ## selection rules
  methods = do.call(c, list(lapply(c("raicc", "aicc", "rcp", "cp"), function(xx){c("ic", xx)}),
                            lapply(c("tenfold", "loo"), function(xx){c("cv", xx)}),
                            lapply(c("sp", "cphat", "gcv"), function(xx){c("ic", xx)}),
                            lapply(c("bic"), function(xx){c("ic", xx)})))
  # methods_name = c("RAICc", "AICc", "RCp", "Cp", "10FCV", "LOOCV", "Sp", "Cptilde", "GCV", "BIC")
  methods_name = c('"RAICc"', '"AICc"', "RC[p]", "C[p]", '"10FCV"', '"LOOCV"', "S[p]", "tilde(C)[p]", '"GCV"', '"BIC"')
  groups = c(1, 1, 2, 2, 3, 3, 4, 4, 4, 5)
  
  ## create the data for the plot
  pp = list()
  for(i in 1:length(n)){
    for(j in 1:length(snr_name)){
      ## simulation results
      filename = paste0(type, "_n", n[i], "_p", p, "_", snr_name[j], "_rho", gsub("[.]","",as.character(rho)))
      result = readRDS(paste0(base_results, "/", type.x, "/", type.problem, "/", filename, ".rds"))
      result_target = lapply(methods, function(xx){result[[xx]]})
      tmp_function <- function(xx){
        xx$loglossF = log( xx$lossF )
        xx$loglossR = log( xx$lossR )
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
          p1 = p1 + geom_signif( comparisons = list(c('"RAICc"', '"AICc"')), test = "wilcox.test", map_signif_level = FALSE,
                                 test.args=list(alternative = "two.sided", paired=TRUE, exact=FALSE), step_increase = .1 )
          y_mean = min(df_toplot[,metrics[k]]) - sd(df_toplot[,metrics[k]])
        }else{
          p1 = p1 + geom_jitter(size = 1.2)
          y_mean = -1
          p1 = p1 + scale_y_continuous(
            breaks = trunc( seq(0, max(df_toplot$nres)) ) )
        }
        # p1 = p1 + geom_signif( comparisons = list(c("RAICc", "AICc"), c("RAICc", "RCp"), c("RAICc", "LOO CV"), c("RAICc", "GCV")), test = "wilcox.test", map_signif_level = FALSE, 
        #                        test.args=list(alternative = "two.sided", paired=TRUE, exact=FALSE), step_increase = .1 )
        # p1 = p1 + coord_cartesian(ylim = ylim1)
        
        p1 = p1 + geom_text(data = means, aes_string(label = metrics[k], y = y_mean), size = 5)
        p1 = p1 + scale_x_discrete(name = "",
                                   breaks = methods_name,
                                   labels = parse(text = methods_name)
                                   #labels = expression("RAICc", "AICc", "RC"[p], "C"[p], "10FCV", "LOOCV", "S"[p], tilde(C)[p], "GCV", "BIC")
        )
        p1 = p1 + xlab("") + ylab(metrics_names[k])
        p1 = p1 + theme(legend.position = "none")
        p1 = p1 + ggtitle(paste0("n=", n[i], ", ", names(snr_name)[j]))
        p1 = p1 + theme(strip.text = element_text(size=20)) + theme(plot.title = element_text(size = 30, face = "bold", hjust=0.5), axis.text=element_text(size=15), axis.title=element_text(size=15))
        pp[[paste0("n=", n[i], ", ", names(snr_name)[j], ", ", metrics_brief[k])]] = p1
      }
    }
  }
  ## print to a file
  setEPS()
  postscript(file=paste0(base_plots, "/main/", type.x, "/", type.problem, "/Ex1.eps"), height=24, width=24)
  pp0 = grid.arrange(arrangeGrob(pp$`n=10, High signal, loss`,
                                 pp$`n=10, High signal, KL`,
                                 pp$`n=10, High signal, nres`,
                                 pp$`n=500, High signal, loss`,
                                 pp$`n=500, High signal, KL`,
                                 pp$`n=500, High signal, nres`,
                                 pp$`n=10, Low signal, loss`,
                                 pp$`n=10, Low signal, KL`,
                                 pp$`n=10, Low signal, nres`,
                                 pp$`n=500, Low signal, loss`,
                                 pp$`n=500, Low signal, KL`,
                                 pp$`n=500, Low signal, nres`,
                                 ncol=3))
  print(pp0)
  dev.off()
}
plot.generalrestriction.maintext <- function(){
  plot.generalrestriction.allsnr("fixedx")
  plot.generalrestriction.allsnr("randomx")
}
plot.generalrestriction.maintext()

## Boxplot in the supplemental material, for subset selection problem
plot.subsetselection.allmetrics.allp <- function( type.x, type, n, snr_name, rho ){
  ## parameters
  p = c(12, n/2, n-4)
  p_name = c("small", "medium", "large")
  type.problem = "subset_selection"
  nrep = 1000
  
  ## evaluation metrics
  metrics = c("loglossF", "loglossR", "logKLF", "logKLR", "nvar")
  metrics_names = c("Loss for fixed-X in log scale (loglossF)", "Loss for random-X in log scale (loglossR)", 
                    "KL for fixed-X in log scale (logKLF)", "KL for random-X in log scale (logKLR)",
                    "Size of the Selected Subset")
  metrics_brief = c("lossF", "lossR", "KLF", "KLR", "nvar")
  ## selection rules
  methods = do.call(c, list(lapply(c("raicc", "aicc", "rcp", "cp"), function(xx){c("ic", xx)}),
                            lapply(c("tenfold", "loo"), function(xx){c("cv", xx)}),
                            lapply(c("sp", "cphat", "gcv"), function(xx){c("ic", xx)}),
                            lapply(c("bic"), function(xx){c("ic", xx)})))
  # methods_name = c("RAICc", "AICc", "RCp", "Cp", "10FCV", "LOOCV", "Sp", "Cptilde", "GCV", "BIC")
  methods_name = c('"RAICc"', '"AICc"', "RC[p]", "C[p]", '"10FCV"', '"LOOCV"', "S[p]", "tilde(C)[p]", '"GCV"', '"BIC"')
  groups = c(1, 1, 2, 2, 3, 3, 4, 4, 4, 5)
  
  ## create the data for the plot
  pp = list()
  for(i in 1:length(p)){
    ## simulation results
    filename = paste0(type, "_n", n, "_p", p[i], "_", snr_name, "_rho", gsub("[.]","",as.character(rho)))
    result = readRDS(paste0(base_results, "/", type.x, "/", type.problem, "/", filename, ".rds"))
    result_target = lapply(methods, function(xx){result[[xx]]})
    tmp_function <- function(xx){
      xx$loglossF = log( xx$lossF )
      xx$loglossR = log( xx$lossR )
      xx$logKLF = log( xx$KLF )
      xx$logKLR = log( xx$KLR )
      xx
    }
    result_target = lapply(result_target, tmp_function)
    ## create the data frame for plot
    df_toplot = data.frame( do.call(cbind, lapply(metrics, function(metric){unlist(lapply(result_target, function(xx){xx[[metric]]}))})) )
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
      if(metrics[k] != "nvar"){
        p1 = p1 + geom_boxplot() 
        p1 = p1 + geom_signif( comparisons = list(c('"RAICc"', '"AICc"')), test = "wilcox.test", map_signif_level = FALSE,
                               test.args=list(alternative = "two.sided", paired=TRUE, exact=FALSE), step_increase = .1 )
        y_mean = min(df_toplot[,metrics[k]]) - sd(df_toplot[,metrics[k]])
      }else{
        p1 = p1 + geom_jitter(size = 1.2)
        y_mean = -round(p[i]/20)
        p1 = p1 + scale_y_continuous(
          breaks = trunc( seq(0, max(df_toplot$nvar), length.out=5) ) )
      }
      # p1 = p1 + geom_signif( comparisons = list(c("RAICc", "AICc"), c("RAICc", "RCp"), c("RAICc", "LOO CV"), c("RAICc", "GCV")), test = "wilcox.test", map_signif_level = FALSE, 
      #                        test.args=list(alternative = "two.sided", paired=TRUE, exact=FALSE), step_increase = .1 )
      # p1 = p1 + coord_cartesian(ylim = ylim1)
      
      p1 = p1 + geom_text(data = means, aes_string(label = metrics[k], y = y_mean), size = 5)
      p1 = p1 + scale_x_discrete(name = "",
                                 breaks = methods_name,
                                 labels = parse(text = methods_name)
                                 #labels = expression("RAICc", "AICc", "RC"[p], "C"[p], "10FCV", "LOOCV", "S"[p], tilde(C)[p], "GCV", "BIC")
      )
      p1 = p1 + xlab("") + ylab(metrics_names[k])
      p1 = p1 + theme(legend.position = "none")
      if(k==1){
        p1 = p1 + ggtitle(paste0("p=", p[i]))
      }
      p1 = p1 + theme(strip.text = element_text(size=20)) + theme(plot.title = element_text(size = 30, face = "bold", hjust=0.5), axis.text=element_text(size=15), axis.title=element_text(size=15))
      pp[[paste0(p_name[i], "p, ", metrics_brief[k])]] = p1
    }
  }
  ## print to a file
  setEPS()
  postscript(file=paste0(base_plots, "/supplement/", type.x, "/", type.problem, "/", type, "_n", n, "_", snr_name, "_rho", gsub("[.]","",as.character(rho)), ".eps"), height=28, width=24)
  pp0 = grid.arrange(arrangeGrob(pp[["smallp, lossF"]],
                                 pp[["mediump, lossF"]],
                                 pp[["largep, lossF"]],
                                 pp[["smallp, lossR"]],
                                 pp[["mediump, lossR"]],
                                 pp[["largep, lossR"]],
                                 pp[["smallp, KLF"]],
                                 pp[["mediump, KLF"]],
                                 pp[["largep, KLF"]],
                                 pp[["smallp, KLR"]],
                                 pp[["mediump, KLR"]],
                                 pp[["largep, KLR"]],
                                 pp[["smallp, nvar"]],
                                 pp[["mediump, nvar"]],
                                 pp[["largep, nvar"]],
                                 ncol=3))
  print(pp0)
  dev.off()
}
plot.subsetselection.supplement <- function(){
  ## Make plots in the supplement for the subset selection problem
  typex_name = c("Fixed-X", "Random-X")
  names(typex_name) = c("fixedx", "randomx")
  snr_fullname = c("high signal", "medium signal", "low signal")
  names(snr_fullname) = c("hsnr", "msnr", "lsnr")
  texfile = file(paste0(base, "/paper/plot_supplement_subsetselection.tex"), "w")
  count = 0
  for(type.x in c("fixedx", "randomx")){
    for(type in c("Sparse-Ex1", "Sparse-Ex2", "Dense-Ex1")){
      for(n in c(40, 200, 1000)){
        for(snr_name in c("hsnr", "msnr", "lsnr")){
          for(rho in c(0, 0.5, 0.9)){
            count = count + 1
            if(count %% 10 == 0){
              print(paste0(count, "/162"))
            }
            plot.subsetselection.allmetrics.allp( type.x, type, n, snr_name, rho )
            write(c("\\begin{figure}[!ht]",
                    "\\centering",
                    paste0("\\includegraphics[width=\\textwidth]{figures/supplement/", type.x,"/subset_selection/", type, "_n", n, "_", snr_name, "_rho", gsub("[.]","",as.character(rho)), ".eps}"),
                    paste0("\\caption{Variable selection, ", typex_name[names(typex_name)==type.x], ", ", type, ", ", "$n=", n, "$, ", snr_fullname[names(snr_fullname)==snr_name], ", and $\\rho=", rho, "$.}"),
                    "\\end{figure}",
                    "\\clearpage"
            ), file = texfile, append=TRUE)
          }
        }
      }
    }
  }
  close(texfile)
}
plot.subsetselection.supplement()
## Boxplot in the supplemental material, for general restriction problem
plot.generalrestriction.allmetrics.allrho <- function( type.x, type, n, snr_name ){
  ## parameters
  p = 6
  rho = c(0, 0.5, 0.9)
  rho_name = c("small", "medium", "large")
  type.problem = "general_restriction"
  nrep = 1000
  
  ## evaluation metrics
  metrics = c("loglossF", "loglossR", "logKLF", "logKLR", "nres")
  metrics_names = c("Loss for fixed-X in log scale (loglossF)", "Loss for random-X in log scale (loglossR)", 
                    "KL for fixed-X in log scale (logKLF)", "KL for random-X in log scale (logKLR)",
                    "Number of restrictions for the selected model")
  metrics_brief = c("lossF", "lossR", "KLF", "KLR", "nres")
  ## selection rules
  methods = do.call(c, list(lapply(c("raicc", "aicc", "rcp", "cp"), function(xx){c("ic", xx)}),
                            lapply(c("tenfold", "loo"), function(xx){c("cv", xx)}),
                            lapply(c("sp", "cphat", "gcv"), function(xx){c("ic", xx)}),
                            lapply(c("bic"), function(xx){c("ic", xx)})))
  # methods_name = c("RAICc", "AICc", "RCp", "Cp", "10FCV", "LOOCV", "Sp", "Cptilde", "GCV", "BIC")
  methods_name = c('"RAICc"', '"AICc"', "RC[p]", "C[p]", '"10FCV"', '"LOOCV"', "S[p]", "tilde(C)[p]", '"GCV"', '"BIC"')
  groups = c(1, 1, 2, 2, 3, 3, 4, 4, 4, 5)
  
  ## create the data for the plot
  pp = list()
  for(i in 1:length(rho)){
    ## simulation results
    filename = paste0(type, "_n", n, "_p", p, "_", snr_name, "_rho", gsub("[.]","",as.character(rho[i])))
    result = readRDS(paste0(base_results, "/", type.x, "/", type.problem, "/", filename, ".rds"))
    result_target = lapply(methods, function(xx){result[[xx]]})
    tmp_function <- function(xx){
      xx$loglossF = log( xx$lossF )
      xx$loglossR = log( xx$lossR )
      xx$logKLF = log( xx$KLF )
      xx$logKLR = log( xx$KLR )
      xx
    }
    result_target = lapply(result_target, tmp_function)
    ## create the data frame for plot
    df_toplot = data.frame( do.call(cbind, lapply(metrics, function(metric){unlist(lapply(result_target, function(xx){xx[[metric]]}))})) )
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
      if(metrics[k] != "nres"){
        p1 = p1 + geom_boxplot() 
        p1 = p1 + geom_signif( comparisons = list(c('"RAICc"', '"AICc"')), test = "wilcox.test", map_signif_level = FALSE,
                               test.args=list(alternative = "two.sided", paired=TRUE, exact=FALSE), step_increase = .1 )
        y_mean = min(df_toplot[,metrics[k]]) - sd(df_toplot[,metrics[k]])
      }else{
        p1 = p1 + geom_jitter(size = 1.2)
        y_mean = -1
        p1 = p1 + scale_y_continuous(
           breaks = trunc( seq(0, max(df_toplot$nres)) ) )
      }
      # p1 = p1 + geom_signif( comparisons = list(c("RAICc", "AICc"), c("RAICc", "RCp"), c("RAICc", "LOO CV"), c("RAICc", "GCV")), test = "wilcox.test", map_signif_level = FALSE, 
      #                        test.args=list(alternative = "two.sided", paired=TRUE, exact=FALSE), step_increase = .1 )
      # p1 = p1 + coord_cartesian(ylim = ylim1)
      
      p1 = p1 + geom_text(data = means, aes_string(label = metrics[k], y = y_mean), size = 5)
      p1 = p1 + scale_x_discrete(name = "",
                                 breaks = methods_name,
                                 labels = parse(text = methods_name)
                                 #labels = expression("RAICc", "AICc", "RC"[p], "C"[p], "10FCV", "LOOCV", "S"[p], tilde(C)[p], "GCV", "BIC")
      )
      p1 = p1 + xlab("") + ylab(metrics_names[k])
      p1 = p1 + theme(legend.position = "none")
      if(k==1){
        p1 = p1 + ggtitle(parse(text = paste0("rho==", rho[i])))
      }
      p1 = p1 + theme(strip.text = element_text(size=20)) + theme(plot.title = element_text(size = 30, face = "bold", hjust=0.5), axis.text=element_text(size=15), axis.title=element_text(size=15))
      pp[[paste0(rho_name[i], "rho, ", metrics_brief[k])]] = p1
    }
  }
  ## print to a file
  setEPS()
  postscript(file=paste0(base_plots, "/supplement/", type.x, "/", type.problem, "/", type, "_n", n, "_", snr_name, ".eps"), height=28, width=24)
  pp0 = grid.arrange(arrangeGrob(pp[["smallrho, lossF"]],
                                 pp[["mediumrho, lossF"]],
                                 pp[["largerho, lossF"]],
                                 pp[["smallrho, lossR"]],
                                 pp[["mediumrho, lossR"]],
                                 pp[["largerho, lossR"]],
                                 pp[["smallrho, KLF"]],
                                 pp[["mediumrho, KLF"]],
                                 pp[["largerho, KLF"]],
                                 pp[["smallrho, KLR"]],
                                 pp[["mediumrho, KLR"]],
                                 pp[["largerho, KLR"]],
                                 pp[["smallrho, nres"]],
                                 pp[["mediumrho, nres"]],
                                 pp[["largerho, nres"]],
                                 ncol=3))
  print(pp0)
  dev.off()
}
plot.generalrestriction.supplement <- function(){
  ## Make plots in the supplement for the general restriction problem
  typex_name = c("Fixed-X", "Random-X")
  names(typex_name) = c("fixedx", "randomx")
  snr_fullname = c("high signal", "medium signal", "low signal")
  names(snr_fullname) = c("hsnr", "msnr", "lsnr")
  m0 = c(4, 2, 3)
  names(m0) = paste0("Ex", seq(1,3))
  texfile = file(paste0(base, "/paper/plot_supplement_generalrestriction.tex"), "w")
  count = 0
  for(type.x in c("fixedx", "randomx")){
    for(type in c("Ex1", "Ex2", "Ex3")){
      for(n in c(10, 100, 500)){
        for(snr_name in c("hsnr", "msnr", "lsnr")){
          count = count + 1
          if(count %% 10 == 0){
            print(paste0(count, "/54"))
          }
          plot.generalrestriction.allmetrics.allrho( type.x, type, n, snr_name )
          write(c("\\begin{figure}[!ht]",
                  "\\centering",
                  paste0("\\includegraphics[width=\\textwidth]{figures/supplement/", type.x,"/general_restriction/", type, "_n", n, "_", snr_name, ".eps}"),
                  paste0("\\caption{General restriction, ", typex_name[names(typex_name)==type.x], ", $n=", n, "$, and ", snr_fullname[names(snr_fullname)==snr_name], ".}"),
                  "\\end{figure}",
                  "\\clearpage"
          ), file = texfile, append=TRUE)
        }
      }
    }
  }
  close(texfile)
}
plot.generalrestriction.supplement()

old <- function(){
  ## Plot an example where BIC overfits
  plot.bic.overfitting <- function(){
    n = c(200, 2000)
    pnratio = 0.98
    rho = 0
    snr = 8.5
    names(snr) = "hsnr"
    type = "Sparse-Ex1"
    
    df_toplot = list()
    for(i in 1:length(n)){
      p = round( n[i]*pnratio )
      filename = paste0(type, "_n", n[i], "_p", p, "_", names(snr), "_rho", gsub("[.]","",as.character(rho)))
      result = readRDS(paste0(base_results, "/", filename, ".rds"))
      
      rep = 10
      beta_covmatrix = gen.beta.covmatrix(p, rho, type)
      data = gen.data(n[i], snr, type, beta_covmatrix$beta, beta_covmatrix$beta0, beta_covmatrix$covmatrix, seed=rep)
      
      betahat = coef.nest(data$x, data$y)
      muhat = cbind(1, data$x) %*% betahat
      df = 1:(p+1)
      df_toplot[[i]] = data.frame( val = c(unlist( lapply(calc.ic.all(muhat, data$y, df, data$sigma)[c("aic", "aicc", "bic")], as.numeric) ),
                                           c(2 * (df + 1) + n[i], n[i] * (n[i] + df)/(n[i] - df - 2), log(n[i]) * (df + 1) + n[i])),
                                   k = rep(0:p, 6),
                                   ic.type = factor( rep(rep(c("AIC", "AICc", "BIC"), each = p+1), 2) ),
                                   type = factor( rep(paste0(c("IC, n=", "Penalty, n="), n[i]), each = 3*(p+1)) )
      )
    }
    df_toplot = do.call(rbind, df_toplot)
    
    p1 = ggplot() + geom_point(data=df_toplot[df_toplot$ic.type%in%c("AIC", "BIC"),], aes(y = val, x = k, colour = ic.type, shape = ic.type), stat="identity",size=2, stroke=1)
    p1 = p1 + scale_colour_manual(name  = "",
                                  breaks=c("AIC", "BIC"),
                                  labels=c("AIC", "BIC"),
                                  values=c("#F8766D", "#00BA38")) +
      scale_shape_manual(name  = "",
                         breaks=c("AIC", "BIC"),
                         labels=c("AIC", "BIC"),
                         values=c(5, 3))
    p1 = p1 + facet_wrap(.~ type, ncol=2, scales="free")
    p1 = p1 + theme(legend.text=element_text(size=20))
    p1 = p1 + theme(plot.title = element_text(size = 20, face = "bold", hjust=0.5), axis.text=element_text(size=10), axis.title=element_text(size=10,face="bold"))
    print(p1)
    # setEPS()
    # postscript(file=paste(base_plots, '/numvar_bs_lbs.eps',sep=""), height=3, width=6)
    # print(p1)
    # dev.off()
  }
  
  ## Compare RAICc with AICc
  plot.raicc.aicc.rcp.cp <- function(){
    n = 50
    type = "Sparse-Ex3"
    p = 20
    rho = 0.5
    snr = c(8.5, 0.2)
    names(snr) = c("hsnr", "lsnr")
    nrep = 1000
    beta_covmatrix = gen.beta.covmatrix(p, rho, type)
    
    methods_toplot = list(c("raicc", "aicc"), c("rcp", "cp"))
    labels_toplot = list(c("RAICc", "AICc"), c("RCp", "Cp"))
    
    pp = list()
    for(k in 1:length(snr)){
      ic_all = list()
      # r_square = c()
      for(rep in 1:nrep){
        data = gen.data(n, snr[k], type, beta_covmatrix$beta, beta_covmatrix$beta0, beta_covmatrix$covmatrix, seed=rep)
        # r_square = c(r_square, summary(lm(as.numeric(data$y) ~ data$x[, which(data$beta!=0)]))$r.squared)
        
        betahat = coef.nest(data$x, data$y)
        muhat = cbind(1, data$x) %*% betahat
        # sigma_sq_hat = colSums((data$y - muhat[, ncol(muhat), drop=FALSE])^2) / (n - sum(betahat[, ncol(muhat)]!=0))
        # use true sigma^2 here for Cp and RCp
        sigma_sq_hat = data$sigma^2
        ic_all[[rep]] = calc.ic.all(muhat, data$y, df=colSums(betahat!=0), sigma.sq = sigma_sq_hat) 
      }
      for(ii in 1:length(methods_toplot)){
        val = do.call(c, lapply(methods_toplot[[ii]], function(method){ colMeans(do.call(rbind, lapply(ic_all, function(xx){xx[[method]]})))  }))
        df_toplot = data.frame( val = val,
                                ic.type = factor( rep(methods_toplot[[ii]], each=p+1) ),
                                k = rep(0:p, length(allmethods))
        )
        p1 <- ggplot() + geom_point(data=df_toplot, aes(y = val, x = k, colour = ic.type, shape = ic.type), stat="identity",size=1.5)
        p1 <- p1 + xlab("Subset size") + ylab("")
        p1 <- p1 + scale_colour_manual(name = "",
                                       breaks = methods_toplot[[ii]],
                                       labels = labels_toplot[[ii]],
                                       values = c("#F8766D", "#00BA38")) +
          scale_shape_manual(name = "",
                             breaks = methods_toplot[[ii]],
                             labels = labels_toplot[[ii]],
                             values = c(5, 3))
        
        p1 <- p1 + theme(legend.text=element_text(size=5),legend.position="bottom") + theme(legend.title=element_blank())
        p1 <- p1 + theme(plot.title = element_text(size = 10, face = "bold", hjust=0.5),axis.text=element_text(size=10), axis.title=element_text(size=10,face="bold"))
        p1 <- p1 + ggtitle(names(snr)[k])
        pp[[ paste(methods_toplot[[ii]][2], names(snr)[k], sep=", ") ]] = p1
      }
    }
    mylegend = g_legend(pp$`aicc, hsnr`)
    setEPS()
    postscript(file=paste(base_plots, '/raicc_aicc.eps',sep=""), height=3, width=6)
    pp0 <- grid.arrange(arrangeGrob(pp$`aicc, hsnr` + theme(legend.position="none"),
                                    pp$`aicc, lsnr` + theme(legend.position="none"),
                                    nrow=1,ncol=2,heights=3),
                        mylegend,nrow=2,heights=c(9,1))
    print(pp0)
    dev.off()
    
    mylegend = g_legend(pp$`cp, hsnr`)
    setEPS()
    postscript(file=paste(base_plots, '/rcp_cp.eps',sep=""), height=3, width=6)
    pp0 <- grid.arrange(arrangeGrob(pp$`cp, hsnr` + theme(legend.position="none"),
                                    pp$`cp, lsnr` + theme(legend.position="none"),
                                    nrow=1,ncol=2,heights=3),
                        mylegend,nrow=2,heights=c(9,1))
    print(pp0)
    dev.off()
  }
  plot.raicc.aicc.rcp.cp()
}


