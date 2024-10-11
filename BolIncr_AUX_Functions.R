# METADATA 
########################################################################################################
# Project: MUNI/ DISSINET 
#
# Related manuscript : Incriminations in the inquisition register of Bologna (1291–1310) 
# Authors of the related manuscript : David Zbíral; Katia Riccardo; Tomáš Hampejs; Zoltan Brys
#
# R-code of auxiliary functions.
#
# R-code: Zoltan Brys and David Zbíral 
########################################################################################################


# FUNCTION: descriptive table
########################################################################################################
# this function create a descriptive table of categorical variables
# and calculates binary variance 
# and based on a treshold calculates to include/exclude that variable
# inputs:  
#   data_frame = the data frame
#   include =  the column names to analyze
#   binaryvarmax = decision threshold based on binary variables (max)
#   pv = passing descrp string variable
# if binaryvarmax < 0 , than function does not calculate it

descriptives_df <- function(data_frame, include = names(data_frame) , binaryvarmax = 12, pv = "") 
{
  #check the input
  if (!class(data_frame) == "data.frame") stop("Input data table is not a data frame!")
  if (is.null(include)) stop("No column names defined")
  if (sum(is.na.data.frame(data_frame))>0) stop("NAs in the data frame")
  if (sum(!(include %in% names(data_frame)))>0) stop("One include variable name is not a col name")
  
  #subset int columns and factorize them
  data_frame <- subset(data_frame, select = include)
  data_frame <- as.data.frame(lapply(data_frame, factor))
  
  #create an empty output table
  descriptive_table <- data.frame(
    variable = character(0),
    category = character(0),
    frequency = numeric(0),
    percentage = numeric(0),
    descr = character(0)
  )
  
  #cycle for each column  
  for (col in names(data_frame)) 
  {
    freq_table <- table(data_frame[[col]])
    percentages <- prop.table(freq_table) * 100
    
    variable_table <- data.frame(
      variable = col,
      category = names(freq_table),
      frequency = as.vector(freq_table),
      percentage = as.vector(percentages)
    )
    
    descriptive_table <- rbind(descriptive_table, variable_table)
  } #col cycle
  
  #adding variance and exluded
  #if binarymax is -1 or lower, then do not calculate
  if (binaryvarmax > -1) 
  {
    n1 <- dim(data_frame)[1]
    descriptive_table$binvar <- n1 * (descriptive_table$percentage/100) *
      (1 - (descriptive_table$percentage/100))
    descriptive_table$binvar_include <- (descriptive_table$binvar > binaryvarmax)
  }
  
  #adding pv descr.  
  descriptive_table$desc <-  rep(pv, dim(descriptive_table)[1] )
  
  return(descriptive_table)
}
#end of decriptive table function
########################################################################################################


# FUNCTION: Jaccard similarity calculation for all pairs
########################################################################################################
#Jaccard similarity coefficient calculation between two sets
jaccardc <- function (x, y) 
{
  if (anyNA(x)) stop("NAs in x")
  if (anyNA(y)) stop("NAs in y")
  if (length(x) != length(x)) stop("Lenghts are non-equal!")
  
  c11 = sum(x == 1 & y == 1)
  c10 = sum(x == 1 & y == 0)
  c01 = sum(x == 0 & y == 1)
  
  return (c11 / (c11 + c10 + c01))
}

#Jaccard similarity coefficient calculation between all pairs of a matrix
jaccard_matrixc <- function(data_frame, include = names(data_frame))
{
  #check the input
  if (!class(data_frame) == "data.frame") stop("Input data table is not a data frame!")
  if (is.null(include)) stop("No column names defined")
  if (sum(is.na.data.frame(data_frame))>0) stop("NAs in the data frame")
  if (sum(!(include %in% names(data_frame)))>0) stop("One include variable name is not a col name")
  
  #subset int columns
  data_frame <- subset(data_frame, select = include)
  
  #basic params
  no_vars <- dim(data_frame)[2] #number of variables
  at_names <- colnames(data_frame)
  
  #define result df
  df_res = as.data.frame(matrix(NA, nrow = length(at_names), ncol = length(at_names)))
  rownames(df_res) <- at_names
  colnames(df_res) <- at_names
  
  #calculate jaccard for every row/column 
  for (rows1 in 1:no_vars) for (cols1 in 1:no_vars) 
  {
    if (rows1 == cols1) {df_res[rows1,cols1] = -1} else 
    {df_res[rows1,cols1] = jaccardc(data_frame[,rows1], data_frame[,cols1])}
  }
  return(df_res)
}
########################################################################################################


# FUNCTION: network descriptive table 
########################################################################################################
#input: igraph object
#output: concaneted table with basic descriptives
descriptives_graph <- function(g_binc) 
{
  
  #check inputs
  if (!class(g_binc) == "igraph") stop("Input graph is not an igraph object!")
  
  #check graph
  if (any_loop(g_binc))     stop("Graph object has loops!")
  if (any_multiple(g_binc)) stop("Graph object multiple edges!")
  if (is_directed(g_binc)==FALSE) stop("Graph object is undirected!")
  
  #nodal parameteres
  no_nodes = igraph::vcount(g_binc)
  no_isolated_nodes = sum(igraph::degree(g_binc)==0)
  no_iso_c <- paste0(no_isolated_nodes," (", round((no_isolated_nodes/no_nodes), 4)*100 ,"%)")
  
  #edge parameters
  no_edges =  igraph::ecount(g_binc)
  mut_edges = sum(which_mutual(g_binc)) #number of mutual edges
  mut_edgesc <- paste0(mut_edges," (", round((mut_edges/no_edges), 4)*100 ,"%)")
  
  #degree parameters
  indeg =  igraph::degree(g_binc, mode="in")
  outdeg = igraph::degree(g_binc, mode="out")
  
  avg_indeg   = mean(indeg)
  med_indeg   = median(indeg)
  Q1_indeg    = quantile(indeg)[2]
  Q3_indeg    = quantile(indeg)[4]
  inmedQ1Q3   <- paste0(as.character(med_indeg), " [", as.character(Q1_indeg), "—"
                        , as.character(Q3_indeg), "]")
  
  avg_outdeg  = mean(outdeg)
  med_outdeg  = median(outdeg)
  Q1_outdeg   = quantile(outdeg)[2] 
  Q3_outdeg   = quantile(outdeg)[4]
  outmedQ1Q3  <- paste0(as.character(med_outdeg), " [", as.character(Q1_outdeg), "—"
                        , as.character(Q3_outdeg), "]")
  
  #structural  
  dens      = igraph::edge_density( g_binc, loops=FALSE)
  recip     = igraph::reciprocity(g_binc)
  recip_cor = netUtils::reciprocity_cor(g_binc)
  
  trans     = transitivity(as.undirected(g_binc))
  trans_cor = transitivity(as.undirected(g_binc))
  
  larg_com = igraph::decompose(g_binc, 
                               mode = "weak", 
                               min.vertices = max(igraph::components(g_binc)$csize) )[[1]]
  
  no_comp         = igraph::components(g_binc)$no
  diam_larg_comp = igraph::diameter(graph=larg_com)
  
  mean_pth_len = igraph::mean_distance(g_binc)
  
  #names
  table1_parameters <-
    c(
      "NODE CHARACTERISTICS",
      "Number of nodes \n (persons involved in the incrimination process)",
      "Number of isolated nodes",
      "EDGE CHARACTERISTICS",
      "Number of edges \n (representing incrimination of somebody else)",
      "Number of mutual edges",
      "DEGREE CHARACTERISTICS",
      "Average indegree",
      "Median [IQR] of indegrees",
      "Average outdegree",
      "Median [IQR] of outdegrees",
      "TOPOLOGICAL CHARACTERISTICS",
      "Density",
      "Reciprocity",
      "Reciprocity correlation coefficient",
      "Number of components",
      "Diamater of the largest component",
      "Mean path length"
    ) #end of table 1 parameter names
  
  
  #initiating the table
  table1_values <- 
    c(
      "",
      as.character(no_nodes),
      as.character(no_iso_c),
      "",
      as.character(no_edges) ,
      as.character(mut_edgesc ) ,
      "",
      as.character(round(avg_indeg,4)),
      as.character(inmedQ1Q3),
      as.character(round(avg_outdeg,4)),
      as.character(outmedQ1Q3),
      "",
      as.character(round(dens,4)), 
      as.character(round(recip,4)),
      as.character(round(recip_cor,4)),
      as.character(no_comp),
      as.character(diam_larg_comp),
      as.character(round(mean_pth_len,4))
    ) #end of table1 values
  
  #create table1
  table1 <- data.frame(table1_parameters=table1_parameters, table1_values=table1_values)
  colnames(table1) <- c("Parameters", "Values")
  
  return(table1)
}
########################################################################################################


#FUNCTIONS: custom summary which reports 10th percentila and 90th percentile instead of Q1, Q3
########################################################################################################
#input: series
#output: min, p10, mean, meadian, p90, max
csummary <- function(x) {
  c(
    min = min(x, na.rm = TRUE),
    p10 = quantile(x, 0.10, na.rm = TRUE),
    mean = mean(x, na.rm = TRUE),
    median = median(x, na.rm = TRUE),
    p90 = quantile(x, 0.90, na.rm = TRUE),
    max = max(x, na.rm = TRUE)
  )
}
########################################################################################################


#FUNCTION: ERGM-model (result object) evaulation
########################################################################################################
#input:
# ergm_mod = ergm model
# VIFc = calculation of VIF 
# MEc = calculation of Marginal Effects
# vp = passing variable
eval_ergm <- function(ergm_mod, VIFc = TRUE, MEc = TRUE, vp = "")
{
  #checking input
  if (!is(ergm_mod, "ergm")) {stop("ergm_mod is not an ergm object")}
  if (!is(VIFc, "logical")) {stop("VIFc is not an logical object")}
  if (!is(MEc, "logical")) {stop("MEc is not a logical object")}
  if (!is(vp, "character")) {stop("vp is not a string object")}
  
  #clean vp
  vp = as.character(vp)
  vp <- gsub(" ", "", vp)
  
  #basic summaries
  ergm_terms <- names(coef(ergm_mod))
  summary_ergm  <- summary(ergm_mod)
  coef_stat     <- coef(summary_ergm) 
  n_terms <- length(ergm_terms) #number of terms
  
  #performance measure
  AIC1    <- summary_ergm$aic[1]
  
  #vif
  vif_ergm <- as.data.frame(rep(NA, n_terms))  #vifs set to NA
  if (n_terms < 3 ) {VIFc <- FALSE} #no VIF can be calculated if there are less than two itmes
  if (VIFc) {vif_ergm <- vif.ergm(ergm_mod) #only if VIFc is true, calculate the VIF values
  vif_ergm <- t(data.frame(NA, vif_ergm))}
  rownames(vif_ergm)[1] <- ergm_terms[1]
  colnames(vif_ergm)[1] <- "vif"
  
  #marginal effects 
  ames  <- as.data.frame(rep(NA, n_terms)) #average marginal main effects set to NA
  amesp <- as.data.frame(rep(NA, n_terms)) #average marginal main effects set to NA
  
  if (n_terms < 2 ) {MEc <- FALSE}
  
  #cycle for MEC
  if (MEc) for (c1 in 1:n_terms)
  { ame_vs <- ergm.AME(ergm_mod,ergm_terms[c1])
  
  ames[c1,1]   <- ame_vs[1]
  amesp[c1,1]  <- ame_vs[4]
  }
  
  #setting colnames
  colnames(ames)[1]   <- "ames"
  colnames(amesp)[1]  <- "amesp"
  
  #adding to one table
  desc = rep(vp, n_terms)
  
  res_table <- data.frame(desc = desc,
                          terms = ergm_terms,
                          coef_stat, 
                          aic = as.data.frame(rep(AIC1, n_terms)), 
                          vif_ergm,
                          ames,
                          amesp
  )
  colnames(res_table) <- c("model_desc", "ergm_term", "B", "std", "MCMC", 
                           "z", "p", "AIC", 
                           "VIF", "AME", "AMEp")
  rownames(res_table) <- c(1:dim(res_table)[1])
  
  return(res_table)
}
#caution: this function uses ergMargins which uses bergm/tergm/ergMargins whose simply overwrite gof
########################################################################################################


#FUNCTIONS: Evaulation of the result of the sensitivity analysis
########################################################################################################
#df = input data frame of eval tables
#jp10mp90 = subsetting it to p10, median, p90
eval_sens_res <- function(res_df, jp10mp90 = TRUE)
{
  #checking input
  if (!is(res_df, "data.frame")) {stop("MEc is not a logical object")}

  #terms
  tmp_ergm_terms <- unique(res_df$ergm_term)
  
  sum_res <- NULL #result of sum
  for (i in tmp_ergm_terms)
  {
    sum_res <- rbind(sum_res, c(i, "Beta",csummary(res_df[res_df$ergm_term==as.character(i), ]$B)) )
    sum_res <- rbind(sum_res, c(i, "z", csummary(res_df[res_df$ergm_term==as.character(i), ]$z)) )
    sum_res <- rbind(sum_res, c(i, "p", csummary(res_df[res_df$ergm_term==as.character(i), ]$p))    )              
    sum_res <- rbind(sum_res, c(i, "AME", csummary(res_df[res_df$ergm_term==as.character(i), ]$AME)) )
  }
  
  sum_res <- as.data.frame(sum_res)
  colnames(sum_res) <- c("ergm_term", "value", "min", "p10", "mean", "med", "p90", "max")
  sum_res$min <- as.numeric(sum_res$min)
  sum_res$P10 <- as.numeric(sum_res$p10)
  sum_res$med <- as.numeric(sum_res$med)
  sum_res$mean <- as.numeric(sum_res$mean)
  sum_res$P90 <- as.numeric(sum_res$p90)
  sum_res$max <- as.numeric(sum_res$max)
  
  if (jp10mp90 == TRUE) {sum_res <- subset(sum_res, 
                               select = c("ergm_term", "value", "p10", "med", "p90"))}
  
  return(sum_res)
}
########################################################################################################
