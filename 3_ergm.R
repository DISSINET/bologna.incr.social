# METADATA 
########################################################################################################
# Project: DISSINET / BolIncr / ERGM
# Related manuscript : Incriminations in the inquisition register of Bologna (1291–1310) 
# Authors of the related manuscript : David Zbíral; Katia Riccardo; Zoltán Brys; Tomáš Hampejs
#
# ERGMs - Exponential random graph models V2.0 serial incriminators included
#
# R-Code: Zoltán Brys and David Zbíral 
#
# Description: this R code 
#         1 prepares the environment
#         2 reads/checks input tables
#         3 define ergm-eval function
#         4 define network object
#         5 define ergm config (formula, constrains, control)
#         6 ergm full net
#         7 sensitivity ergms
########################################################################################################


# CODING CONVENTION
########################################################################################################
#variable naming convention
# at_ denotes attribute/feature/term description vectors
# fn_ denotes file names (fn_inp_ is input file names, fn_out_ is output file names)
# df_ denotes general data frames
# cx_ denotes unique cycle var (x = x + 1)
# ergm_ denotes ergm models
# eval_ denotes ergm model evaulation data frames
# inside cycles variable denotes are uncoventional
# temporary variable naming is also unconventional
########################################################################################################


# 1 ENVIRONMENT PACKAGES INPUT/OUTPUT FILENAMES AND MAIN FUNCTION
########################################################################################################
#environment
rm(list = ls())
if (as.numeric(gsub(".*:(\\s*)(\\d+)(\\s+)\\d+.*", "\\2", (system("free -m", intern = TRUE)[2])))<4192) 
  stop("Memory are likely not will be enough for ERGM!") #checking free memory
if (!("stats" %in% (.packages()) )) stop("R Environment is not fully loaded!")

#libraries
library("parallel")
library("Matrix")
library("igraph")
library("network") 
library("ergm") 
library("statnet.common")
library("sna") 
library("tergm")
library("networkDynamic")
library("statnet")
library("intergraph")
library("ergMargins")

#input filenames (fn_)
fn_inp_incr_nodes <- paste0(getwd(), "/data/df_nodes.tsv")
fn_inp_incr_edges <- paste0(getwd(), "/data/df_edges.tsv")
#environment prepared, filenames (varaibles starting with fn_) set. 
########################################################################################################


# 2 LOADING INPUT DATA 
########################################################################################################
#reading the two graphs, data.frame, df_
df_incr_nodes <- read.delim(fn_inp_incr_nodes, sep="\t", header=TRUE, fileEncoding = "UTF-8")
df_incr_edges <- read.delim(fn_inp_incr_edges, sep="\t", header=TRUE, fileEncoding = "UTF-8")

#quick check
if (!exists("df_incr_nodes")) stop("Input data table of nodes is not loaded!")
if (!exists("df_incr_edges")) stop("Input data table of edges is not loaded!")

if (!class(df_incr_nodes) == "data.frame") stop("Input data table of nodes is not a data frame!")
if (!class(df_incr_edges) == "data.frame") stop("Input data table of edges is not a data frame!")

#check if there are loops 
if (dim(subset(df_incr_edges, from == to))[1] != 0) stop("There are self-loops!")

#check if there are multiple edges
if (!identical(unique(df_incr_edges), df_incr_edges)) stop("There are multiple edges!")
#input data loaded and checked
########################################################################################################


# 3 FUNCTIONS: ERGM RESULTS EVAULATION
########################################################################################################
#general frequentist ergm evaulation function
#with optional calculation of VIF and Marginal Effects
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

#evaulation the sensitivity analysis data frame for Table 5
#df = input data frame of eval tables
#vp = variable passing
eval_sens_res <- function(res_df, vp = "")
{
  #checking input
  if (!is(res_df, "data.frame")) {stop("MEc is not a logical object")}
  if (!is(vp, "character")) {stop("vp is not a string object")}
  
  #terms
  tmp_ergm_terms <- unique(res_df$ergm_term)
  
  sum_res <- NULL #result of sum
  for (i in tmp_ergm_terms)
  {
    sum_res <- rbind(sum_res, c(i, "Beta",  summary(res_df[res_df$ergm_term==as.character(i), ]$B)) )
    sum_res <- rbind(sum_res, c(i, "z",  summary(res_df[res_df$ergm_term==as.character(i), ]$z)) )
    sum_res <- rbind(sum_res, c(i, "p",  summary(res_df[res_df$ergm_term==as.character(i), ]$p))    )              
    sum_res <- rbind(sum_res, c(i, "AME",  summary(res_df[res_df$ergm_term==as.character(i), ]$AME)) )
  }
  
  sum_res <- as.data.frame(sum_res)
  colnames(sum_res) <- c("ergm_term", "value", "min", "Q1", "med", "mean", "Q3", "max")
  sum_res$min <- as.numeric(sum_res$min)
  sum_res$Q1 <- as.numeric(sum_res$Q1)
  sum_res$med <- as.numeric(sum_res$med)
  sum_res$mean <- as.numeric(sum_res$mean)
  sum_res$Q3 <- as.numeric(sum_res$Q3)
  sum_res$max <- as.numeric(sum_res$max)
  
  return(sum_res)
}

#subfunction Jaccard calulactions
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

#jaccard similarity calculation for all pairs
jaccard_matrixc <- function(data_frame, include = names(data_frame))
{
  #check the input
  if (!class(data_frame) == "data.frame") stop("Input data table is not a data frame!")
  if (is.null(include)) stop("No column names defined")
  if (sum(is.na.data.frame(data_frame))>0) stop("NAs in the data frame")
  if (sum(!(include %in% names(data_frame)))>0) stop ("One include variable name is not a column name")
  
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


# 4 DEFINE NETWORK
########################################################################################################
#creating a binary variable for involved under Florius Vicenza or under Bonifacius of Ferrara
df_incr_nodes$inq_FV_or_inq_BdF = df_incr_nodes$inq_BdF + df_incr_nodes$inq_FV

#post analysis of JI of the new var
pji <- jaccard_matrixc(df_incr_nodes, include = c(  "sex" , 
                                             "churchperson", 
                                             "middling" ,  
                                             "cathar_aff", 
                                             "apostle_aff", 
                                             "redeponent" ,
                                             "ever_summoned" , 
                                             "ever_pledged" ,
                                             "inq_FV_or_inq_BdF")
                  )
#maximum ji of the involved bin attrs
max(pji)

#define the full graph
#previous analysis were done in igraph and to maintain consistency we first use igraph
g_binc <- graph_from_data_frame( d = df_incr_edges , 
                                 directed = TRUE ,
                                 vertices = df_incr_nodes)

#data transformation to network format
net_incr <- intergraph::asNetwork(g_binc)

#attributes not used in the ergm
at_nodal_unused <- c("name" ,  
                     "label" , 
                     "ever_incarcerated", 
                     "ever_tortured",
                     "inq_FV",
                     "inq_BdF",
                     "inq_GP",
                     "inq_GV",
                     "non_id_aff",
                     "other_heterodoxy_aff"
)

delete.vertex.attribute(net_incr, at_nodal_unused)

#check 
print.network(net_incr)
########################################################################################################


# 5 ERGM CONFIG
########################################################################################################
#ergm samping space constraints
maxdeg <- network.size(net_incr) - 1 #maxdeg
deponent_boolean_matrix <- matrix(c(TRUE, FALSE, TRUE, FALSE), nrow = 2, byrow = TRUE) #only deponents can accuse

constraint_ergm <-  ( ~ bd(maxout = maxdeg , maxin = maxdeg ) + 
                        blocks(attr = ~deponent, levels2 = deponent_boolean_matrix ))

control_ergm <- control.ergm(MCMC.maxedges = maxdeg^3, parallel=2, parallel.type="PSOCK") #parallel proc

#null modells 
null_mod_form <- formula(net_incr ~ edges)

#logical matrix for setting femmale->female, male->male as thetas of interest
mm_boolean_matrix <- matrix(c(TRUE, FALSE, FALSE, TRUE), nrow=2, by=2)


#model
full_mod_form <- formula(net_incr ~ 
                           
                           #TOPOLOGICAL CONTROL      
                           edges + 
                           #include and control serial incriminators - reviewer's 2 request
                           F(~gwodegree(decay = 0.7, fixed = TRUE), ~nodefactor("inq_FV_or_inq_BdF") == 0) +
                           F(~gwodegree(decay = 0.7, fixed = TRUE), ~nodefactor("inq_FV_or_inq_BdF") == 1) +
                           
                           #DYADIC CONTROL
                           mutual(by="deponent", levels=2) + 
                           (nodematch("cathar_aff") : nodeofactor("deponent"))+
                           (nodematch("apostle_aff") : nodeofactor("deponent"))+
                           
                           #DYADIC INPUT
                           (nodematch("family_id") : nodeofactor("deponent")) +
                           (nodemix("sex", levels2 = mm_boolean_matrix ) : nodeofactor("deponent")) +
                           
                           #NODAL CONTROL 
                           (nodeofactor("redeponent") : nodeofactor("deponent")) + 
                           (nodeofactor("ever_summoned") : nodeofactor("deponent")) +
                           (nodeofactor("ever_pledged") : nodeofactor("deponent")) + 
                           
                           #NODAL INPUT
                           (nodeofactor("churchperson") : nodeofactor("deponent")) +
                           F(~(nodeifactor("middling") : nodeofactor("deponent")), ~nodefactor("cathar_aff") == 1) +
                           F(~(nodeifactor("middling") : nodeofactor("deponent")), ~nodefactor("apostle_aff") == 1)
) 

#disable warning about ill-defines loglik due to sample constrains
options(ergm.loglik.warn_dyads=FALSE)
########################################################################################################


# 6 FULL ERGM
########################################################################################################

#network level statistic, specifically the number of edges meeting the ERGM terms condition (Table 4).
summary(full_mod_form)

#null model
ergm_full_null <- ergm(null_mod_form, 
                       constraints =  constraint_ergm,
                       control = control_ergm, 
)

eval_full_null <- eval_ergm(ergm_full_null, 
                            VIFc = FALSE, 
                            MEc = FALSE, 
                            vp = "full_null_model" 
)

#ergm model
ergm_full <- ergm(full_mod_form , 
                  constraints =  constraint_ergm,
                  control = control_ergm)

par(mar = c(1, 1, 1, 1))
mcmc_ergm_full <- ergm::mcmc.diagnostics(ergm_full, vars.per.page = 1)
gof_ergm_full <- ergm::gof(ergm_full)  

eval_full <- eval_ergm(ergm_full, 
                       VIFc = TRUE, 
                       MEc = TRUE, 
                       vp = "ergm_full_modell" 
)
########################################################################################################


# 7 SENSITIVITY ANALYSIS
########################################################################################################
#start parameters
base_net_incr <- net_incr #save baseline net_incr

at_n_edges <- sum(sna::degree(base_net_incr))/2  #number of edges baseline
at_ten_p_edges <- round(0.1 * at_n_edges) #10% of the edges number
val_eids <- network::valid.eids(base_net_incr) #valid edge ids

df_res_sens <- NULL #results df of the sensitivity analysis

#repeated ERGM with 10% rand removed edges
for (c1 in 1:100)
{
  #reseting baseline
  net_incr <- base_net_incr 
  
  #remove 10% of the edges
  smp_edges <- sample(val_eids, size = at_ten_p_edges) #10% random
  net_incr <- network::delete.edges(net_incr, eid = smp_edges) #deleting 10% of the edges
  
  #ergm null
  ergm_tmp <- ergm(full_mod_form,
                   constraints = constraint_ergm,
                   control = control_ergm
  )
  #eval
  eval_tmp <- eval_ergm(ergm_tmp, 
                        VIFc = TRUE, 
                        MEc = TRUE, 
                        vp = paste0(as.character(c1), "_", as.character(sum(sna::degree(net_incr))/2)) #cycle_edges
  )
  #save
  df_res_sens <- rbind(df_res_sens, eval_tmp) 
  Sys.sleep(1)
  cat(" \n----------------------------------", c1 , "-------------------------------- \n")
}

#resulting tables
summary(df_res_sens$AIC)
df_raw_table5 <- eval_sens_res(df_res_sens)
########################################################################################################
