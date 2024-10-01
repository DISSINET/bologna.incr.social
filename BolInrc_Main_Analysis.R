# METADATA 
########################################################################################################
# Project: MUNI/ DISSINET 
# Related manuscript : Incriminations in the inquisition register of Bologna (1291–1310) 
# Authors of the related manuscript : David Zbíral; Katia Riccardo; Tomáš Hampejs; Zoltan Brys
#
# R-Code: Zoltan Brys and David Zbíral 
#
# Sections:
#         1  setting up the environment
#         2  load and check input datasets
#         3  characteristics of the trial subjects (Table 2.)
#         4  Jaccards (S1 Table)
#         5  define graph 
#         6  plot graph (Fig 1.)
#         7  network descriptives (Table 3)
#         8  aux graph stats (S2 Table, S3 Table, S4Fig, S5 Table)
#         9  generating new nodal attribute (involved under Bonifacius of Ferrara or Florius Vicenza)
#         10 ERGM preparation
#         11 ERGM configruations
#         12 ERGM
#         13 sensitivity analysis
#         14 save results
########################################################################################################


# 1 SETTING UP THE ENVIRONMENT
########################################################################################################
#memory check
  rm(list = ls()) #deleting the memory
  if (as.numeric(gsub(".*:(\\s*)(\\d+)(\\s+)\\d+.*", "\\2", (system("free -m", intern = TRUE)[2])))<2048) 
    stop("Memory is likely not will be enough for running an ERGM!") #checking free memory

#R-check
  if (!("stats" %in% (.packages()) )) stop("R Environment is not fully loaded!") #check R environment

#libraries
  library("parallel")
  library("Matrix")

  library("igraph")
  library("netseg")
  library("netUtils")
  library("network") 

  library("ergm") 
  library("statnet.common")
  library("sna") 
  library("tergm")
  library("networkDynamic")
  library("statnet")
  library("intergraph")
  library("ergMargins")

#special functions
  source("BolIncr_Special_Functions.R")

#input filenames (fn_)
  fn_inp_incr_nodes <- paste0(getwd(), "/data/df_nodes.tsv")
  fn_inp_incr_edges <- paste0(getwd(), "/data/df_edges.tsv")
#environment prepared, filenames (varaibles starting with fn_) set. 
########################################################################################################

  
# 2 LOAD AND CHECK INPUT DATASETS 
########################################################################################################
#reading
  df_incr_nodes <- read.delim(fn_inp_incr_nodes, sep="\t", header=TRUE, fileEncoding = "UTF-8")
  df_incr_edges <- read.delim(fn_inp_incr_edges, sep="\t", header=TRUE, fileEncoding = "UTF-8")

#check if the reading was OK
  if (!exists("df_incr_nodes")) 
    stop("Input data table of nodes is not loaded!")
  if (!exists("df_incr_edges")) 
    stop("Input data table of edges is not loaded!")

#check if both are a data frame
  if (!class(df_incr_nodes) == "data.frame") 
    stop("Input data table of nodes is not a data frame!")
  if (!class(df_incr_edges) == "data.frame") 
    stop("Input data table of edges is not a data frame!")

#check if there are loops 
  if (dim(subset(df_incr_edges, from == to))[1] != 0) 
    stop("There are self-loops!")

#check if there are multiple edges
  if (!identical(unique(df_incr_edges), df_incr_edges)) 
    stop("There are multiple edges!")

#check if there are nodes defined multiple times
  if (!identical(unique(df_incr_nodes$name), df_incr_nodes$name)) 
    stop("Attributes of one node defined multiple times")

#brief check of NAs
  if (sum(is.na.data.frame(df_incr_edges))>0) stop ("Edge list contains NA(s)")
  if (sum(is.na.data.frame(df_incr_nodes))>0) stop ("Node list contains NA(s)")
#input data loaded and checked
########################################################################################################


# 3 CHARACTERISTICS OF TRIAL SUBJECTS (TABLE 2)
########################################################################################################  
#defining binary, categorical and chr variable name vectors
  at_node_var_chr <- c("name" ,  
                       "label")
  at_node_var_cat <- c("family_id")
  at_node_var_bin <- c(  "sex" , 
                         "churchperson", 
                         "middling" ,  
                         "cathar_aff", 
                         "apostle_aff", 
                         "other_heterodoxy_aff",
                         "non_id_aff",
                         "deponent"  , 
                         "redeponent" ,
                         "ever_summoned" , 
                         "ever_pledged" , 
                         "ever_incarcerated", 
                         "ever_tortured",                      
                         "inq_FV", 
                         "inq_GV", 
                         "inq_GP", 
                         "inq_BdF")

#generate descriptives
  Table2_binary <- descriptives_df(data_frame = df_incr_nodes, 
                                    include = at_node_var_bin , 
                                    binaryvarmax = 30,
                                    pv = "full_dataset_binary") 
  
  Table2_family<- descriptives_df(data_frame = df_incr_nodes, 
                                       include = "family_id" , 
                                       binaryvarmax = -1,
                                       pv = "full_dataset_familiy_id")  
  
  Table2_nokinship =  sum(Table2_family$frequency==1) # No other kinship group member in the data
  Table2_kinship    = sum(Table2_family$frequency[Table2_family$frequency>1]) #At least one
#characteristics of trial subjects generated
########################################################################################################  

  
# 4 JACCARDS (S1 Table.)
########################################################################################################  
# calcualte Jaccards for all pairs after binary variance based selection
  S1_Table <- jaccard_matrixc(data_frame = df_incr_nodes, 
                                          include = at_node_var_bin) 
########################################################################################################  


# 5 DEFINE GRAPH 
########################################################################################################
#define the graph
  g_binc <- graph_from_data_frame( d = df_incr_edges , 
                                   directed = TRUE ,
                                   vertices = df_incr_nodes)
########################################################################################################

 
# 6 PLOT GRAPH (Fig.1)
########################################################################################################
#define the graph
  #Figure 1 - network vis
  #define TIFF
  tiff(filename = "Fig1.tiff",
       width = 33, height = 33, units = "cm", 
       compression = "lzw",
       bg = "white",
       res = 600
  )
  
  # Fruchterman-Reingold layout
  layout1 <- layout.fruchterman.reingold(g_binc)
  
  # color vector based on the "sex" attribute
  node_colors <- ifelse(V(g_binc)$sex == "1", "blue", "orange")
  
  # shapes based on the "deponent" attribute
  node_shapes <- ifelse((V(g_binc)$deponent == 1), "square", "circle")
  
  # set node sizes proportional to indegree
  node_sizes <- log(igraph::degree(g_binc, mode = "in")+3)
  
  #plot Fig1
  plot(
    g_binc, 
    layout = layout1, 
    vertex.label = NA, 
    vertex.color = node_colors, 
    vertex.shape = node_shapes, 
    vertex.size = node_sizes,
    edge.arrow.size = 0.3
  )
  
  #write Figure1.tiff (Fig 1.)
  dev.off()
#graph plotted
########################################################################################################

  
# 7 NETWORK DESCRIPTIVE (Table 3.)
########################################################################################################
#Table 3. Characteristics of the directed incrimination network.
  Table3 <- descriptives_graph(g_binc)
########################################################################################################

    
# 8 AUX GRAPH STATS (S2 Table, S3 Table, S4Fig, S5Table)
########################################################################################################
#S2 Table
#Indegree distribution of the incrimination network.
  S2_Table <- table(igraph::degree(g_binc, mode = "in"))
  S2_Table <- as.data.frame(S2_Table)
  colnames(S2_Table) <- c("degree", "indeg_freq")

  
#S3 Table
#Outdegree distribution of the incrimination network. “Involved under” characteristics are not disjunct.
#outdeg
  outdeg <-  as.data.frame(table(igraph::degree(g_binc, mode = "out")))
  rownames(outdeg) <- outdeg$Var1
  outdeg <- as.data.frame(outdeg)
  colnames(outdeg) <- c("degree","all_out")
  
#outdeg by FV
  outdeg_FV <- as.data.frame(table(igraph::degree(g_binc, mode = "out"), V(g_binc)$inq_FV)[,2])
  colnames(outdeg_FV) <- c("FV_out")
  
#outdeg by GV
  outdeg_GV <- as.data.frame(table(igraph::degree(g_binc, mode = "out"), V(g_binc)$inq_GV)[,2])
  colnames(outdeg_GV) <- c("GV_out")
  
#outdeg by GP
  outdeg_GP <- as.data.frame(table(igraph::degree(g_binc, mode = "out"), V(g_binc)$inq_GP)[,2])
  colnames(outdeg_GP) <- c("GP_out")
  
#outdeg by BdF
  outdeg_BdF <- as.data.frame(table(igraph::degree(g_binc, mode = "out"), V(g_binc)$inq_BdF)[,2])
  colnames(outdeg_BdF) <- c("BdF_out")
  
#outdeg all and by inquisitors, S3 Table
  outdegs <- cbind(outdeg, outdeg_FV, outdeg_GV, outdeg_GP, outdeg_BdF)  
  S3_Table <- outdegs

#S4 Fig
#Visualization of outdegree distribution by “Involved under” variables. Axes are logarithmic.
#define TIFF
  tiff(filename = "S4_Fig.tiff",
       width = 33, height = 33, units = "cm", 
       compression = "lzw",
       bg = "white",
       res = 600,
       pointsize = 24
  )

  outdegs$log_degree <- log(as.numeric(as.character(outdegs$degree)) + 1)
  outdegs$log_all_out <- log(outdegs$all_out + 1)
  outdegs$log_FV_out <- log(outdegs$FV_out + 1)
  outdegs$log_GV_out <- log(outdegs$GV_out + 1)
  outdegs$log_GP_out <- log(outdegs$GP_out + 1)
  outdegs$log_BdF_out <- log(outdegs$BdF_out + 1)

#plot S4 Fig.
  plot(outdegs$log_degree, outdegs$log_all_out, type = "l", col = "blue", xlab = "log(outdegree+1)", ylab = "log(freqency+1)", lwd = 3, cex = 1)
  lines(outdegs$log_degree, outdegs$log_FV_out, col = "red", lwd = 2)
  lines(outdegs$log_degree, outdegs$log_GV_out, col = "orange", lwd = 2)
  lines(outdegs$log_degree, outdegs$log_GP_out, col = "green", lwd = 2)
  lines(outdegs$log_degree, outdegs$log_BdF_out, col = "purple", lwd = 2)
  legend("topright", legend = c("All outdegree", 
                                "Florius Vicenza", 
                                "Guido Vicentinus", 
                                "Guido Parmensis", 
                                "Bonifacius de Feraria"), col = c("blue", 
                                                                  "red", 
                                                                  "orange" , 
                                                                  "green", 
                                                                  "purple"), lty = 1)
  
#write TIFF (S4 Fig.)
  dev.off()

#S5 Table
#Triad census of the observed graph and the median values of 10,000 random graphs of similar size.
  triad_cens <- NULL
  triad_cens <- as.data.frame(igraph::triad_census(g_binc))
  triad_nms <- c("003",
                 "012",
                 "102",
                 "021D",
                 "021U",
                 "021C",
                 "111D",
                 "111U",
                 "030T",
                 "030C",
                 "201",
                 "120D",
                 "120U",
                 "120C",
                 "210",
                 "300")
  
  rownames(triad_cens) <- triad_nms
  
  #median triad cencus of 10000 generated similar Erdos-Renyi graph
  rnd_triad_cens <- NULL
  for (c1 in 1:10000)
  {
    tmp_random_graph <- igraph::erdos.renyi.game( n= vcount(g_binc), p.or.m = ecount(g_binc), type="gnm", directed = TRUE ) 
    tmp_triad_cens <- as.data.frame(triad_census(tmp_random_graph))
    rnd_triad_cens <- rbind(rnd_triad_cens, t(tmp_triad_cens))
  }
  
  colnames(rnd_triad_cens) <- triad_nms 
  triad_medians <- as.data.frame(apply(rnd_triad_cens, 2, median))
  rownames(triad_medians) <- triad_nms
  
  #adding the results of 1000 random graph triad census to observed graph triad census, S4 Table
  triad_cens <- cbind(triad_cens, triad_medians)
  colnames(triad_cens) <- c("observed", "random")
  triad_cens$rat <- triad_cens$observed / triad_cens$random
  
  triad_cens$triadc_id <- rownames(triad_cens)
  rownames(triad_cens) <- c(1:dim(triad_cens)[1])
  triad_cens$rat[is.infinite(triad_cens$rat)] <- -1  
  
  S5_Table = triad_cens
#S2 Table, S3 Table, S4Fig, S5Table generated
########################################################################################################

  
# 9 GENERATING NEW NODAL ATTRIBUTE (involved under Bonifacius of Ferrara or Florius Vicenza)
########################################################################################################
#Visual analysis indicated a difference between the outdegree distribution 
# of trial subjects involved under the inquisitor Bonifacius of Ferrara or Florius Vicenza 
# and the outdegree distribution of trial subjects not involved under. 
  
#-> We add a generate a node property for the combined var
  V(g_binc)$inq_FV_or_inq_BdF <- df_incr_nodes$inq_BdF + df_incr_nodes$inq_FV
  
#-> We check Jaccard too
  df_incr_nodes$inq_FV_or_inq_BdF = df_incr_nodes$inq_BdF + df_incr_nodes$inq_FV
  S1_Table_extended = jaccard_matrixc(data_frame = df_incr_nodes, 
                                      include = c(at_node_var_bin, "inq_FV_or_inq_BdF"))
#fine  
########################################################################################################


# 10 ERGM PREPARATION
########################################################################################################
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
                       "non_id_aff"
                       )

  delete.vertex.attribute(net_incr, at_nodal_unused)

#check 
  Net_attrs = print.network(net_incr)
########################################################################################################

  
# 11 ERGM CONFIGURATIONS
########################################################################################################
#ergm samping space constraints
  maxdeg <- network.size(net_incr) - 1 #maxdeg
  deponent_boolean_matrix <- matrix(c(TRUE, FALSE, TRUE, FALSE), 
                                    nrow = 2, byrow = TRUE) #only deponents can accuse

  constraint_ergm <-  ( ~ bd(maxout = maxdeg , maxin = maxdeg ) + 
                          blocks(attr = ~deponent, levels2 = deponent_boolean_matrix ))

  control_ergm <- control.ergm(MCMC.maxedges = maxdeg^3,
                               parallel=2, parallel.type="PSOCK") #parallel proc

#null modells 
  null_mod_form <- formula(net_incr ~ edges)

#logical matrix for setting femmale->female, male->male as thetas of interest
  mm_boolean_matrix <- matrix(c(TRUE, FALSE, FALSE, TRUE), nrow=2, by=2)


#model
full_mod_form <- formula(net_incr ~ 
                           
                           #TOPOLOGICAL CONTROL      
                           edges + 
                           #outdegree in two subsets
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


# 12 ERGM
########################################################################################################

#network level statistic, specifically the number of edges meeting the ERGM terms condition (Table 4).
Table4_N = as.data.frame(summary(full_mod_form))

#null model
  ergm_full_null <- ergm(null_mod_form, 
                         constraints =  constraint_ergm,
                         control = control_ergm
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
  
#evaulation of ergm full
  par(mar = c(1, 1, 1, 1))
  mcmc_ergm_full <- ergm::mcmc.diagnostics(ergm_full, vars.per.page = 1) #S6 Document. MCMC
  gof_ergm_full <- ergm::gof(ergm_full)  # S7 Document. Goodness of Fit diagnostic of the main model.

  eval_full <- eval_ergm(ergm_full, 
                         VIFc = TRUE, 
                         MEc = TRUE, 
                         vp = "ergm_full_modell"
                         )
  Table4 = eval_full #plus NTable4 
########################################################################################################


# 13 SENSITIVITY ANALYSIS
########################################################################################################
#start parameters
  base_net_incr <- net_incr #save baseline net_incr
  
  at_n_edges <- sum(sna::degree(base_net_incr))/2  #number of edges baseline
  at_ten_p_edges <- round(0.1 * at_n_edges) #10% of the edges number
  val_eids <- network::valid.eids(base_net_incr) #valid edge ids
  
  df_res_sens <- NULL #results df of the sensitivity analysis

#repeated ERGM with 10% rand removed edges
  
#NOTE:
# random 10% removal of edges can remove too much edges for
#  F(nodefactor("cathar_aff")==1)~nodeifactor.middling.1:nodeofactor.deponent.1
#  F(nodefactor("apostle_aff")==1)~nodeifactor.middling.1:nodeofactor.deponent.1
# to be evaulated, this case evaulation of the model stop
# this case c1 to be adapted and continued
#
# also results can slightly differ in the paper due to random

for (c1 in 1:2)
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
}

#resulting tables
Table5_AIC = summary(df_res_sens$AIC)
Table5 <- eval_sens_res(df_res_sens, jp10mp90 = TRUE)
Table5 <- subset(df_raw_table5, 
                        select = c("ergm_term", "value", "p10", "med", "p90"))
########################################################################################################


# 14 SAVE RESULTS
########################################################################################################
#plots are saves as TIFFS.

#R-objects
save(Table2_binary,
     Table2_family,
     Table2_kinship,
     Table2_nokinship,
     
     Table3,
     
     Table4_N,
     Table4,
     
     Table5_AIC,
     Table5,
     
     S1_Table,
     S1_Table_extended,
     S2_Table,
     S3_Table,
     S5_Table,
     S6,
     S7,
     
     file="Bologna_Incriminations_Results.RData"
     )
########################################################################################################
