# METADATA 
########################################################################################################
# Project: DISSINET / BolIncr / whole network approach / ERGM
# Related manuscript : Incriminations in the inquisition register of Bologna (1291–1310) 
# Authors of the related manuscript : David Zbíral; Katia Riccardo; Zoltán Brys; Tomáš Hampejs
#
# NETWORK DESCRIPTIVES and VIS V2.0
#
# R-Code: Zoltán Brys and David Zbíral 
#
# Description: this R code 
#         1 prepares the environment
#         2 reads and checks input tables
#         3 define two main functions (directed graph descriptives, flag measures, )
#         4 define the graph and calculate graph statistics (desc, indeg, outdeg, triad)
#         5 create a visualization
########################################################################################################


# CODING CONVENTION
########################################################################################################
#flagvable names
# at_ denotes attribute/feature/term description vectors
# fn_ denotes file names (fn_inp_ is input file names, fn_out_ is output file names)
# df_ denotes data frames
# cx_ denotes unique cycle var (x = x + 1)
########################################################################################################


# 1 ENVIRONMENT PACKAGES INPUT/OUTPUT FILENAMES AND MAIN FUNCTION
########################################################################################################
#environment
  rm(list = ls()) #deleting the memory
  if (as.numeric(gsub(".*:(\\s*)(\\d+)(\\s+)\\d+.*", "\\2", (system("free -m", intern = TRUE)[2])))<2048)
    stop("Memory is likely not enough for ERGM!") #checking free memory
  if (!("stats" %in% (.packages()) )) stop("R Environment is not fully loaded!")

#libraries
  library(netUtils)
  library(netseg)
  library(igraph)

#input filenames (fn_)
  fn_inp_incr_nodes <- paste0(getwd(), "/data/df_nodes.tsv")
  fn_inp_incr_edges <- paste0(getwd(), "/data/df_edges.tsv")
#environment prepared, filenames (varaibles starting with fn_) set.
########################################################################################################

  
# 2 LOAD AND CHECK INPUT DATA 
########################################################################################################
#reading
  df_incr_nodes <- read.delim(fn_inp_incr_nodes, sep="\t", header=TRUE, fileEncoding = "UTF-8")
  df_incr_edges <- read.delim(fn_inp_incr_edges, sep="\t", header=TRUE, fileEncoding = "UTF-8")

#check if the reading was OK
  if (!exists("df_incr_nodes")) stop("Input data table of nodes is not loaded!")
  if (!exists("df_incr_edges")) stop("Input data table of edges is not loaded!")

#check if both are a data frame
  if (!class(df_incr_nodes) == "data.frame") stop("Input data table of nodes is not a data frame!")
  if (!class(df_incr_edges) == "data.frame") stop("Input data table of edges is not a data frame!")

#check if there are loops 
  if (dim(subset(df_incr_edges, from == to))[1] != 0) stop("There are self-loops!")

#check if there are multiple edges
  if (!identical(unique(df_incr_edges), df_incr_edges)) stop("There are multiple edges!")

#check if there are nodes defined multiple times
  if (!identical(unique(df_incr_nodes$name), df_incr_nodes$name)) stop("Attributes of one node defined multiple times")

#brief check of NAs
  if (sum(is.na.data.frame(df_incr_edges))>0) stop ("Edge list contains NA(s)")
  if (sum(is.na.data.frame(df_incr_nodes))>0) stop ("Node list contains NA(s)")
#input data loaded and checked
########################################################################################################


# 2 NETWORK DESC FUNCTION
########################################################################################################
#function of directed graph, basic char of a graph object
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


# 4 DEFINE GRAPH AND CREATE TABLES
########################################################################################################
#define the graph
  g_binc <- graph_from_data_frame( d = df_incr_edges , 
                                   directed = TRUE ,
                                   vertices = df_incr_nodes)
  
#create table3
  table3 <- descriptives_graph(g_binc)

#indeg, S2 Table
  indeg <- table(igraph::degree(g_binc, mode = "in"))
  indeg <- as.data.frame(indeg)
  colnames(indeg) <- c("degree", "indeg_freq")
  
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
  

#triad census of the observed graph and median random graphs
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
#graph defined, table3, indeg, outdeg, triad_cens calculated
########################################################################################################

     
# 5 Figures
########################################################################################################
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
  node_colors <- ifelse(V(g_binc)$sex == "m", "blue", "orange")
  
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
  
  #write Figure1.tiff
  dev.off()

#Supporting Information Figure 4, S4 Fig
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
  
  #write TIFF
  dev.off()
########################################################################################################  
