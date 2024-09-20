# METADATA 
########################################################################################################
# Project: DISSINET / BolIncr / whole network approach / ERGM
# Related manuscript : Incriminations in the inquisition register of Bologna (1291–1310) 
# Authors of the related manuscript : David Zbíral; Katia Riccardo; Tomáš Hampejs; Zoltan Brys
#
# DESCRIPTIVES V2.0
#
# CHANGES FROM V1.0 to V2.0 
# Based on reviewers' requests outlier-removal was deleted 
# and also (infered) non-identified religious affiliation is evaulated.
#
# R-Code: Zoltan Brys and David Zbíral 
#
# Description: this R code 
#         1 prepares the environment
#         2 reads and checks input tables
#         3 defines descriptive table function and Jaccard calculation function
#         4 generate Table 2. 
#         5 generate Supporting Information Table 1.
########################################################################################################


# CODING CONVENTION
########################################################################################################
#only use the absolutely neccessary packages
#variable naming:
  # at_ denotes attribute/feature/term description vectors
  # fn_ denotes file names (fn_inp_ is input file names, fn_out_ is output file names)
  # df_ denotes data frames
  # cx_ denotes unique cycle var (x = x + 1)
########################################################################################################


# 1 ENVIRONMENT PACKAGES INPUT/OUTPUT FILENAMES
########################################################################################################
#environment
  rm(list = ls()) #deleting the memory
  if (as.numeric(gsub(".*:(\\s*)(\\d+)(\\s+)\\d+.*", "\\2", (system("free -m", intern = TRUE)[2])))<4096) 
    stop("Memory is likely not will be enough for ERGM!") #checking free memory
  if (!("stats" %in% (.packages()) )) stop("R Environment is not fully loaded!")

#libraries
 #none

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

  
# 3 FUNCTIONS
########################################################################################################
#two functions are defined: descriptive tables and jaccard matrix
  
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
    if (sum(!(include %in% names(data_frame)))>0) stop ("One include variable name is not a column name")
    
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
      descriptive_table$binvar <- n1 * (descriptive_table$percentage/100) * (1 - (descriptive_table$percentage/100))
      descriptive_table$binvar_include <- (descriptive_table$binvar > binaryvarmax)
      }
    
    #adding pv descr.  
    descriptive_table$desc <-  rep(pv, dim(descriptive_table)[1] )
    
    return(descriptive_table)
  }
#end of decriptive table function

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

  
# 4 DESCRIPTIVE TABLE 2
########################################################################################################  
#defining binary, categorical and chr variable name vectors
at_node_var_chr <- c("name" ,  "label")
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
  table2_binary <- descriptives_df(data_frame = df_incr_nodes, 
                                    include = at_node_var_bin , 
                                    binaryvarmax = 30,
                                    pv = "full_dataset") 
  
  table2_family<- descriptives_df(data_frame = df_incr_nodes, 
                                       include = "family_id" , 
                                       binaryvarmax = -1,
                                       pv = "full_dataset") 

########################################################################################################  

  
# 6 DESCRIPTIVE AND JACCARD SIMILARITY MEASURE FOR BINARY VARIABLES
########################################################################################################  
# calcualte Jaccards for all pairs after binary variance based selection
  SP1_table_binary_jaccard <- jaccard_matrixc(data_frame = df_incr_nodes, 
                                          include = at_node_var_bin) 
########################################################################################################  
