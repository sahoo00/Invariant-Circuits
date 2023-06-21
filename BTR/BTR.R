
#-------------------------- 
rm(list=ls())
set.seed(0)  #use to ensure reproducibility. remove in actual use.
setwd("/BTR")
# (1) Setup paths and environment.
library(BTR)

# If intending to use parallel processing, uncomment the following lines.
library(doParallel) #num_core = 4 #specify the number of cores to be used.
# doParallel::registerDoParallel(cores=num_core)
#Setting up for parallel processing.
# num_cores = 10
# path="/BTR"
# setwd(path)
# registerDoParallel(cores=num_cores) #this automatically calls mclapply() when no cl is given.

# (2) Load data.
data = read.table("expr-btr.txt", header=TRUE, sep="\t")
print(dim(data))
data1 = data[,2:16]
rownames(data1) = data[,1]
print(dim(data1))
cdata = initialise_raw_data(data1, max_expr = "low")

fcdata = cdata


int_state <- read.csv(file = 'initial.csv', header = TRUE, row.names = 1)
int_M <- read.csv(file = 'initialB.csv', header = TRUE, row.names = 1)
bmodel = initialise_model(int_M)
colnames(fcdata) = tolower(colnames(fcdata))
colnames(int_state) = tolower(colnames(int_state))
# , self_loop = TRUE
final_model = model_train(cdata = fcdata , bmodel=bmodel, istate = int_state, max_varperrule = 10, self_loop = TRUE, verbose = T)



# (5) Inferring Boolean model.
saveRDS(final_model, "final_model_btr-T.rds")
outgraph_model(final_model)
writeBM(final_model, 'final_model_btr_rule-T.txt')



# 
# 
dysfunctional=c("ACTA2","TAGLN",'BMPR1A',"PDGFRB","MYL9","CSTK","WISP2",'ACRVL1')
funct=c("IHH","FN1","FERMT1",'SOX9',"PHOSPHO1",'ACRBP')
shared=c("EGFR",'FBXO7')



name2="BTR_result"
# install.packages("visNetwork")
require(visNetwork, quietly = TRUE)
library(data.table)
edges <- fread("edges.csv")
nodes <- fread("nodes.csv")

colnames(edges)


colnames(nodes)
names(edges)[names(edges) == "Source"] <- "from"
names(edges)[names(edges) == "Target"] <- "to"

list_str <- sort(unique(nodes$Id))
list_int <- 1:length(list_str)
names(list_int) = list_str
nodes$id = list_int[nodes$Id]

recode_str_int_nodes <- function(df) {
  df[, shape := ifelse(grepl('^and_', df$Id), "circularImage", "circle")]
  df[, image := ifelse(grepl('^and_', df$Id), "http://hegemon.ucsd.edu/~sataheri/KZC.svg", "")]
  df[, value := ifelse(grepl('^and_', df$Id), "10", "80")]
  df[, ty := ifelse(grepl('^and_', df$Id), "and", df$Id)]
  df$color.border = rep("black", times = nrow(df))
  df$borderWidth = rep("0.2", times = nrow(df))
  df$widthConstraint = rep("30", times = nrow(df))
  # df[, color.background := ifelse(df$ty == 'and', "white",'skyblue')]
  df[, label := ifelse(grepl('^and_',df$Id), '',  toupper(df$Id))]
  df[, shadow.size := ifelse(grepl('^and_', df$Id), "5", "15")]
  df[, group := ifelse(df$ty == 'and',"",ifelse(toupper(df$Id) %in% toupper(funct),
                       'Funct',ifelse(toupper(df$Id) %in% toupper(shared),'Shared', "Dysfunct")))]
  df[, color.background := ifelse(df$ty == 'and',"",ifelse(toupper(df$Id) %in% toupper(funct) ,
                                                'tomato',ifelse(toupper(df$Id) %in% toupper(shared),'green', "grey")))]


  df
#
} # recode_str_int and  set other features
#
nodes_m <- recode_str_int_nodes(nodes)
#
#
# #-------------
library(dplyr)
library(plyr)
library(stringr)
library(tidyr)
j = dplyr::rename(dplyr::rename(dplyr::rename(edges, "Id"="to") %>%inner_join(nodes, by = "Id"), "to"="id"),"to_name"="Id")
j = dplyr::rename(dplyr::rename(dplyr::rename(j, "Id"="from") %>%inner_join(nodes, by = "Id"), "from"="id"),"from_name"="Id")
# #-------------------
recode_str_int_edges <- function(df) {


  df$arrows  = rep("to", times = nrow(df))
  df[, color := ifelse(Directed == 'activates',
                       "darkblue",
                       ifelse(Directed == 'inhibits', "red", 'rest'))]
  # df$dashes = c(TRUE, FALSE)
  df[, dashes := ifelse(grepl('activates', df$Directed), FALSE, TRUE)]


  # df[, dashes:= ifelse(Directed=='activates',"FALSE",ifelse(Directed=='inhibits',"TRUE",'rest'))]

  df$shadow  = rep("TRUE", times = nrow(df))

#   df$arrows.to.type  = rep("arrow", times = nrow(df))
#   df$arrows.to.enabled  = rep("TRUE", times = nrow(df))




#

  # 3. Result
  df
} # recode_str_int

edges_m <- recode_str_int_edges(setDT(j))
colnames(edges_m)
edges_m=edges_m[,c("from","to","arrows","color","dashes")]
write.csv(edges_m, file =  paste0("edges_m", "_", name2, ".csv"))
write.csv(nodes_m, file = paste0("nodes_m", "_", name2, ".csv"))
#------------------------------

label_vec <- c(red = "activate", blue = "inhibits")


legendNodes <- data.frame(
  label = c( "Fucntional","shared","dysfunctional",'AND'),
  color.background = c( "#FF6347",'green','grey','white'),
  color.border = c("black", "black", "black","black"),
  shape = c( "dot","dot","dot","image"),
  image = c("","","" ,"http://hegemon.ucsd.edu/~sataheri/KZC.svg")
)
net <- visNetwork(nodes_m, edges_m, width = "100%") %>%
  visLegend(
    useGroups = FALSE,
    addNodes = legendNodes,

    addEdges = data.frame(
      label = c("activate", 'inhibit'),
      color = c("darkblue", 'red'),
      dashes = c(FALSE, TRUE)
    )
  ) %>%
    visOptions(collapse = TRUE)  %>%
    visOptions(highlightNearest = TRUE,selectedBy = "label") %>%
    visPhysics(solver = "forceAtlas2Based",
                forceAtlas2Based = list(gravitationalConstant = -28)) %>%
    visLayout(randomSeed = 114)
    

#visEdges(arrows = list(to = list(enabled = TRUE, type = "bar"))) %>% 

name2='visnet'
net %>% visSave(file =  paste0("BTR_NET", "_", name2, ".html"))
# rm(list=ls())







#
#pdf("result.pdf")
#
## Creating a plot
#plotBM(final_model)
#
## Closing the graphical device
#dev.off() 
