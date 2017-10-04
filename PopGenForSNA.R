#  "PopGenForSNA_ESM.R"   (see "eff_deg_code.R" for precursor code)
#    Use four dashes to create folds in RStudio (small arrows to rt of line nums)
#  Script analyzes a .csv symmetric (undirected) WEIGHTED adjacency matrix ####
#              such as "AdjMat156.csv" (156 male LTM)
#        Row and column names are the ID names of the nodes (Long-tailed Manakins)
#        Weights represent the total number of interactions (2-hr Observation periods)
#            over a 17-yr period from 1983 to 1998, from the Formal Obs database
#
#   ***   VITAL FIRST STEP   ****
#   *** USE A PROJECT FOLDER WITH AT LEAST THE FOLLOWING FIVE SUBFOLDERS ***:
#      1) "Data"   2) "Scripts"  3) "Outputs"  4) "Figures and plots" 
#       *** if you don't, many functions will break ***  
#
#  Imports: {AdjacencyMatrix}.csv ## these are data files such as "AdjMat156.csv"
#
#  Depends: 
#      1) iGraph 1.0 R package: http://cran.r-project.org/web/packages/igraph/index.html
#      2) Matrix R package (eigenanalysis)
#      3) Moments R package (skew calculations)
#      4) tnet R package (Tore Opsahl's metrics for weighted networks)
# 
#  Where to find Eqns 1 to 5 of the manuscript:
#      Precursor edge-weight proportions
#         Function NodeMetricsFn in "EdgeWeightFunctions_ESM.R"
#      Eqn 1 (Observed edge-weight diversity)
#         Function NodeMetricsFn in "EdgeWeightFunctions_ESM.R"
#      Eqn 2 (Expected edge-weight diversity)
#         Function NodeMetricsFn in "EdgeWeightFunctions_ESM.R"
#      Eqn 3 (Effective degree)
#         Function NodeMetricsFn in "EdgeWeightFunctions_ESM.R"
#      Eqn 4 (Concentration)
#         Function NodeMetricsFn in "EdgeWeightFunctions_ESM.R"
#      Eqn 5 (F-statistics)
#         Function FxNStatsFn in "EdgeWeightFunctions_ESM.R"
#   created 8-Sep-2014 by D.B. McDonald
#
#   last modified 4-Nov-2017 by D.B. McDonald dbmcd@uwyo.edu 
#    Developed using RStudio.app 1.0.153 for Mac OSX Sierra 10.12.6  
#
#  If you use this R script (and the associated EdgeWeightFunctions_ESM.R script)
#   please cite 
#      McDonald, D.B., and E.A. Hobson. 2017. Edge-weight variance: population genetic metrics for social network analysis. 
#         Anim. Behav. (in press)   
#                
################################## End header do file input 

# Set wd() and file input ----------------------------------------------------

  ##   SET WORKING DIRECTORY  
  #     all  calls are relative to this directory, so to perform analyses simply
  #     change "wd" to point to YOUR project directory, then all the rest will work.
  #     On a Mac, you can just take out everything between the quotes, then drag
  #     the folder name into the quotes; the path will fill in automatically
  #
setwd("/Users/davidmcdonald/Desktop/October 4 versions/ESM R scripts and .csv files/")
  #  Check, if necessary, with getwd() 

  # load packages and source scripts
require(igraph)
require(Matrix)
require(moments)   ## For skewness() function
  ## Use some routines from Tore Opsahl's tnet R package in source -- "EdgeWeightFunctions.R" (set of functions)
require(tnet)
  
source("./Scripts/EdgeWeightFunctions.R")

  ## Import the adjacency matrix (perhaps best to have its individuals ordered in some logical way)
AdjFile <- file.choose()
AdjMat  <- read.csv(AdjFile,header=TRUE,row.names=1,check.names=FALSE)  # import adj mat str(AdjMat)
  ## extract the name of the input file (prior to outputting the .csv)
InputcsvFileName<-gsub(paste(getwd(),"/Data/", sep=""),"", AdjFile)
InputcsvFileName<-gsub(".csv","", InputcsvFileName)
  ######### End file input start initial background calcs 

# Initial calculations ----------------------------------------------------

AdjMat <- as.matrix(AdjMat) ## iGraph requires a matrix str(AdjMat) max(AdjMat)
MainGraph <- graph.adjacency(AdjMat,mode="undirected",weighted=T,diag=FALSE)
  ## MainGraph <- graph.adjacency(AdjMat,mode="undirected",weighted=NULL,diag=FALSE)
  ## Dimsn <- sqrt(length(AdjMat))
   ## BinaryAdjMat <- AdjMat
    ## for (i in 1:Dimsn){
    ## for (j in 1:Dimsn) {
    ## { ifelse(BinaryAdjMat[[i,j]]==0, 
    ##         BinaryAdjMat[[i,j]] <- BinaryAdjMat[[i,j]],
    ##       BinaryAdjMat[[i,j]] <- BinaryAdjMat[[i,j]]/BinaryAdjMat[[i,j]]) }
    ## }}   max(BinaryAdjMat)
    ## write.csv(BinaryAdjMat, "./Outputs/BinaryAdj.csv")

  ## CommntyDtctFstGrdyFn() is a function in "EdgeWeightFunctions.R" 
     ## It uses Newman's fast-greedy algorithm to assign the nodes to Communities
       ## Adds $Color and $FstGrdyCommnty to MainGraph,for Cols 9 & 10 of NodeStats (node metrics) output

FGDetectOutput <- CommntyDtctFstGrdyFn(MainGraph)  ## E(MainGraph)$weight
MainGraph <- FGDetectOutput$Graph                 ## MainGraph$Color has node colors;  V(MainGraph)$FstGrdyCommnty has community membership numbers
QScoreShrt <- round(FGDetectOutput$Qscr, digits=4)
print(paste("Modularity Q-score of the fast-greedy community detection function is ", QScoreShrt))
      ## E(MainGraph)$weight    head(E(MainGraph))    edge_attr(MainGraph)
  ## TEMP FIX FOR Rachel raccoons
    ## V(MainGraph)$FstGrdyCommnty[16] <- 4
    ## V(MainGraph)$FstGrdyCommnty[14] <- 4
    ## V(MainGraph)$FstGrdyCommnty[c(8:10,12,13,15)] <- 1
NoadLst <- unlist(V(MainGraph)$name)     ## Names/IDs (chr/string variable) of the nodes
EjCownt <- gsize(MainGraph)              ## Number of (binary) edges in undirected graph
EDjList <- as_edgelist(MainGraph)        ## head(EDjList)
EdgeIndex <- rep(0,EjCownt)              ## Used to color inter-community edges in plots (not implemented here)
for (i in 1:EjCownt){
  Ind1 <- which(NoadLst[] == EDjList[i,1])
  Ind2 <- which(NoadLst[] == EDjList[i,2])
  if(V(MainGraph)$FstGrdyCommnty[Ind1] != V(MainGraph)$FstGrdyCommnty[Ind2]) 
  {EdgeIndex[i] <- 1}
}   ##  sum(EdgeIndex) = 232 in LTM156 w/ fast-greedy = number community-bridging edges
PercentBridges <- sum(EdgeIndex)/EjCownt   ## str(EdgeIndex)
BridgeList <- which(EdgeIndex>0)
PrcntBrdgWt <- sum(E(MainGraph)$weight[BridgeList])/sum(E(MainGraph)$weight)
MeanBrdgWt <- sum(E(MainGraph)$weight[BridgeList])/sum(EdgeIndex)
StDevBrdgWt <- sd(E(MainGraph)$weight[BridgeList])
OnlyInCommWts <- E(MainGraph)$weight[-c(BridgeList)]   ## sum(OnlyInCommWts)
MeanInCommWt <- mean(OnlyInCommWts)
StDevInCommWts <- sd(OnlyInCommWts)

  ## YOU CAN SEE WHETHER your ModularityQ score is highest for fast-greedy using the lines below
   ## ModQ = 0.5855 for LTM data -- as good as any of the alternatives (Optimal intractable with this large a matrix)
  ## EdgBtwnComms <- cluster_edge_betweenness(MainGraph)
  ## modularity(EdgBtwnComms)    Only 0.4077 for LTM156 data 
   ## Walktrap <- cluster_walktrap(MainGraph)
   ## modularity(Walktrap)    0.5597 for LTM156 data 
      ## Spinglass <- cluster_spinglass(MainGraph)
      ## modularity(Spinglass)    0.1816 for LTM156 data 
        ## Optimal <- cluster_optimal(MainGraph)
        ## modularity(Optimal)    0.2060 for monk parakeet data (vs. 0.2012 for fast-greedy)
         ## InfoMap <- cluster_infomap(MainGraph)
         ## modularity(InfoMap)    0.5605 for LTM156 data 
           ## LabelProp <- cluster_label_prop(MainGraph)
           ## modularity(LabelProp)    0.5211 for LTM156 data
              ## LdEign <- cluster_leading_eigen(MainGraph)
              ## modularity(LdEign)    0.4762 for LTM156 data
        ## Louvain <- cluster_louvain(MainGraph)
        ## modularity(Louvain)    0.5853 for LTM156 data 

  ## MODIFIED SEP-2016 to incorporate communities with ≥ 3 nodes
BigCommNames <- FGDetectOutput$NameList           ## List of larger Communities for later use
NumBigComms <- length(BigCommNames)

  ## RandWeightGraphFn() is a function in "EdgeWeightFunctions.R"
    ## randomly redistributes weights of MainGraph over all edges
RandWtOutput <- RandWtGraphFn(MainGraph)
RandWtMain <- RandWtOutput$Graph

#################### End file input and initial calculations section 

################### Node-level analyses ######################################################### 
  ## Here's where Eqns 1 to 4 in Anim Behav paper are called
  ## NodeMetricsFn() is a function in "EdgeWeightFunctions.R" NODE-LEVEL METRICS
    ## Calculates $NodeID [$Degree] $EffDgr $DERatio $TotWt $EffWt $ObsEjWtDiv $ExpEjWtDiv  $Conc    
        # Cols    1       2         3       4        5      6      7           8            9
         ##       $FGCmnty $Color  $NrstNbrWtDgr $WtClstr [$WClsns] $Btwns  [$EignC]  $StClsns
         # Cols   10       11      12            13        14       15      16        17  
      ## and outputs the NodeStats dataframe to a .csv file
NdMtrOutput  <- NodeMetricsFn(MainGraph,InputcsvFileName) ## names(NdMtrOutput)
  ## NodeStats needed as argument for FxNStatsFn()
NodeStats    <- NdMtrOutput$NodeStatsFrame   ## head(NodeStats) head(NodeStats$Btwns)
V(MainGraph)$Btwns <- NodeStats$Btwns ## Added 8-Feb-17 to allow plotting for unweighted graphs (e.g., Edelman ERGM graphs)

RandName       <- paste("Rand", InputcsvFileName, sep="")
NdMtrOutput2   <- NodeMetricsFn(RandWtMain, RandName)
RandNodeStats  <- NdMtrOutput2$NodeStatsFrame

################### END  Node-level analyses 

############ Calculations at the Egonet level  ###########################
  ## Output Egonet stats as dataframe $EgoFrame and to .csv files w/ fn "EgoNetStatsFn"
      ## $EgoNdeCnt,EdgCnt,Dnsty,GlobClust,TotWt,WtMean,WtVar,WtSkew,ObsEdgWtDiv,ExpEdgWtDiv,EgoConc
OutputEgSts <- EgoNetStatsFn(MainGraph,InputcsvFileName)
    ## return is EgStsTable <- list(EgoFrame = EgostatsDF, FIE = Fie)
EgoStats    <- OutputEgSts$EgoFrame   ## head(EgoStats)
   ## ## FIE by EED as mean EEDego, weighted by Egonet total weight
     ## and OED over all individuals in all Egonets (thus, some individuals will apear many times)
        ## with weighting by their node edge-weight total (in the local context, which may sometimes involve pruning)
FIE         <- OutputEgSts$FIE  

OutputRandEgSts <- EgoNetStatsFn(RandWtMain,RandName)
RandEgoStats <- OutputRandEgSts$EgoFrame   ## head(RandEgoStats)
RandFIE <- OutputRandEgSts$FIE
################### END  Egonet-level analyses 

########################## Calculations at the community level     ########################## 
  ## Set up a dataframe for graph-level and community-level metrics such as density, edge-weight var and skew, and Kemeny's constant
GraphStats <- as.data.frame(matrix(ncol = 23, nrow = NumBigComms+2))
colnames(GraphStats) <- c("NdeCnt","EdgCnt","EffDgr","DEratio","Dnsty","GlobClust","Diam","TotWt","WtVar","WtSD","WtCoV","WtSkew","KemenyC","ModEignval2","ObsEdgWtDiv","ExpEdgWtDiv","Conc","FIE","FIC","FIN","FEN","FCN","ModularityQ")
rownames(GraphStats) <- c(InputcsvFileName,"Rand",BigCommNames)
GraphStats[] <- NA

  ## Output Community-level stats for .csv file GraphStats w/ fn "CommStatsFn()"
OutputCommSts <- CommStatsFn(MainGraph, GraphStats, NumBigComms)
   ## return is CommStsTable <- list(FIC = Fic, GrStats = GrphStats)
FIC <- OutputCommSts$FIC   ## FIC based on mean OED of each of the nodes in every community, weighted by their total edge-weight
GraphStats <- OutputCommSts$GrStats
################### END  Community-level analyses 

################### Full main graph-level analysis    ################################################################

  ## Fill GraphStats for MainGraph and RandWtMain graphs
   ## add 1st 17 weight-stats to the network-level dataframe "GraphStats using function "GraphStatsFn()"
GrStsOutput         <- GraphStatsFn(MainGraph)    ## Outputs vector $GrphStVctr
MainVector          <- GrStsOutput$GrphStVctr
GraphStats[1, 1:17] <- MainVector                 ## 1:17 because GraphStatsFn calculates 1st 17 cols
GraphStats[1, 23]   <- QScoreShrt

GrStsOutput2        <- GraphStatsFn(RandWtMain)
RandVector          <- GrStsOutput2$GrphStVctr
GraphStats[2, 1:17] <- RandVector                 ## 1:17 because GraphStatsFn calculates 1st 17 cols
GraphStats[1,"FIE"] <- FIE; GraphStats[2,"FIE"] <- RandFIE; GraphStats[1,"FIC"] <- FIC

  ## Note that I send NodeStats as an argument
FxNOutput           <- FxNStatsFn(MainGraph, NumBigComms)
GraphStats[1,"FIN"] <- FxNOutput$Fin; GraphStats[1,"FEN"] <- FxNOutput$Fen; GraphStats[1,"FCN"] <- FxNOutput$Fcn

FxNOutputRand       <- FxNStatsFn(RandWtMain, NumBigComms)
GraphStats[2,"FIN"] <- FxNOutputRand$Fin; GraphStats[2,"FEN"] <- FxNOutputRand$Fenrand

  ##### Output GraphStats
StatsName <- paste("./Outputs/",InputcsvFileName,"GraphStats.csv", sep="")
write.csv(GraphStats, StatsName, row.names=T)

  ## Test the fit of the negative binomial distribution to the edge weights, and output summary and plots
   ## [N.B.      Negative binomial analysis **** MAY FAIL FOR SOME INPUT MATRICES!!! **** ]
NegBinomFitFn(MainGraph, InputcsvFileName)

################# PLOTTING SCRIPT NEEDS WORK IF YOU want some provisional scripts, email Dave or Liz ###################### 

############## END SCRIPT ####