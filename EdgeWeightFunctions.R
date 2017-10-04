############################################################################## 
#  "EdgeWeightFunctions.R"   (source for analyses in "PopGenForSNA.R") 
#     N.B. in RStudio functions fold automatically (small arrow to rt of line num)
#  ----
#  This script runs functions on weighted graphs to assess structure and function 
#             of weighted networks
#      Set up as functions because will be used repeatedly for inputs such as 
#      1) Original graph
#      2) Graph with weights randomly scattered over the same edge structure as original
#      3) Graph with weights reallocated to enhance differences among communities
#               detected by one or more community-detection algorithms
#
#     All but Eqn 5 are embedded in Function NodeMetricsFn
#
#      List of functions (with their arguments), where n is number of nodes in inGraph
#       Function name                      Output
#         RandWtGraphFn(inGraph)             graph with same edge-structure as original, but weights redistributed randomly
#         CommntyDtctFstGrdyFn(inGraph)      list of communities (subgraphs) in original graph (by Newman fast-greedy method)
#         OpsahlWtMetricsFn(inGraph)         Opsahl weighted closeness (col 14 in NodeStats output)
#         KemenyFn(inGraph)                  scalar of Kemeny's constant for inGraph
#         ModEignVal2Fn(inGraph)             scalar of modulus 2nd eigenvalue (Shader metric for flow in weighted networks)
#         NodeMetricsFn(inGraph,fylname)     n rows -by- 17 column matrix for global NodeStats, also output as a .csv
#            *** NodeMetricsFn() has code for Eqns 1 to 4 in the AnimBehav paper
#         EgoNetStatsFn(inGraph)             1) n-elelment list of Egonet metrics; 2) FIE 
#         CommStatsFn                        1) List of Community metrics; 2) FIC    
#         GraphStatsFn(inGraph)              19-column list of graph-level weight metrics (part of "GraphStats" dataframe)
#                                               Number of rows is inGraph, RandGraph + Number of Communities
#         FxNStatsFn                         3 variables for GraphStats FIN, FEN, FCN (may add FIx where x=E,C,N)
#         RescaleEdgeWidthsFn(inGraph,a,b)   n-element list of line-widths for plotting (scales fr width a to max width b)     
#         BridgeColorFn                      Function to color inter-community edges orange (vs. black) 
#         NegBinomFitFn(inGraph, infilename) Function to fit a negative binomial distribution to the edge weights
#                                               and output various related plots
#   created 17-Sep-2014 by D.B. McDonald
#
#       last modified 23-Aug-2017 by D.B. McDonald dbmcd@uwyo.edu 
#             -- recalibrated " RescaleEdgeWidthsFn(inGraph,a,b)"      
##############################################################################   

  # Creates graph with weight randomly distributed over the original edges
RandWtGraphFn = function(inGraph)  {
  randWtGrph <- inGraph
  SumWt <- sum(E(randWtGrph)$weight)
  EdgeSet <- rep(1,gsize(randWtGrph))
  EdgeSelector <- c(1:gsize(randWtGrph))
  RandWts <- SumWt - gsize(randWtGrph)
  for (i in 1:RandWts) {
    AddPlace <- sample(EdgeSelector,1)
    EdgeSet[AddPlace] <- EdgeSet[AddPlace]+1 }
  E(randWtGrph)$weight <- EdgeSet
  RandWtGraphOutput <- list(Graph=randWtGrph)
    ## This is the output of the function. I.e., a graph differing from input in having 
     # weights randomly redistributed over the same edge-structure
  return(RandWtGraphOutput) 
}   # End function RandWtGraphFn

  ##  Adds fast-greedy community membership (and color) to the vertex attributes of the input graph 
   ## Revised 9-Sep-16 to update the iGraph function (was  edge_betweenness_community()   )
CommntyDtctFstGrdyFn = function(inGraph){
  NodCnt <- gorder(inGraph)
         ### TRIAL of old edge-between approach    
         ## fastgreedycomm <- cluster_edge_betweenness(inGraph, weights = E(inGraph)$weight, directed = F,
         ## edge.betweenness = TRUE, merges = TRUE, bridges = TRUE, modularity = TRUE, membership = TRUE)
  
    ## iGraph's fast-greedy community detection algorithm finds communities
    ## Notice that we account for weights, so membership depends on weight distribution
    ## Adds V(g)$FstGrdyCommnty & $Color to inGraph for Cols 9, 10 of NodeStats output
   fastgreedycomm <- cluster_fast_greedy(inGraph, merges = TRUE, modularity = TRUE,membership = TRUE, weights = E(inGraph)$weight)
    ## Modularity Q-score of Clauset, Newman, Moore. Finding community structure in very large networks, Phys Rev E 2004
  Qscore <- modularity(fastgreedycomm, weights = E(inGraph)$weight)  ## Row 1, Col 24 of "InputFilename"GraphStats.csv output
  V(inGraph)$FstGrdyCommnty <- fastgreedycomm$membership  ## Col 9 of NodeStats 
    ## find the size of each of the communities 
  FGCommSizeTable <- table(V(inGraph)$FstGrdyCommnty)  
  NumComms <- length(FGCommSizeTable)
    ## CHANGE THIS IF DIFFERENT THRESHOLD DESIRED
  NmBigCom <- length(which(FGCommSizeTable>2))
  BigCommNums <- unlist(dimnames(FGCommSizeTable[which(FGCommSizeTable>3)]))
    ##  sep="" avoids spaces in row names
  BigCommNames <- rep(NA,NmBigCom)
  for (i in 1:NmBigCom) {BigCommNames[i] <- paste("Comm",i,sep="")}
    ## Color vector as V(g)$ attribute for plotting by community membership
  ColorFullSet <- c("yellow","cyan", "lightblue","red", "orange", "purple", "green", "blue", "tan3", "black")
  UsedColors <-c(ColorFullSet[1:NmBigCom])
  VertColSet <- rep(NA,NodCnt)
  
  for (i in 1:NodCnt) { 
    Placer <- which(BigCommNums==V(inGraph)$FstGrdyCommnty[i])
      ## If not in set BigCommNums, returns named integer(0) which has length 0
    ifelse(length(Placer)==0,VertColSet[i] <-"gray", VertColSet[i] <-UsedColors[Placer])
    }
    ## Add the color vector to V(g)$
  V(inGraph)$Color <- VertColSet        ## Col 10 of NodeStats
  FGCommOutput <- list(Graph = inGraph,NameList=BigCommNames, Qscr = Qscore)
  return(FGCommOutput)
}    # End function CommntyDtctFstGrdyFn

  ## Function adds $WClsns, Col 14 of NodeStats output    # inGraph <- MainGraph 
OpsahlWtMetricsFn = function(inGraph) {
  NoadList <- unlist(V(inGraph)$name)   ## Names/IDs (chr/string variable) of the nodes
  NodeOrderList <- which(NoadList %in% V(inGraph)$name==T)  ## str(NodeOrderList)
  WClsnsVector <- rep(NA, length(NoadList))
  EjCnt <- gsize(inGraph)
  EjList <- as_edgelist(inGraph)     ## head(EjList)
     ## Use some routines from Tore Opsahl's tnet R package
  EdgeListW <- cbind(rep(0,EjCnt),rep(0,EjCnt),rep(0,EjCnt))
  for (i in 1:EjCnt){
    EdgeListW[i,1] <- which(NoadList==EjList[i,1])
    EdgeListW[i,2] <- which(NoadList==EjList[i,2])
    EdgeListW[i,3] <- E(inGraph)$weight[i]
  }    ## head(EdgeListW)
    ## Input for undirected has to specify both edge directions
  EdgeList2W <- suppressWarnings(symmetrise_w(EdgeListW))   # str(EdgeList2W) head(EdgeList2W) 
  ClosenessWeights <- suppressWarnings(closeness_w(EdgeList2W))  ## Using Oppsahl's tnet package function "closeness_w()"
  Switch <- identical(as.numeric(ClosenessWeights[,1]),as.numeric(NodeOrderList))
    ## Add the closeness metric of Oppsahl to inGraph          ## Col 13 of NodeStats
  ifelse(Switch, WClsnsVector  <- ClosenessWeights[,2], WClsnsVector  <- NA)
  OpsOutput <- list(WtClsVector = WClsnsVector)
  return(OpsOutput)
}  # End function "OpsahlWtMetricsFn"
  
  # Function to assess Kemeny's constant
    ## calculated by Eqn 4 of Levene & Loizou. 2012. Kemeny's constant and the random surfer. Am.Math.Monthly 109: 741.
    ## K(A,lambda)= sum{i=2 to n} (1/(lambda - lambda-i))
KemenyFn = function(inGraph) {
    ## Note the coercion to class "matrix"
  AdjMat <- as.matrix(get.adjacency(inGraph,type="both",attr="weight"))
  EigVals <- eigen(AdjMat,only.values=T) ## The ordered set of eigenvalues
  OtherEigs <- EigVals$values[-1]
  Lambda <- EigVals$values[1]
  Kmny <- sum(1/(Lambda-OtherEigs)) ## This is the output  <- a scalar
  KmnyOutput <- list(KmnyConst=Kmny)
  return(KmnyOutput)
 }  ## End "KemenyFn" function

  # Function to assess (via Bryan Shader) the eigenvalue modulus metric 
ModEignVal2Fn = function(inGraph) {
  NodCwnt <- gorder(inGraph)
  AdjMat <- as.matrix(get.adjacency(inGraph,type="both",attr="weight"))
  AdjMatNorm <- AdjMat
  AdjMatRowSum <- rowSums(AdjMat)
    ## Now row-normalize 
  for (i in 1:NodCwnt) {AdjMatNorm[i,] <- AdjMatNorm[i,]/AdjMatRowSum[i]
  }
  AdjMatNormEigenVals <- eigen(AdjMatNorm,only.values=T) ## The ordered set of eigenvalues
  ModSet <- Mod(AdjMatNormEigenVals$values)
  ModEigValOut <- list(ModEigVal = ModSet[2])
  return(ModEigValOut)
}  ## End ModEignVal2Fn

## Function NodeMetricsFn     # inGraph <- MainGraph   i <- 1   fylname <- InputcsvFileName
  # Contains the code for Eqns 1 to 4 in the Anim Behav paper
  ## Calculates $NodeID $Degree $EffDgr $DERatio $TotWt $EffWt $ObsEjWtDiv $ExpEjWtDiv  $Conc    
      # Cols    1       2       3       4        5      6      7           8            9
    ##          $FGCmnty $Color  $NrstNbrWtDgr  $WtClstr  [$WClsns]  $Btwns  $EignC  $StClsns
      # Cols    10       11      12             13        14         15      16      17  
   ## and outputs the NodeStats dataframe to a .csv file
NodeMetricsFn = function(inGraph, fylname) {
  ndeCnt         <- gorder(inGraph)
  nodeStats <- as.data.frame(matrix(ncol= 17, nrow= ndeCnt))
    ## Give the NodeStats dataframe row and column names 
    # EffWt is Edge-weight sum/Effective degree, Concentration is (ObsEWD-ExpEWD)/ExpEWD; 
  colnames(nodeStats) <- cbind("NodeID","Degree", "EffDgr", "DEratio","TotWt","EffWt","ObsEdgWtDiv","ExpEdgWtDiv", "Conc", "FGCmnty","Color", "NrstNbrWtDgr","WtClstr","WClsns", "Btwns","EignC","StClsns")
    ## Initialize cells at 0
    nodeStats[]=0
    nodeStats$Degree           <- degree(inGraph)            ## col 2 of nodeStats, used below 
  ndeList         <- unlist(V(inGraph)$name)
  ejList          <- as_edgelist(inGraph)   ## This needs to be local to accommodate Community graphs
  ObsEdgWtDivSet  <- rep(NA,ndeCnt)
  ExpEdgWtDivSet  <- rep(NA,ndeCnt)
  EffDgrSet       <- rep(NA,ndeCnt)
  TotWtSet        <- rep(NA,ndeCnt)
  ConcSet         <- rep(NA,ndeCnt)  
  for (i in 1:ndeCnt) {
     ## EdgeSet is list of edges in which Node i is involved (to and from; by presence in 1st or 2nd column)
    EdgeSet     <- which(ejList[,1] == ndeList[i] | ejList[,2] == ndeList[i])
    WtSet       <- E(inGraph)$weight[EdgeSet] 
    TotWtSet[i] <- unname(strength(inGraph, vids=V(inGraph)[i], loops=F))  ## strength() in iGraph sums the weights of a node's edges
      # Basis for calculating Eqns 1 to 3 (Obs and Exp edge-wt diversity and Effective degree)
    WtPrpn      <- WtSet/sum(WtSet)
      # This is the core code for calculating Eqn 1, Observed edge-weight diversity
    ObsEdgWtDivSet[i] <- 1-sum(WtPrpn^2)   ## Just the edges in which node is involved
      # Core code for Eqn 2, Expected edge-weight diversity
      # Exp edg-wt diversity = (k-1)/k, where k is the degree.
    ExpEdgWtDivSet[i] <- (nodeStats$Degree[i]-1)/(nodeStats$Degree[i])
      # Core code for Eqn 3, Effective degree
    EffDgrSet[i]      <- 1/sum(WtPrpn^2)
      # Core code for Eqn 4, Concentration
      ### To avoid NaN results, the code handles case of ExpEdgWtDiv = 0, to make Conc = 1 (all weight concentrated on the lone edge) 
    ifelse(ExpEdgWtDivSet[i]==0,ConcSet[i] <- 1, ConcSet[i] <- (ExpEdgWtDivSet[i]-ObsEdgWtDivSet[i])/ExpEdgWtDivSet[i])
    ifelse(0.0000001>ConcSet[i]&&ConcSet[i]>-0.0000001,ConcSet[i] <- 0,ConcSet[i] <- ConcSet[i])
      
    } ## End for (i in 1:ndeCnt)
    
  nodeStats$NodeID           <- V(inGraph)$name                        ## col 1  fr main script (upon creation of graph)
  nodeStats$Degree           <- degree(inGraph)                        ## col 2     
  nodeStats$EffDgr           <- EffDgrSet                              ## col 3  
  nodeStats$DEratio          <- nodeStats$Degree/EffDgrSet             ## col 4  
  nodeStats$TotWt            <- strength(inGraph,mode="total")         ## col 5  
  nodeStats$EffWt            <- TotWtSet/EffDgrSet                     ## col 6  
  nodeStats$ObsEdgWtDiv      <- ObsEdgWtDivSet                         ## col 7  
  nodeStats$ExpEdgWtDiv      <- ExpEdgWtDivSet                         ## col 8
  nodeStats$Conc             <- ConcSet                                ## col 9
  nodeStats$FGCmnty          <- V(inGraph)$FstGrdyCommnty              ## col 10  fr function CommntyDtctFstGrdyFn
  nodeStats$Color            <- V(inGraph)$Color                       ## col 11  fr function CommntyDtctFstGrdyFn
    NrstNbrWtDgr                   <- knn(inGraph)                        ## Nearest neighbor weighted degree of Barrat et al. 2005 PNAS
  nodeStats$NrstNbrWtDgr     <- NrstNbrWtDgr$knn                       ## col 12  
      ## iGraph calls cluster cffcnt "transitivity" fr Barrat et al. 2005 PNAS
  nodeStats$WtClstr          <- transitivity(inGraph,type="barrat")    ## col 13
     ## OpsahlWtMetricsFn(inGraph) adds Cols 14 to nodeStats$WtCloseness
   OpsOutput <- OpsahlWtMetricsFn(inGraph)  
  nodeStats$WClsns           <- OpsOutput$WtClsVector                  ## col 14  fr function OpsahlWtMetricsFn
  nodeStats$Btwns            <- betweenness(inGraph, directed = F)     ## col 15
  nodeStats$EignC            <- unname(eigen_centrality(inGraph)[[1]]) ## col 16 (note that eigen_centrality() function ACCOUNTS for weights)
  nodeStats$StClsns          <- closeness(inGraph)                     ## col 17
  
    ##### Output nodeStats  
  filename  <- paste("./Outputs/NodeStats", fylname, sep="")
  filename <- paste(filename,".csv",sep="")
  write.csv(nodeStats, filename, row.names=FALSE)
  
  NdMtrOutput <- list(NodeStatsFrame = nodeStats)    ## nodeStats get sent to FxNStatsFn as an argument
  return(NdMtrOutput)
}  # End function "NodeMetricsFn"

  ## Function to calculate Egonet-level stats   ## inGraph <- MainGraph   infileName <- InputcsvFileName
EgoNetStatsFn = function(inGraph, infileName) {
  NodCnt <- gorder(inGraph)
    ## Create a dataframe to store and output Egonet stats
  EgostatsDF <- as.data.frame(matrix(ncol= 13, nrow= NodCnt))
    rownames(EgostatsDF) <- unlist(V(inGraph)$name)
    colnames(EgostatsDF) <- c("EgoNdeCnt", "EdgCnt","EffDgr","Dnsty","GlobClust","TotWt","WtMean","WtVar","WtSkew","ObsEdgWtDiv", "ExpEdgWtDiv", "EgoConc", "ComBrdg")
    EgostatsDF[] <- NA   ## head(EgostatsDF)
    ## For each Individual in each Egonet, calculate OED and total weight of its edges
    ## NOTE that some individuals will have "pruned" edges that extend outside the Egonet
      ## we will sum OEDIXBit and divide by OEDWtSumBit to get mean EgoOED
  CommBridge <- rep(0, NodCnt)
  OEDIXBit <- rep(NA,NodCnt)  ## OED piece 1 for FIE -- OED(I)*Wt(I)
  OEDWtSumBit <- rep(NA,NodCnt) ## OED piece 2 for FIE -- Total wt of egonet
     
  EgoGraphSet <- make_ego_graph(inGraph, 1, mode = "all")  # str(EgoGraphSet)   i <- 1
     ## EgoGraphSet creates NdeCnt number of Egonet graphs (NdeCnt is # of vertices in MainGraph)
      ## Note that some vertices appear in many different Egonets and that Alters can have pruned edges
  
  for (i in 1:NodCnt) {
               ## Call the ith Egonet via EgoGraphSet[[i]]   i <- 1 
    EgoNet <- EgoGraphSet[[i]]  ## str(EgoNet)  class(EgoGraphSet[[i]])
        ## Assign numeric "names" to the nodes in the Egonet (by order in MainGraph, thus by input adj mat order)
      NodList <- unlist(V(inGraph)$name)     ## Names/IDs (chr/string variable) of the nodes
      NdNameAsNumEgo   <- which(NodList %in% V(EgoNet)$name==T)  ## str(NdNameAsNumEgo)
      EjListEgo        <- as_edgelist(EgoNet)                     ## head(EjListEgo)
      EgoNdCnt         <- gorder(EgoNet)
      EjCntEgo         <- gsize(EgoNet)                           ## this is k (degree)
      EgoTotWt         <- sum(E(EgoNet)$weight)
      EgoWtPrptns      <- E(EgoNet)$weight/EgoTotWt
      ObsEdgWtDiv           <- 1-sum(EgoWtPrptns^2)  
         ## Expected collapses to (k-1)/k where k is edge count (degree)
      ExpEdgWtDiv <- (EjCntEgo -1 )/EjCntEgo
     
       ## Initialize OEDIXBit[i] and OEDWtSumBit[i]
    OEDIXBit[i] <- 0; OEDWtSumBit[i] <- 0 
    EdgeListEgoW     <- cbind(rep(0,EjCntEgo),rep(0,EjCntEgo),rep(0,EjCntEgo))
    for (j in 1:EjCntEgo){
      EdgeListEgoW[j,1] <- which(NodList==EjListEgo[j,1])
      EdgeListEgoW[j,2] <- which(NodList==EjListEgo[j,2])
      EdgeListEgoW[j,3] <- E(EgoNet)$weight[j]
    }           ## head(EdgeListW)     
       ## Calculate OED and WtSum for each node in EgoNet and add to OEDIXBit et al.      ## k <- 1
    for (k in 1:EgoNdCnt) {
        ## Find the edges in which kth node of Egonet is involved
      AlterEdgeSet     <- which(EdgeListEgoW[,1]==NdNameAsNumEgo[k] | EdgeListEgoW[,2]==NdNameAsNumEgo[k])
      AlterWtSet       <- EdgeListEgoW[AlterEdgeSet,3]
      AlterWtPrptn     <- AlterWtSet/sum(AlterWtSet) ## Proportion of node's total 1k wt in kth node's edges
      AlterOED         <- 1-sum(AlterWtPrptn^2)
      OEDIXBit[i]      <- OEDIXBit[i] + AlterOED*sum(AlterWtSet)
      OEDWtSumBit[i]   <- OEDWtSumBit[i] + sum(AlterWtSet)
         ## Flags nodes that have community-bridging edges with a 1 (0 otherwise)
      if( V(inGraph)$FstGrdyCommnty[i]  != V(inGraph)$FstGrdyCommnty[NdNameAsNumEgo[k]] ) CommBridge[i] <- 1
      }    ## End for k in 1:EgoNdCnt
    
    EgostatsDF[i,"EgoNdeCnt"]     <- EgoNdCnt
    EgostatsDF[i,"EdgCnt"]        <- EjCntEgo
    EgostatsDF[i,"EffDgr"]        <- 1/sum(EgoWtPrptns^2)
    EgoMaxEdges                      <- ((EgoNdCnt-1)*EgoNdCnt)/2   ## = half of (NdCnt-1)*NdCnt 
    EgostatsDF[i,"Dnsty"]         <- EjCntEgo/EgoMaxEdges
    EgostatsDF[i,"GlobClust"]     <- transitivity(EgoNet)
    EgostatsDF[i,"TotWt"]         <- EgoTotWt
    EgostatsDF[i,"WtMean"]        <- EgoTotWt/EjCntEgo     ## EjCntEgo is k (degree)
    EgostatsDF[i,"WtVar"]         <- var(EdgeListEgoW[,3])
    EgostatsDF[i,"WtSkew"]        <- skewness(EdgeListEgoW[,3])
    EgostatsDF[i,"ObsEdgWtDiv"]   <- ObsEdgWtDiv
    EgostatsDF[i,"ExpEdgWtDiv"]   <- ExpEdgWtDiv
         ### To avoid NaN results, the code handles case of ExpEdgWtDiv = 0, to make Conc = 1 (wt concentrated on 1 edge)
       ## also partly prevented by using the mean degree as the basis for expectation (for egonet and higher)?
    ifelse(ExpEdgWtDiv==0,Conc <- 1, Conc <- (ExpEdgWtDiv-ObsEdgWtDiv)/ExpEdgWtDiv)
    ifelse(0.0000001>Conc&&Conc>-0.0000001,Conc <- 0,Conc <- Conc)
    EgostatsDF[i,"EgoConc"]       <- Conc
    EgostatsDF[i,"ComBrdg"]       <- CommBridge[i]
  } ## End for (i in 1:NdeCnt)
    ##### Output EgostatsDF
  filename  <- paste("./Outputs/EgoStats",infileName,sep="")
  filename <- paste(filename,".csv",sep="")
  write.csv(EgostatsDF, filename, row.names=T)
  
  MeanOEDIx     <- sum(OEDIXBit)/sum(OEDWtSumBit)
  MeanEgoEED    <- sum(EgostatsDF$ExpEdgWtDiv*EgostatsDF$TotWt)/sum(EgostatsDF$TotWt) 
  Fie           <- (MeanEgoEED - MeanOEDIx)/MeanEgoEED    
  
  EgStsTable <- list(EgoFrame = EgostatsDF, FIE = Fie)
  return(EgStsTable)
}    ### END function EgoNetStatsFn


  ## Calculates Community-level metrics and FIC     inGraph <- MainGraph GrphStats <- GraphStats
CommStatsFn = function(inGraph, GrphStats, NumBigComs){      
  EdgCntByComm     <- rep(NA, NumBigComs)  ## Used in calc of FIC, near end of function
  TotWtbyComm      <- rep(NA, NumBigComs)  ## Used in calc of FIC, near end of function
  OEDIXBitComm     <- rep(NA, NumBigComs)  ## Used in calc of FIC, near end of function
  OEDWtSumBitComm  <- rep(NA, NumBigComs)  ## Used in calc of FIC, near end of function

  for (i in 1:NumBigComs){    
      ## Next 2 lines create the ith Comm graph    i <- 1
    CommList    <- which(V(inGraph)$FstGrdyCommnty==i)  # length(CommList) Numeric names of nodes in Comm i
    CommGraph   <- induced.subgraph(inGraph, CommList, impl="auto")
    AdjMatCheck <- as.matrix(get.adjacency(CommGraph,type="both",attr="weight"))
    CommGraph   <- graph.adjacency(AdjMatCheck,mode="undirected",weighted=T,diag=FALSE)
    
       ## Assign numeric "names" to the nodes in CommGraph (by order in MainGraph, thus by input adj mat order)
    EjListComm        <- as_edgelist(CommGraph)    ## head(EjListComm)
    NoadLst          <- unlist(V(CommGraph)$name)   ## Names/IDs (chr/string variable) of the nodes
    CommNdCnt         <- gorder(CommGraph)
    EjCntCom          <- gsize(CommGraph)
    EdgCntByComm[i]   <- EjCntCom                  ## this is the number of edges, k (degree)
    CommTotWt         <- sum(E(CommGraph)$weight)
    TotWtbyComm[i]    <- CommTotWt                 ## sum is over k weights
    CommWtPrptns      <- E(CommGraph)$weight/CommTotWt
    ObsEdgWtDivComm   <- 1-sum(CommWtPrptns^2)  
    ExpEdgWtDivComm   <- (EjCntCom -1 )/EjCntCom
    
      ## This outputs CommNodeStats as a.csv with CommGraph as input
          ## N.B. Degree differs from overall, because some nodes lose  "external" edges 
    infileName <- paste(InputcsvFileName,"Comm",i,sep="")
    NdMtrcsCommOuput <- NodeMetricsFn(CommGraph, infileName)
         ## could return NodeStatsFrame, but we don't need it here
      ## CommNodeStats <- NdMtrcsCommOuput$NodeStatsFrame      
    
       ## Initialize OEDIXBit[i] and OEDWtSumBit[i]
    OEDIXBitComm[i] <- 0; OEDWtSumBitComm[i] <- 0 
    EdgeListCommW          <- cbind(rep(0,EjCntCom),rep(0,EjCntCom),rep(0,EjCntCom))
    for (j in 1:EjCntCom){
      EdgeListCommW[j,1] <- which(NoadLst==EjListComm[j,1])
      EdgeListCommW[j,2] <- which(NoadLst==EjListComm[j,2])
      EdgeListCommW[j,3] <- E(CommGraph)$weight[j]
    }    
      ## Calculate OED and WtSum for each node in CommGraph and add to OEDIXBit et al.      ## k <- 1
      ## Pieces for calculating FIC
    for (k in 1:CommNdCnt) {
        ## Find the edges in which kth node of CommGraph is involved
      CommIEdgeSet     <- which(EdgeListCommW[,1] == CommList[k] | EdgeListCommW[,2]==CommList[k])
      CommIWtSet       <- EdgeListCommW[CommIEdgeSet,3]
      CommIWtPrptn     <- CommIWtSet/sum(CommIWtSet)
      CommIOED         <- 1-sum(CommIWtPrptn^2)
      OEDIXBitComm[i]  <- OEDIXBitComm[i] +  CommIOED*sum(CommIWtSet)
      OEDWtSumBitComm[i] <- OEDWtSumBitComm[i] + sum(CommIWtSet) 
      }     ## End for k in 1:CommNdCnt
    
    GrStsOutput <- GraphStatsFn(CommGraph)   ## Outputs vector $GrphStVctr
    CommVector <- GrStsOutput$GrphStVctr
    RowPlacer <- i+2 ## 1st Comm in row 3 etc.
    GrphStats[RowPlacer, 1:17] <- CommVector  ## 1:17 because GraphStatsFn calculates 1st 17 cols
    
      ## Save the graph plot to the output file
    a <- 2; b <- 20
    RescaleOutput <- RescaleEdgeWidthsFn(CommGraph, a, b)   ##  RescaleOutput <- list(Graph = inGraph)
    CommGraph <- RescaleOutput$Graph                        ## head(E(CommGraph)$EdgPlotWdth)
    
    CommPlotFilename <- paste(getwd(), "/Outputs/",BigCommNames[i],".pdf", sep="")
    pdf(CommPlotFilename)
    par(bg="gray63")
    plot.igraph(CommGraph,vertex.size=V(CommGraph)$VrtSz, vertex.color=V(CommGraph)$Color, vertex.label.cex=0.5, vertex.label.dist=0.25, vertex.label.degree=-pi/2, edge.color="black",edge.width=E(CommGraph)$EdgPlotWdth)
    dev.off()  ## turn off pdf device
  } ## End for(i in 1:NumBigComs)
    
    ## Finish FIC calculations
  MeanOEDComm <- sum(OEDIXBitComm)/sum(OEDWtSumBitComm)
    End <- NumBigComs + 2
    CommRows <- c(3:End)    ## Index for rows to pull from GraphStats
      ## Now take (EEDcomm vector * CommTotWt vector)/(sum(CommTotWt vector)), the weighted mean EED (weighted by TotWt of Comms)
  MeanCommEED    <- sum(GrphStats[CommRows,"ExpEdgWtDiv"]*GrphStats[CommRows,"TotWt"])/sum(GrphStats[CommRows,"TotWt"]) 
  Fic <- (MeanCommEED - MeanOEDComm)/MeanCommEED
  
  CommStsTable <- list(FIC = Fic, GrStats = GrphStats)
  return(CommStsTable)
} ## End function CommStatsFn

  ## Function to calculate 1st 17 Network-level metrics to the dataframe "GraphStats"
    # for each community and FullOriginal and Rand graphs
    # "NdeCnt","EdgCnt","EffDgr","DEratio","Dnsty","GlobClust","Diam","TotWt","WtVar","WtSD","WtCoV","WtSkew","Kmny","ModEign", ObsEdgWtDiv",
        #   "ExpEdgDiv","Conc"
GraphStatsFn=function(inGraph) {
  ndeCnt            <- gorder(inGraph)                   ## Used later
  edgCnt            <- gsize(inGraph)                    ## Used later
  MaxEdges          <- (ndeCnt*(ndeCnt-1))/2              ## Half of (NdCnt-1)*NdCnt; i.e., number of cells in upper diag of AdjMat
  TotWt             <- sum(E(inGraph)$weight)            ## Used later
  EjWtPrptns        <- E(inGraph)$weight/TotWt           ## Used later
  
  graphStatVector <- rep(NA,17)
  graphStatVector[1]  <- gorder(inGraph)                          ## 1  NdeCnt
  graphStatVector[2]  <- edgCnt                                   ## 2  EdgCnt
  graphStatVector[3]  <- 1/sum(EjWtPrptns^2)                      ## 3  EffDgr
  graphStatVector[4]  <- edgCnt/(1/sum(EjWtPrptns^2))             ## 4  DERatio
  graphStatVector[5]  <- edgCnt/MaxEdges                          ## 5  Dnsty
  graphStatVector[6]  <- transitivity(inGraph)                    ## 6  GlobClust; called "transitivity" in iGraph
  graphStatVector[7]  <- diameter(inGraph)                        ## 7  Diam     ; might use to correct Kemeny and ModEigenval metrics                               
  graphStatVector[8]  <- TotWt                                    ## 8  TotWt
  graphStatVector[9]  <- var(E(inGraph)$weight)                   ## 9  WtVar
  graphStatVector[10] <- sd(E(inGraph)$weight)                    ## 10 WtSD
    NtwrkWtMean         <- sum(E(inGraph)$weight)/edgCnt                           ## intermediate calculation 
  graphStatVector[11] <- sd(E(inGraph)$weight)/NtwrkWtMean        ## 11 WtCoV
  graphStatVector[12] <- skewness(E(inGraph)$weight)              ## 12 WtSkew   ; requires package "moments"
    KmnyOutput          <- KemenyFn(inGraph)                                       ## calls function "KemenyFn"
  graphStatVector[13] <- KmnyOutput$KmnyConst                     ## 13 Kmny
    ModEigValOutput     <- ModEignVal2Fn(inGraph)                                  ## calls function "ModEignVal2Fn" 
  graphStatVector[14] <- ModEigValOutput$ModEigVal                ## 14 ModEign
    NtwrkWtObsDiv     <- 1 - sum(EjWtPrptns^2)
  graphStatVector[15] <- NtwrkWtObsDiv                            ## 15 ObsEdgWtDiv
    NtwrkWtExpDiv     <- (edgCnt - 1)/edgCnt
  graphStatVector[16] <- NtwrkWtExpDiv                            ## 16 ExpEdgDiv
               ## Prevent dividing by 0               
    ifelse(NtwrkWtExpDiv==0, NtwrkWtConc<-1, NtwrkWtConc <-(NtwrkWtExpDiv - NtwrkWtObsDiv)/NtwrkWtExpDiv)
  graphStatVector[17] <- NtwrkWtConc                              ## 17 Conc

  GrStsTable <- list(GrphStVctr = graphStatVector)
  return(GrStsTable)
  }   ## End function "GraphStatsFn"

  ## Calculates the F-statistics for inGraph (cols 18-22 of GraphStats output)
  ## NodeStats, EgoStats and GraphStats are global calls
FxNStatsFn = function(inGraph, NumBgComs){
     ### FIN -- edge-wt variance in Individuals (nodes) compared to Network -- analog of FIS
  OEDnodeVector   <- NodeStats$ObsEdgWtDiv                                     ## NodeStats$ObsEdgWtDiv
  NodeTotWtVector <- NodeStats$TotWt                                           ## head(NodeStats$TotWt)
    ## Note that Edge weights are the analogs for Ni size correction (in PopGen F-stats)
  MeanOEDnode     <- sum(OEDnodeVector * NodeTotWtVector)/sum(NodeTotWtVector)  ## Weighted mean, where weighting is by total weight of the edges on a node
  EEDGlobal       <- (gsize(inGraph) - 1)/gsize(inGraph)                        ## gsize(graph) is count of edges ExpGlobal = (TotalEdges -1)/TotalEdges
  FIN             <- (EEDGlobal - MeanOEDnode)/EEDGlobal                        ## FIN = 0.2133 for LTM156 

    ### Calculate FEN -- variance in Egonets compared to Network 
  OEDegoVector    <- EgoStats$ObsEdgWtDiv
  EgoWtVector     <- EgoStats$TotWt
  MeanOEDego      <- sum(OEDegoVector * EgoWtVector)/sum(EgoWtVector)
    ## Have already calculated global expected diversity as EEDGlobal
  FEN           <- (EEDGlobal - MeanOEDego)/EEDGlobal                           ## FEN = 0.1374 for LTM 156
    ## Now do random graph
  RandOEDvctr     <- RandEgoStats$ObsEdgWtDiv
  RandEgoWtVctr   <- RandEgoStats$TotWt
  RandMeanOEDego  <- sum(RandOEDvctr * RandEgoWtVctr)/sum(RandEgoWtVctr)
  FENrand         <- (EEDGlobal - RandMeanOEDego)/EEDGlobal 
  
    ## FCN -- variance in Communities compared to Network -- analog of FST
  CommEnder     <- NumBgComs+2
  ObsDivInComms <- GraphStats[3:CommEnder,"ObsEdgWtDiv"]
  CommWt        <- GraphStats[3:CommEnder,"TotWt"]
  MeanOEDcomm   <- sum(ObsDivInComms * CommWt)/sum(CommWt)
  CommEdges     <- GraphStats[3:CommEnder,"EdgCnt"]
  EEDGlobal     <- (sum(CommEdges)-1)/sum(CommEdges)
  FCN           <- (EEDGlobal - MeanOEDcomm)/EEDGlobal                          ## FCN = 0.0996 for LTM  
 
  FxNOutput     <- list(Fin=FIN, Fen=FEN, Fenrand = FENrand, Fcn=FCN)
  return(FxNOutput)
}   ## End function FxNStatsFn

  ## Function to rescale the edge weights to get line-widths for plotting [not implemented here, but fire away ...]
   ## Rescales from a to b using the formula  X' = a + ( (X - Xmin)*(b-a) )/(Xmax-Xmin)
RescaleEdgeWidthsFn = function(inGraph, a, b) {
  Wts      <- E(inGraph)$weight
  Xmin     <- min(Wts)
  if (Xmin<1) Wts <- Wts*(1/Xmin)  ## Increase wts for low-weight cases
  Xmax     <- max(Wts)
  ifelse (Xmax<10, MaxWdth <- 10, MaxWdth <- (b - a) )    ## MaxWdth is (b-a), but set to 10 for low-weight cases
  EdgWdths <- a +( (Wts-1)*MaxWdth)/(Xmax-1)
    ## N.B. simpler way to set edge or vertex attributes than set.edge.attribute 
  E(inGraph)$EdgPlotWdth <- EdgWdths
  RescaleOutput <- list(Graph = inGraph)
  return(RescaleOutput)
}   # End function RescaleEdgeWidthsFn

  ## Function to color inter-community edges orange (vs. black) [not implemented here, but fire away ...]
BridgeColorFn = function(inGraph) {
  NoadLst <- unlist(V(inGraph)$name)  ## Names/IDs (chr/string variable) of the nodes
  EjCownt <- gsize(inGraph)           ## Number of (binary) edges in undirected graph
  EDjList <- as_edgelist(inGraph)     ## head(EDjList)
  EdgeIndex <- rep(0,EjCownt)  ## Used to color bridges btwn Communities in plots
  for (i in 1:EjCownt){
    Ind1 <- which(NoadLst[] == EDjList[i,1])
    Ind2 <- which(NoadLst[] == EDjList[i,2])
    if(V(inGraph)$FstGrdyCommnty[Ind1] != V(inGraph)$FstGrdyCommnty[Ind2]) 
    {EdgeIndex[i] <- 1}
  }   ##  sum(EdgeIndex) = 232 in LTM156 w/ fast-greedy = number community-bridging edges  
  ColorIndex <- which(EdgeIndex==1)      ## sum(EdgeIndex)
  E(inGraph)$Color <- rep("black",EjCownt)  ## initialize to black
  E(inGraph)$Color[ColorIndex] <- "gray"
  ColorOutput <- list(Graph = inGraph)
  return(ColorOutput)
}   ## End function BridgeColorFn

  ## Function to calculate fit of edge weights to a negative binomial distribution [MAY CRASH FOR SOME NETWORKS!!!!]
NegBinomFitFn = function(inGraph, infileName){  
  library(fitdistrplus)
  library(stats)
  NodCnt <- gorder(inGraph)
  NumEjs <- length(E(MainGraph)$weight)
  EjList <- get.data.frame(inGraph,what="edges")  
    ## Calculate the max possible number of (binary) edges
  MaxWt <- max(EjList[,3])
  WtVals <- E(inGraph)$weight
  MaxPossEdges <- (NodCnt*(NodCnt - 1)/2) ## Max possible number edges is n*(n-1)/2
    ## Create a list of the NON-ZERO edges with weights from 1 to MaxWt
  IntegerTable <- tabulate(EjList[,3],nbins=MaxWt)    ## head(IntegerTable)
    ## Find the number of edges with 0 weight
  NumMissing <- MaxPossEdges - NumEjs   ## Number of edges of weight 0  log(NumMissing)       log(390)
    ## Find the "missing" frequency of zero-weight edges
  ZeroFreq <- NumMissing/MaxPossEdges
    ## Add the count of zeros to the count of edge weights ≥1
  IntegerZeroTable <- c(NumMissing,IntegerTable)   ## head(IntegerZeroTable)   length(IntegerZeroTable)
    ## Now create a full (from 0 to MaxWt) list of weight proportions 
  FullWtPrptnList <- IntegerZeroTable/MaxPossEdges  ## sum(FullWtPrptnList) length(FullWtPrptnList) head(FullWtPrptnList) 
  
  WtValsZero    <- c(rep(0,NumMissing),WtVals) ## length(WtValsZero)    str(WtValsZero)
  MeanWtValZero <- mean(WtValsZero)            ## Mean weight, INCLUDING the zeros
  VarWtValZero  <- var(WtValsZero)             ## sd(WtValsZero)
  SkewWtValZero <- skewness(WtValsZero)        ## Requires {moments} package 
  
  NegBinomFit <- fitdist(WtValsZero, "nbinom") ## str(NegBinomFit) 
  
   ## qqPlotFilename <- paste(getwd(), "/Outputs/","qqPlot.pdf", sep="")
   ## pdf(qqPlotFilename)
   ## XwZero <- 0:MaxWt                          ## Used in the plotting        
   ## k <- unname(NegBinomFit$estimate[1])       ## This is the size parameter, k
   ## MeanFit <- unname(NegBinomFit$estimate[2]) ## This is the mean parameter, mu
   ## Ybnm <-dnbinom(XwZero, mu = MeanFit, size = k)
   ## qqplot(Ybnm, FullWtPrptnList)
   ## abline(0,1)                                ## a 45-degree reference line is plotted
   ## dev.off()  ## turn off pdf device
  
  ## NEEDS WORK  ppcomp() assesses fitdist fit attempts
  ## ppPlotFilename <- paste(getwd(), "/Outputs/","ppPlot.pdf", sep="")
    ## pdf(ppPlotFilename)
    ## par(bg="gray63")
    ## ppcomp(NegBinomFit)     ## class(NegBinomFit$data)
  ## dev.off()  ## turn off pdf device
  
  NegBinomFitPlotFilename <- paste(getwd(), "/Outputs/","NegBinomFitPlot.pdf", sep="")
    pdf(NegBinomFitPlotFilename)
    plot(NegBinomFit, demp=T)
  dev.off()  ## turn off pdf device
  
  cdfPlotFilename <- paste(getwd(), "/Outputs/","cdfNgBnmPlot.pdf", sep="")
    pdf(cdfPlotFilename)
    cdfcomp(NegBinomFit, xlim=c(0,60))
  dev.off()  ## turn off pdf device
  
    ## Text output of summary(NegBinomFit) plus a few metrics fr edge-wt distribution
  filename  <- paste("./Outputs/NegBnmFitSumm",infileName,sep="")
  filename <- paste(filename,".txt",sep="")
  NegBinomSumm <- summary(NegBinomFit)  ## str(NegBinomSumm)
  capture.output(NegBinomSumm, file = filename)
    ## Note that each line of code gets a new line in output (desirable) 
  VarWtVal <- paste("Var of weights (incl. zeros) is ", round(VarWtValZero,3))
  NumZeros <- paste("Number of zero-weight edges is ", NumMissing)
  SkewOut  <- paste("Skew of edge weight distr. is ", SkewWtValZero)
  write(c(VarWtVal,NumZeros, SkewOut), file = filename, append=TRUE)
}   ## END function NegBinomFitFn

############################################ END  