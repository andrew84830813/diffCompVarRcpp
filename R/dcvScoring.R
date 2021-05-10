#' Computes the Differential Compositional Variation (DCV) Scoring ALgoritmn
#'
#'
#' Scores log-ration using the DCV metric (Hinton (2021))
#'
#' @importFrom magrittr %>%
#' @param logRatioMatrix a n-sample by p-logratio matrix
#' @param includeInfoGain should the infromation gain score be added
#' @param nfolds number of parition used to compute and average scores
#' @param numRepeats number of repeats on the nfold paritions
#' @param seed_ random seed control
#' @param rankOrder should the score be rank ordered by score along with computing the running number of distinct features (Adds non-trivial computational time).  
#'
#' @useDynLib diffCompVarRcpp
#'
#' @examples
#' finalDCV()
#'
#' @references
#' Hinton, A.L., Mucha, P.J., (2021). Simultaneous variable selection and group association testing in sparse high dimensional omics data. XXXX.
#'
#' @return A list containing:\tabular{ll}{
#'    \code{lrs} \tab a p-row log-ratio scoring matrix aggregated across nfolds and numRepeats  \cr
#'    \tab \cr
#'    \code{rawDCV} \tab The raw log-ratio scoring matrix with scores from each metric (Hinton(2021)) unaggreagted \cr
#' }
#' use_pipe(export = TRUE)
#'
#' @export
dcvScores <-
  function(logRatioMatrix,includeInfoGain = T,nfolds = 5,numRepeats = 1,seed_ = 08272008,rankOrder = T){
    
    cvDCV = data.frame()
    for(r in 1:numRepeats){
      
      set.seed(r)
      foldData = caret::createFolds(y = logRatioMatrix[,1],k = nfolds)
      
      for(f in 1:nfolds){
        
        rows = foldData[[f]]
        ####################################################
        ##  select fold data for train and test splilt
        trainData = logRatioMatrix[rows,-1]
        ytrain = factor(logRatioMatrix[rows,1])
        classes = as.character(unique(ytrain))
        #####################################################
        ## compute metrics
        overallMedian = column_median(as.matrix(trainData))
        N_p = nrow(trainData) - n_distinct(classes)
        #Group 1 
        g1 = trainData[ytrain==classes[1],]
        g1Medians = column_median(as.matrix(g1))
        g1Means = colMeans(g1)
        g1Var = matrixStats::colVars(as.matrix(g1))
        names(g1Var) = names(g1Means)
        n1 = nrow(g1)
        #Group 2
        g2 = trainData[ytrain==classes[2],]
        n2 = nrow(g2)
        g2Medians = column_median(as.matrix(g2))
        g2Means = colMeans(g2)
        g2Var = matrixStats::colVars(as.matrix(g2))
        names(g2Var) = names(g2Means)
        n2 = nrow(g2)
        
        ## Brown Forsyth
        num = n1*(g1Medians- overallMedian)^2 + n2*(g2Medians- overallMedian)^2
        denom = column_medianVar(as.matrix(g1),g1Medians) + column_medianVar(as.matrix(g2),g2Medians)
        medianF = N_p*(num / denom)
        medianF = scale(as.vector(medianF))
        
        ## Welchs Tstat
        tstat = abs( (g1Means-g2Means) / sqrt( (g1Var/n1) + (g2Var/n2) ) )
        tstat = scale(tstat)
        
        ## F-ratio
        sm = colMeans(trainData)
        expVar =n1*(g1Means - sm)^2 + n2*(g2Means - sm)^2
        unexpVar = g2Var + g1Var
        fRatio = expVar / unexpVar
        fRatio = scale(fRatio)
        
        ## KS Statistic
        y = data.frame(rbind(g1,g2))
        yy = as.matrix(y)
        levs = unique(ytrain)
        mat1 = yy[ytrain==levs[1],]
        mat2 = yy[ytrain==levs[2],]
        KS.df = data.frame(KS = K_S(mt = mat1,mt2 = mat2))
        KS.df = scale(KS.df)
        
        #Combine
        dfc = cbind(tstat,fRatio,KS.df,medianF)
        #dfc = cbind(tstat,fRatio,KS.df)
        
        
        #Information Gain
        if(includeInfoGain){
          ig = FSelectorRcpp::information_gain(x = trainData,y = ytrain,type = "gainratio" ,discIntegers = T)$importance
          ig[is.na(ig)]=0
          ig = scale(ig)
          dfc = cbind(dfc,ig)
        }
        
        #Scale
        dfc[is.na(dfc)]=0
        dfc[is.infinite(dfc)]=0
        rowmean = rowSums(dfc,na.rm = T)
        dfc = cbind(dfc,rowmean)
        dfc = cbind.data.frame(dfc,f)
        dfc = cbind(dfc,colnames(trainData))
        colnames(dfc) = c("welch_T","classical_F","KS","brownForsythe_F","IG","rowmean","fold","Ratio")
        #colnames(dfc) = c("welch_T","classical_F","KS","IG","rowmean","fold","Ratio")
        
        cvDCV = rbind(cvDCV,dfc)
        message("Fold ", f, " of ", nfolds)
      }
      
    }  
    
    
    if(rankOrder){
      ## Aggregate
      dcv = cvDCV %>% 
        dplyr::group_by(Ratio) %>% 
        dplyr::summarise_all(.funs = mean) 
      
      dcv = data.table::setDT(dcv)[order(-rowmean)]
      
      
      #node strength
      el = data.frame(Ratio = dcv$Ratio,Score = dcv$rowmean)
      el = tidyr::separate(data = el,col = 1,into = c("num","denom"),sep = "___",remove = F)
      nDistinct = c()
      for(i in 1:nrow(el)){
        num = unique(el$num[1:i])
        den = unique(el$denom[1:i])
        nDistinct[i] = dplyr::n_distinct(c(num,den))
      }
      dcv$nDistinct = nDistinct
      
      return(list(lrs = dcv[,c("Ratio","rowmean","nDistinct")],rawDCV = dcv))
    }else{
      dcv = cvDCV %>% 
        dplyr::group_by(Ratio) %>% 
        dplyr::summarise_all(.funs = mean)
      return(list(lrs = dcv[,c("Ratio","rowmean")],rawDCV = dcv))
    }
    
    
    
  }