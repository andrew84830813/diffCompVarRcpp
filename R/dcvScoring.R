#' Computes the Differential Compositional Variation (DCV) Scoring ALgoritmn
#'
#'
#' Scores log-ration using the DCV metric (Hinton (2021))
#'
#' @importFrom magrittr %>%
#' @importFrom data.table .SD
#' @importFrom foreach %dopar%
#' @param logRatioMatrix a n-sample by p-logratio matrix
#' @param includeInfoGain should the infromation gain score be added
#' @param nfolds number of parition used to compute and average scores
#' @param numRepeats number of repeats on the nfold paritions
#' @param seed_ random seed control
#' @param rankOrder should the score be rank ordered by score along with computing the running number of distinct features (Adds non-trivial computational time).  
#'
#' @useDynLib diffCompVarRcpp
#'
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
  function(logRatioMatrix,includeInfoGain = T,nfolds = 1,numRepeats = 1,seed_ = 08272008,rankOrder = T){
    
    Ratio = NULL
    f = NULL
    label = data.frame(logRatioMatrix)[,1]
    cvDCV = data.frame()
    c = combinat::combn2(as.character(unique(label)))
    for(jz in 1:nrow(c)){
      
      cc = c[jz,]
      ph = logRatioMatrix[label%in%cc,]
      pwComparison = paste0(c[jz,1],"_",c[jz,2])
      
      for(r in 1:numRepeats){
        
        set.seed(r)
        foldData = caret::createFolds(y = ph[,1],k = nfolds)
        message("Compute DCV Scores ",pwComparison,"-(",jz," of ",nrow(c),")")
        cv_DCV = foreach::foreach(f = 1:nfolds,.combine = rbind)%dopar%{
          
          rows = foldData[[f]]
          ####################################################
          ##  select fold data for train and test splilt
          trainData = ph[rows,-1]
          ytrain = factor(ph[rows,1])
          classes = as.character(unique(ytrain))
          #####################################################
          
          #### --------------------*
          ##  compute metrics 
          #### --------------------*
          overallMedian = diffCompVarRcpp::column_median(as.matrix(trainData))
          N_p = nrow(trainData) - dplyr::n_distinct(classes)
          
          ## Group 1 
          g1 = trainData[ytrain==classes[1],]
          g1Medians = diffCompVarRcpp::column_median(as.matrix(g1))
          g1Means = diffCompVarRcpp::column_mean(g1) # g1Means = colMeans(g1)
          g1Var = diffCompVarRcpp::column_medianVar(as.matrix(g1),g1Means)  # g1Var = matrixStats::colVars(as.matrix(g1))
          names(g1Var) = names(g1Means)
          n1 = nrow(g1)
          
          ## Group 2
          g2 = trainData[ytrain==classes[2],]
          n2 = nrow(g2)
          g2Medians = diffCompVarRcpp::column_median(as.matrix(g2))
          g2Means = diffCompVarRcpp::column_mean(g2)  #g2Means = colMeans(g2)
          g2Var = diffCompVarRcpp::column_medianVar(as.matrix(g2),g2Means)  #g2Var = matrixStats::colVars(as.matrix(g2))
          names(g2Var) = names(g2Means)
          n2 = nrow(g2)
          
          ## Brown Forsyth
          num = n1*(g1Medians- overallMedian)^2 + n2*(g2Medians- overallMedian)^2
          denom = diffCompVarRcpp::column_medianVar(as.matrix(g1),g1Medians) + diffCompVarRcpp::column_medianVar(as.matrix(g2),g2Medians)
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
          KS.df = data.frame(KS = diffCompVarRcpp::K_S(mt = mat1,mt2 = mat2))
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
          message("Fold ",f, " of ",nfolds)
          dfc$comp = pwComparison
          dfc
          
        }
        cvDCV = rbind(cvDCV,cv_DCV)
        
      }  
    }
    
    
    
    if(rankOrder){
      ## Aggregate
      dcv = data.table::setDT(cvDCV[,-9])[,by = Ratio,lapply(.SD, mean)]
      dcv = data.table::setDT(dcv)[order(-rowmean)]
      
      message("Compute number of distinct parts")
      #node strength
      el = data.frame(Ratio = dcv$Ratio,Score = dcv$rowmean)
      el = tidyr::separate(data = el,col = 1,into = c("num","denom"),sep = "___",remove = F)
      mxd = dplyr::n_distinct(c(el$num,el$denom))
      nDistinct = rep(mxd,nrow(el))
      for(i in 1:nrow(el)){
        num = unique(el$num[1:i])
        den = unique(el$denom[1:i])
        nDistinct[i] = dplyr::n_distinct(c(num,den))
        if(nDistinct[i]==mxd){
          break
        }
      }
      dcv$nDistinct = nDistinct
      
      return(list(lrs = subset(dcv,select = c("Ratio","rowmean","nDistinct")),rawDCV = dcv))
      
    }else{
      dcv = data.table::setDT(cvDCV[,-9])[,by = Ratio,lapply(.SD, mean)]
      dcv = data.table::setDT(dcv)[order(-rowmean)]
      
      return(list(lrs = dcv[,c("Ratio","rowmean")],rawDCV = dcv,cv_DCV = cvDCV))
    }
    
    
    
  }
