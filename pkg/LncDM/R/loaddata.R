loaddata <-
function(fileDir,is_beta=FALSE,beta_method=c("M/(M+U)","M/(M+U+100)"),groupfile,samplefilter = FALSE,contin=c("ON","OFF"),samplefilterperc = 0.75,XYchrom = c(FALSE, "X","Y", c("X", "Y")),sitefilter = FALSE, sitefilterperc = 0.75,filterDecetP=0.05, normalization  = FALSE,transfm = c(FALSE, "arcsinsqr", "logit"),snpfilter=c(FALSE,"within_10","prob_snp"),gcase="case",gcontrol="control",skip=0,imputation=c("mean","min","knn"),knn.k=10)
{  
  densityBeanPlot<-function(dat,sampGroups=NULL,main = NULL,pal = c("chocolate2","cyan3","darkgoldenrod2","green3"), numPositions = 10000)
  {
    #require(beanplot)
    #require(reshape)
    n <- ncol(dat)
    if (is.null(main))
    main <- "Beta"
    if (is.null(sampGroups))
    sampGroups <- rep(1, n)
    sampGroups <- as.factor(sampGroups)
    col <- lapply(sampGroups, function(x) rep(pal[x], 4))
    if (is.null(numPositions))
    idx <- 1:dim(dat)[1]
    else
     idx <- sample(nrow(dat), numPositions)
     x <- melt(dat[idx, ], varnames = c("cpg", "sample"))
     o <- order(unique(colnames(dat)))
     beanplot(value ~ sample, horizontal = TRUE, what = c(0, 1, 1, 0), log = "", las = 1,  xlab = "Beta",main = main, col = col[o], data = x, cex.lab = 0.9,beanlinewd = 1,border = NA)
     abline(h = 1:(n + 1) - 0.5, lty = 3, col = "grey70")
     legend("topright", legend = levels(sampGroups), text.col = pal,text.width=0.05)
  }
  
  densityPlot<-function (dat, sampGroups=NULL, main = "", xlab = "Beta", pal = c("chocolate2","cyan3","darkgoldenrod2","green3"), xlim, ylim, add = TRUE, tag = TRUE, ...)
  {
    d <- lapply(dat,density, na.rm = TRUE)
    if (missing(ylim))
       ylim <- range(sapply(d, function(i) range(i$y)))
    if (missing(xlim))
       xlim <- range(sapply(d, function(i) range(i$x)))
    if (is.null(sampGroups))
    {
      sampGroups <- rep(1, length(d))
    }
    else if (length(sampGroups) == 1)
    {
      sampGroups <- rep(sampGroups, length(d))
    }
    sampGroups <- as.factor(sampGroups)
    if (add)
    {
      plot(0, type = "n", ylim = ylim, xlim = xlim, ylab = "Density", xlab = xlab, main = main, ...)
      abline(h = 0, col = "grey80")
    }
    for (i in 1:length(d))
    {
      lines(d[[i]], col = pal[sampGroups[i]])
    }
    if (tag & length(levels(sampGroups)) > 1)
      legend("topright", legend = levels(sampGroups), text.col = pal,text.width=0.1)
   }


##############################load annotation########################################
  
  annotation_path <- system.file("extdata", "annotation",package="LncDM")
  load(paste(annotation_path, '/annotation.Rdata', sep=''))   #the annotation of CpG sites     
  load(paste(annotation_path, '/CpGannotation.Rdata', sep=''))
  cat("Total CpG sites without any filtering are:", nrow(CpGannotation),"\n")
  #the sites after filter XYchrom and snp
  if (XYchrom[1] != FALSE)
  {
    chr = CpGannotation[, "CHR_38"]
    index = which(chr %in% XYchrom)
    good_chrom = rownames(CpGannotation)[-index]
    cat(length(index), "sites on chr", XYchrom[-1], "are removed\n")
  }else
  {
    good_chrom = rownames(CpGannotation)
  }

  if(snpfilter[1])
  {
    if (snpfilter[2]=="prob_snp")
    {
      snpsites=unique(which(!CpGannotation[,"Probe_SNPs_37"]==""),which(!CpGannotation[,"Probe_SNPs_10_37"]==""))
    }else
    {
      snpsites= which(!CpGannotation[,"Probe_SNPs_10_37"]=="")
    }
    good_sites = rownames(CpGannotation)[-snpsites]
    cat(length(snpsites), "sites contain snps and removed.\n")
  }else
  {
    good_sites=rownames(CpGannotation)
  }
  good_sites=intersect(good_sites,good_chrom)

####################read the pheno and methylation data##############################

  #load the pheno file   
  group = data.frame(read.delim(groupfile, sep = "\t", header = TRUE))
  group[,1]=as.character(group[,1])
  group[,2]=as.character(group[,2])
  rownames(group)=group[,1]
  #filename
  fileName=paste(as.character(group[,1]),".txt",sep="")
  sampleID= as.character(group[,1])

  fileList=list.files(fileDir)  
  temp=fileName[fileName %in% fileList]
  sampleID=sampleID[fileName %in% temp]
  group=group[sampleID,]
  fileName=temp
  eval(parse(text=paste("dir = paste(\"./",fileDir,"/\",fileName,sep=\"\")",sep="")))

  betaMatrix=matrix(,length(good_sites),)
  avgPval=NA;
  if(!is_beta)
  {
    for (i in 1:length(dir))
    {
      temp=read.table(file = dir[i],header=F,sep="\t",skip=skip)
      #the cols are CpG ID,methylation, unmethylation, Pvalue
      temp[,1]=as.character(temp[,1])
      temp[,2]=as.numeric(as.character(temp[,2]))
      temp[,3]=as.numeric(as.character(temp[,3]))
      temp[,4]=as.numeric(as.character(temp[,4]))
      avgPval=c(avgPval,sum(temp[,4]>filterDecetP)/length(temp[,4])*100)
      if (beta_method == "M/(M+U)")
      {
        beta=temp[,2]/(temp[,2]+temp[,3])
      }else
      {
        beta=temp[,2]/(temp[,2]+temp[,3]+100)
      }
      if(samplefilter | sitefilter)
      {
        beta[temp[,4]>filterDecetP]=NA
      }
      names(beta)=temp[,1]
      beta=beta[good_sites]
      betaMatrix=cbind(betaMatrix,beta)
    }
  }else
  {
    for (i in 1:length(dir))
    {
      temp=read.table(file = dir[i],header=F,sep="\t",skip=skip)
      #the cols are CpG ID and beta
      temp[,1]=as.character(temp[,1])
      temp[,2]=as.numeric(as.character(temp[,2]))
      beta=temp[,2] 
      names(beta)=temp[,1]
      beta=beta[good_sites]
      betaMatrix=cbind(betaMatrix,beta)
    }
  }
  betaMatrix=betaMatrix[,-1]
  avgPval=avgPval[-1]
  colnames(betaMatrix)=sampleID
  cat("Total samples are:", ncol(betaMatrix), "\n")

####################filter the uneffective sample and sites##########################

  #filter samples
  if(samplefilter)
  {
    goodsample = colSums(!is.na(betaMatrix)) >= samplefilterperc * nrow(betaMatrix)  
    betaMatrix = betaMatrix[, goodsample]
    cat(length(goodsample) - ncol(betaMatrix), "samples removed with at least", (1-samplefilterperc) * 100, "percentage sites having pvalue greater than", filterDecetP, "\n")
    group = group[goodsample, ]
    sampleID=sampleID[goodsample]
  }
  
  #filter sites
  linc_cgList=unique(unlist(strsplit(as.character(lincRNA[,4]), ",")))
  processed_cgList=unique(unlist(strsplit(as.character(processed_transcript[,4]), ",")))
  protein_cgList=unique(unlist(strsplit(as.character(protein_coding[,4]), ",")))
  pseudogene_cgList=unique(unlist(strsplit(as.character(pseudogene[,4]), ",")))  
  cgList=c(linc_cgList,processed_cgList,protein_cgList,pseudogene_cgList)
  good_sites=intersect(good_sites,cgList)

  if (sitefilter)
  {
    if (contin=="OFF")
    {
      gcaseBeta=betaMatrix[,sampleID[group[sampleID,2]==gcase]]
      gcontrolBeta=betaMatrix[,sampleID[group[sampleID,2]==gcontrol]]
      gcase_loci = rownames(gcaseBeta)[rowSums(!is.na(gcaseBeta)) >= sitefilterperc * ncol(gcaseBeta)]
      gcontrol_loci = rownames(gcontrolBeta)[rowSums(!is.na(gcontrolBeta)) >= sitefilterperc * ncol(gcontrolBeta)]
      good_loci=intersect(gcase_loci,gcontrol_loci)
    }else
    {
      good_loci = rownames(betaMatrix)[rowSums(!is.na(betaMatrix)) >= sitefilterperc * ncol(betaMatrix)]
    }
    cat(nrow(betaMatrix) - length(good_loci), "sites had at least", (1-sitefilterperc) * 100, "% samples with pvalue great than",filterDecetP, "and are removed\n")
  }else
  {
    good_loci = rownames(betaMatrix)
  }
  good_loci=intersect(good_loci,good_sites)
  betaMatrix=betaMatrix[good_loci,]

  CpGannotation=CpGannotation[rownames(betaMatrix),]
  groupinfo=group
##########################solve the beta matrix######################################
  
  #Fill the missing values
  if(imputation=="knn")
  {
    #require(impute)
    imputed=impute.knn(as.matrix(betaMatrix) ,k = knn.k)
    betaMatrix=imputed$data
  }else
  {
    if(imputation=="mean")
    {
      if(contin=="OFF")
      {
        betaMatrix=t(apply(betaMatrix,1,FUN=function(x){x[is.na(x) & group[sampleID,2]==gcase]=mean(x[!is.na(x) & group[sampleID,2]==gcase])
             return(x)
           }))
        betaMatrix=t(apply(betaMatrix,1,FUN=function(x){x[is.na(x) & group[sampleID,2]==gcontrol]=mean(x[!is.na(x) & group[sampleID,2]==gcontrol])
             return(x)
           }))
      }else
      {
        betaMatrix=t(apply(betaMatrix,1,FUN=function(x){x[is.na(x)]=mean(x[!is.na(x)])
             return(x)
           }))
      }  
    }else
    {
      if(contin=="OFF")
      {
        betaMatrix=t(apply(betaMatrix,1,FUN=function(x){x[is.na(x) & group[sampleID,2]==gcase]=min(x[!is.na(x) & group[sampleID,2]==gcase])
             return(x)
           }))
        betaMatrix=t(apply(betaMatrix,1,FUN=function(x){x[is.na(x) & group[sampleID,2]==gcontrol]=min(x[!is.na(x) & group[sampleID,2]==gcontrol])
             return(x)
           }))
      }else
      {
        betaMatrix=t(apply(betaMatrix,1,FUN=function(x){x[is.na(x)]=min(x[!is.na(x)])
             return(x)
           }))
      } 
    }
    
  }
  
  #Standardization
  if (normalization)
  {
    #require(preprocessCore)
    betaMatrix = normalize.quantiles(as.matrix(betaMatrix))
    colnames(betaMatrix) = colnames(betaMatrix)
    rownames(betaMatrix) = rownames(betaMatrix)
    cat("Quantile normalization Performed\n")
  }

  #transform beta
  if (transfm[1])
  {
    if (transfm[2] == "arcsinsqr")
    {
      betaMatrix = asin(sqrt(betaMatrix))
      cat("Transfer beta matrix by the arcsin square root\n")
    }
    if (transfm[2] == "logit")
    {
      betaMatrix[betaMatrix == 0] <- min(betaMatrix[betaMatrix > 0], 0.001)/10
      betaMatrix[betaMatrix == 1] <- max(betaMatrix[betaMatrix < 1], 0.999) + (1 - max(betaMatrix[betaMatrix< 1], 0.999))/100
      betaMatrix = log(betaMatrix/(1 - betaMatrix))
      cat("Transfer beta matrix by the logit transformation \n")
    }
  }
  colnames(betaMatrix)=sampleID
#############################Data visualization######################################
  if(!is_beta){
  #plot barplots of effective sites percent of samples 
  pdf("detect pvalue.pdf",width=24)
  ylab=paste("% of detect pvalue >",filterDecetP,sep="")
  print(avgPval)
  barplot(avgPval, ylab = ylab,ylim=c(0,100))
  abline(h=(1-samplefilterperc)*100)
  dev.off() 
  }  
  
  #plot heatmap of sites and four transcripts
  #require(gplots)
  
  linc_cgList=intersect(good_loci,linc_cgList)  
  processed_cgList=intersect(good_loci,processed_cgList)
  protein_cgList=intersect(good_loci,protein_cgList)
  pseudogene_cgList=intersect(good_loci,pseudogene_cgList)  
  trans_lable=c("linc_cgList","processed_cgList","protein_cgList","pseudogene_cgList")
  trans_name=c("lincRNA","processed_transcript","protein_coding","pseudogene")
    
  pdf("heatmap of CpG sites.pdf")
  
  mads=apply(betaMatrix,1,mad)
  heatmapData=betaMatrix[rev(order(mads))[1:5000],]
  colnames(heatmapData)= group[colnames(heatmapData), 2]
  col=c(rgb(0,191,0,maxColorValue =255),rgb(35,245,26,maxColorValue =255),rgb(145,249,12,maxColorValue =255),rgb(239,239,0,maxColorValue =255),rgb(249,237,167,maxColorValue =255),rgb(177,177,255,maxColorValue =255),rgb(141,141,255,maxColorValue =255),rgb(104,104,255,maxColorValue =255),rgb(66,66,255,maxColorValue =255),rgb(0,0,255,maxColorValue =255))    
  if (contin=="OFF")
  {
    #all sites
    ColSideColors=rep("white",ncol(heatmapData))  
    ColSideColors[colnames(heatmapData)==gcase]= "chocolate2"
    ColSideColors[colnames(heatmapData)==gcontrol]= "cyan3"
    hv <- heatmap.2(heatmapData, dendrogram="none",col = col, trace="none", margins = c(5,5),xlab = "Samples", ylab =  "CpG Sites", main = "heatmap of CpG sites",labRow=FALSE,labCol=FALSE,ColSideColors=ColSideColors,density.info="density",key.title = "color key",key.xlab = "beta")
    legend("topright",legend=c(gcase,gcontrol),fill=c("chocolate2","cyan3"),bty="n")
    #four transcripts
    for (j in 1:length(trans_lable))
    { 
      eval(parse(text=paste("temp=betaMatrix[",trans_lable[j],",]",sep="")))
      mads=apply(temp,1,mad)
      heatmapData=temp[rev(order(mads))[1:min(5000,nrow(temp))],]
      colnames(heatmapData)= group[colnames(heatmapData), 2]
      ColSideColors=rep("white",ncol(heatmapData))  
      ColSideColors[colnames(heatmapData)==gcase]= "chocolate2"
      ColSideColors[colnames(heatmapData)==gcontrol]= "cyan3"
      hv <- heatmap.2(heatmapData, dendrogram="none",col = col, trace="none", margins = c(5,5),xlab = "Samples", ylab =  "CpG Sites", main =trans_name[j] ,labRow=FALSE,labCol=FALSE,ColSideColors=ColSideColors,density.info="density",key.title = "color key",key.xlab = "beta")
      legend("topright",legend=c(gcase,gcontrol),fill=c("chocolate2","cyan3"),bty="n")
    }
  }else
  {
    ColSideColors=rep("white",ncol(heatmapData))
    Col.col=rainbow(length(unique(group[, 2])))
    sampGroups=as.factor(group[, 2])
    for(i in 1:length(ColSideColors))
    {
      ColSideColors[i]=Col.col[sampGroups[i]]
    }
    hv <- heatmap.2(heatmapData, dendrogram="none",col = col, trace="none", margins = c(5,5),xlab = "Samples", ylab =  "CpG Sites", main = "heatmap of CpG sites",labRow=FALSE,labCol=FALSE,ColSideColors=ColSideColors,density.info="density",key.title = "color key",key.xlab = "beta")
    legend("topright",legend=levels(sampGroups),fill=Col.col,bty="n") 
    for (j in 1:length(trans_lable))
    {
      eval(parse(text=paste("temp=betaMatrix[",trans_lable[j],",]",sep="")))
      mads=apply(temp,1,mad)
      heatmapData=temp[rev(order(mads))[1:min(5000,nrow(temp))],]
      colnames(heatmapData)= group[colnames(heatmapData), 2]
      ColSideColors=rep("white",ncol(heatmapData))
      Col.col=rainbow(length(unique(group[, 2])))
      sampGroups=as.factor(group[, 2])
      for(i in 1:length(ColSideColors))
      {
        ColSideColors[i]=Col.col[sampGroups[i]]
      }
      hv <- heatmap.2(heatmapData, dendrogram="none",col = col, trace="none", margins = c(5,5),xlab = "Samples", ylab =  "CpG Sites", main = trans_name[j],labRow=FALSE,labCol=FALSE,ColSideColors=ColSideColors,density.info="density",key.title = "color key",key.xlab = "beta")
      legend("topright",legend=levels(sampGroups),fill=Col.col,bty="n")
    }
  }
  dev.off()
    
  #density plots 
  if (contin=="OFF")
  {
    pdf("density plot.pdf",width=12)
    sampGroups = unique(groupinfo[,2])
    c1 = groupinfo[, 2] %in% sampGroups[1]
    c2 = groupinfo[, 2] %in% sampGroups[2]
    con=betaMatrix[,c1]
    trt=betaMatrix[,c2]
    dim(con)=c(ncol(con)*nrow(con),1)
    dim(trt)=c(ncol(trt)*nrow(trt),1)    
    eval(parse(text=paste("eset=list(",sampGroups[1],"=con,",sampGroups[2],"=trt)",sep="")))     
    pal=c("chocolate2","cyan3","darkgoldenrod2","green3")
    densityPlot(dat=eset, sampGroups=sampGroups, main = "densityPlot of CpG sites", xlab = "Beta", pal = pal, add = TRUE, tag = TRUE)
    boxplot(eset,ylab="beta",main="box plot",boxwex=0.5,col=c("chocolate2","cyan3"))
     
    #eset=betaMatrix
    #colnames(eset) = group[,2]
    #sampGroups=group[,2]
	con=rowMeans(betaMatrix[,c1])
    trt=rowMeans(betaMatrix[,c2])
    eset=cbind(con,trt)
	colnames(eset)=sampGroups
    densityBeanPlot(dat=eset,sampGroups=sampGroups,main = "densityBeanPlot of CpG sites",pal = c("chocolate2","cyan3","darkgoldenrod2","green3"), numPositions = round(nrow(betaMatrix)/(ncol(eset)*10)))
    dev.off()
  }
  
#########################revise the annotation#######################################

  temp=rownames(betaMatrix)
  linc_cgList=unique(unlist(strsplit(as.character(lincRNA[,4]), ",")))
  processed_cgList=unique(unlist(strsplit(as.character(processed_transcript[,4]), ",")))
  protein_cgList=unique(unlist(strsplit(as.character(protein_coding[,4]), ",")))
  pseudogene_cgList=unique(unlist(strsplit(as.character(pseudogene[,4]), ",")))

  linc_cgList=linc_cgList[linc_cgList %in% temp]
  processed_cgList=processed_cgList[processed_cgList %in% temp]
  protein_cgList=protein_cgList[protein_cgList %in% temp]
  pseudogene_cgList=pseudogene_cgList[pseudogene_cgList %in% temp]
  
  processed=strsplit(as.character(processed_transcript[,4]), ",")
  temp=lapply(processed,FUN=function(x){x[x %in% processed_cgList]})
  temp1=temp[temp != "character(0)"]
  processed_transcript=processed_transcript[temp != "character(0)",]
  processed=unlist(lapply(temp1,FUN=function(x){paste(x,collapse=",")}))
  processed_transcript[,4]=processed

  protein=strsplit(as.character(protein_coding[,4]), ",")
  temp=lapply(protein,FUN=function(x){x[x %in% protein_cgList]})
  temp1=temp[temp != "character(0)"]
  protein_coding=protein_coding[temp != "character(0)",]
  protein=unlist(lapply(temp1,FUN=function(x){paste(x,collapse=",")}))
  protein_coding[,4]=protein

  linc=strsplit(as.character(lincRNA[,4]), ",")
  temp=lapply(linc,FUN=function(x){x[x %in% linc_cgList]})
  temp1=temp[temp != "character(0)"]
  lincRNA=lincRNA[temp != "character(0)",]
  linc=unlist(lapply(temp1,FUN=function(x){paste(x,collapse=",")}))
  lincRNA[,4]=linc

  pseu=strsplit(as.character(pseudogene[,4]), ",")
  temp=lapply(pseu,FUN=function(x){x[x %in% pseudogene_cgList]})
  temp1=temp[temp != "character(0)"]
  pseudogene=pseudogene[temp != "character(0)",]
  pseu=unlist(lapply(temp1,FUN=function(x){paste(x,collapse=",")}))
  pseudogene[,4]=pseu

###########################return the result#########################################
   setClass("LincMethy450", representation(bmatrix = "matrix",annot = "matrix", groupinfo = "data.frame", linc="matrix",processed="matrix",protein="matrix",pseu="matrix"),where = topenv(parent.frame()))
   Methy450 = new("LincMethy450", bmatrix = as.matrix(betaMatrix), annot = as.matrix(CpGannotation),linc=as.matrix(lincRNA),processed=as.matrix(processed_transcript),protein=as.matrix(protein_coding),pseu=as.matrix(pseudogene), groupinfo = groupinfo)
   return(Methy450)
}
