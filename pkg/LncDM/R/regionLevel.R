regionLevel <-
function(data,indexmethod = c("mean", "median"),classes=c("gene","lincRNA","processed_transcript","pseudogene"))
{
  if(length(classes)>1)
    stop("The classes must have only one element")
  bmatrix = data@bmatrix
  annotation = data@annot
  groupinfo = data@groupinfo
  lincRNA =  data@linc
  processed_transcript =  data@processed
  protein_coding =  data@protein
  pseudogene =  data@pseu

###########################get beta matrix of elements###############################
  if(classes=="gene")
    transAnno=protein_coding
  else
    eval(parse(text=paste("transAnno=",classes,sep="")))
  
  transBeta=matrix(NA, nrow = nrow(transAnno), ncol = ncol(bmatrix))
  for(i in 1:nrow(transAnno))
  {
    cg=unlist(strsplit(transAnno[i,4],","))
    temp=bmatrix[cg,]
    if (length(cg) == 1)
    {
      transBeta[i, ] = temp
    }else
    {
      transBeta[i, ] = apply(temp, 2, eval(indexmethod), na.rm = TRUE)
    }
  }
  colnames(transBeta)=colnames(bmatrix)
  RegionName = unique(transAnno[,6])
  cat("There some information about sub_regions of ",classes,":\n")
  for (i in 1:length(RegionName))
  {
    if(sum(transAnno[,6]==RegionName[i])>0)
      cat(RegionName[i]," region contains:",nrow(transAnno[transAnno[,6]==RegionName[i],]),"GENCODE ",classes," transcript \n")
    else
      cat(RegionName[i]," region contains: 0 GENCODE ",classes," transcript \n")
  }  
  
#############################Visualization of data###################################
  NumOfRegion=NA
  NumOfSites=NA
  for (i in 1:length(RegionName))
  {
    temp=transAnno[transAnno[,6]==RegionName[i],]
    if(sum(transAnno[,6]==RegionName[i])>1)
    {
      NumOfRegion=c(NumOfRegion,nrow(temp))
      NumOfSites=c(NumOfSites,length(unique(unlist(strsplit(as.character(temp[,4]), ",")))))
    }
    if(sum(transAnno[,6]==RegionName[i])==1)
    {
      NumOfRegion=c(NumOfRegion,1)
      NumOfSites=c(NumOfSites,length(unique(unlist(strsplit(as.character(temp[4]), ",")))))  
    }    
  }
  NumOfRegion=NumOfRegion[-1]
  NumOfSites=NumOfSites[-1]
  eval(parse(text = paste("pdf(\"",classes,"Region.pdf\",width=12)" , sep = "")))
  par(mar=c(4,4,4,4))
  barplot(NumOfRegion,names.arg=RegionName,ylim=c(0,max(NumOfRegion)*1.1),width=1,space=0,xlim=c(0,length(RegionName)),ylab="number of transcripts",cex.names=0.7)
  par(new=T,mar=c(4,4,4,4),usr=c(0,length(RegionName),0,max(NumOfSites)*1.1))
  plot(x=1:length(RegionName)-0.5,y=NumOfSites,type="o",pch=3,axes=FALSE,xlab="",ylab="",xlim=c(0,length(RegionName)),ylim=c(0,max(NumOfSites)*1.1))
  axis(4)
  mtext("number of CpG sites",4,2)
  legend("topright",pch=3, legend="CpG sites", horiz=T, lty=1,bty="n")
  legend("topleft",legend="transcript",fill="grey",bty="n")

  sampGroups = unique(groupinfo[,2])
  if (length(sampGroups)==2)
  {
    name=c()
    for (i in 1:length(RegionName))
    {
      eset=transBeta[transAnno[,6]==RegionName[i],]
      c1 = groupinfo[, 2] %in% sampGroups[1]
      c2 = groupinfo[, 2] %in% sampGroups[2]
      if (sum(transAnno[,6]==RegionName[i])>1)
      {
        titl=RegionName[i]
        if(titl=="5'UTR")
          titl="5UTR"
        if(titl=="3'UTR")
          titl="3UTR"
        eval(parse(text=paste(sampGroups[1],"_",titl,"= apply(eset[, c1], 1, mean)",sep="")))
        eval(parse(text=paste(sampGroups[2],"_",titl,"= apply(eset[, c2], 1, mean)",sep="")))
        name=c(name,paste(sampGroups[1],"_",titl,sep=""),paste(sampGroups[2],"_",titl,sep=""))
      }
      if (sum(transAnno[,6]==RegionName[i])==1)
      {
        titl=RegionName[i]
        if(titl=="5'UTR")
          titl="5UTR"
        if(titl=="3'UTR")
          titl="3UTR"
        eval(parse(text=paste(sampGroups[1],"_",titl,"= eset[c1]",sep="")))
        eval(parse(text=paste(sampGroups[2],"_",titl,"= eset[c2]",sep="")))
        name=c(name,paste(sampGroups[1],"_",titl,sep=""),paste(sampGroups[2],"_",titl,sep=""))
      } 
    }
    eval(parse(text=paste("boxplot(",paste(name,collapse =","),",col=c(\"chocolate2\",\"cyan3\")",",ylab=\"beta\",xlab=NULL,main=classes)",sep = "")))
    x=c(1:length(name))
    y=par("usr")[3]-0.01
    text(x, y, labels=name, adj=1, srt=45, xpd=TRUE)
    legend("topright", legend = unique(sampGroups), text.col = c("chocolate2","cyan3"))
  }
  dev.off()
################################return result########################################
  setClass("RegionMethy450", representation(groupinfo = "data.frame",annotation = "matrix",transAnno = "matrix",transBeta="matrix"), where = topenv(parent.frame()))
  RegionLevel = new("RegionMethy450", groupinfo = groupinfo,annotation = annotation,transAnno = transAnno,transBeta=transBeta)
  cat("\nA RegionMethy450 class is created and the slotNames are:\n", slotNames(RegionLevel), "\n")
  return(RegionLevel)
}
