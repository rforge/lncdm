dmr <- function(data,contin=c("ON","OFF"),classes=c("lincRNA","gene","processed_transcript","pseudogene"),testmethod = c("wilcox","limma", "t.test", "satterthwaite"), Padj = c("holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr", "none"), gcase = "case", gcontrol = "control", paired = FALSE,rawpcut = 0.05, adjustpcut = 0.05, betadiffcut = 0.3,num,sole=FALSE)
{ 
  test<-function(eset,testmethod = c("wilcox","limma", "t.test", "satterthwaite"), Padj = c("holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr", "none"), groupinfo, gcase = "case", gcontrol = "control", paired = FALSE)
  {
      grouplev = groupinfo[, 2]
      caseind = which(grouplev %in% gcase)
      controlind = which(grouplev %in% gcontrol)

      temp=apply(eset,1,FUN=function(x){length(table(x))})
      eset=eset[temp>1,]
      
      if(paired == TRUE)
      {
         lev1 = caseind[order(groupinfo[caseind, 2])]
         lev2 = controlind[order(groupinfo[controlind, 2])]
      }
      else
      {
           lev1 = caseind
           lev2 = controlind
      }
      eset = eset[, c(lev1, lev2)]


       if (testmethod == "wilcox")
       {
          cat("Performing Wilcoxon test...\n")
          testout = apply(eset, 1, function(x) {
             wilcox.test(x[1:length(lev1)], x[(length(lev1) + 1):ncol(eset)], paired = paired)$p.value
            })
        }

        if (testmethod == "limma")
        {
            #require(limma)
            cat("Performing limma...\n")
            TS = as.factor(c(rep("T", length(lev1)), rep("C", length(lev2))))
            SS = as.factor(rep(1:length(lev1), 2))

            if (paired == FALSE)
            {
              design = model.matrix(~0 + TS)
              rownames(design) = colnames(eset)
              colnames(design) = c("C", "T")
              fit = lmFit(eset, design)
              cont.matrix = makeContrasts(comp = T - C, levels = design)
              fit2 = contrasts.fit(fit, cont.matrix)
              fit2 = eBayes(fit2)
              result1 = topTable(fit2, coef = 1, adjust.method = Padj, number = nrow(fit2))
            }
            else
            {
              design = model.matrix(~SS + TS)
              fit = lmFit(eset, design)
              fit2 = eBayes(fit)
              result1 = topTable(fit2, coef = "TST", adjust.method = Padj, number = nrow(fit2))
            }
            testout = result1[match(rownames(eset), rownames(result1)), "P.Value"]
        }

        if (testmethod == "t.test")
        {
           cat("Performing t.test...\n")
           testout = apply(eset, 1, function(x) {
              t.test(x[1:length(lev1)], x[(length(lev1) + 1):ncol(eset)], var.equal = TRUE, paired = paired)$p.value
            })
        }

        if (testmethod == "satterthwaite")
        {
            cat("Performing satterthwaite t.test...\n")
            testout = apply(eset, 1, function(x) {
                t.test(x[1:length(lev1)], x[(length(lev1) + 1):ncol(eset)],
                  paired = paired)$p.value
            })
        }

        adjustP = p.adjust(testout, method = Padj)
        difb = apply(eset, 1, function(x) {
            mean(x[1:length(lev1)]) - mean(x[(length(lev1) + 1):ncol(eset)])
            })
        out = cbind(testout, adjustP, difb, rowMeans(eset[, 1:length(lev1)]), rowMeans(eset[, (length(lev1) + 1):ncol(eset)]))
        rownames(out) = rownames(eset)
        colnames(out) = c("P-Value", "Adjust Pval", "beta-Difference",paste("Mean", paste(gcase, collapse = "_"), sep = "_"),paste("Mean", paste(gcontrol, collapse = "_"), sep = "_"))

       return(out)
  }

  pvalue_filter<-function(raw_result, rawpcut = 0.05, adjustpcut = 0.05, betadiffcut = 0.3)
  {
    if(is.null(rawpcut)&is.null(adjustpcut)&is.null(betadiffcut))
    {
      cat("All of them are retained")
      out= raw_result
    }else
    {
      if (is.null(rawpcut))
      {
        rawpcutout=raw_result[,"P-Value"]<=1
      }else
      {
        rawpcutout=raw_result[,"P-Value"]<=rawpcut
      }
      if (is.null(adjustpcut))
      {
       adjustpcutout=raw_result[,"Adjust Pval"]<=1
      }else
      {
        adjustpcutout=raw_result[,"Adjust Pval"]<=adjustpcut
      }
      if (is.null(betadiffcut))
      {
        betadiffcutout=abs(raw_result[,"beta-Difference"])<=1
      }else
      {
        betadiffcutout=abs(raw_result[,"beta-Difference"])>=betadiffcut
      }
      out = raw_result[rawpcutout & adjustpcutout & betadiffcutout,]
    }
    return(out)
  }

  linear<-function(eset,Padj = c("holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr", "none"),groupinfo)
  {
     cat("RUN linear regression\n")
     #require(MASS)
     testout = apply(eset, 1, function(x) { temp = summary(lm(x ~ as.numeric(as.character(groupinfo[,2])))) 
         pvalue = temp$coefficients[2, c(1, 4)]     
         return(pvalue) } )
     adjustP = p.adjust(testout[2, ], method = Padj) 
     out = cbind(testout[2, ], adjustP, testout[1, ])
     rownames(out) = rownames(eset)
     colnames(out) = c("P-Value", "Adjust Pval", "Coefficient")
     return(out)
  }

  BoxOfdm <- function(data,DM,contin=c("ON","OFF"),num,groupinfo,region,gcase="g1",gcontrol="g2")
  {
    if (nrow(DM)<num)
    {
      cat("The difference methylation number is less than ",num,"\n")
      num=nrow(DM)
    }
    if (num>0)
    {
      DM=DM[order(DM[,2]),]
      name=rownames(DM)[1:num]
      data=data[name,]
      if (contin=="ON")
      {
        groupinfo=groupinfo[order(groupinfo[,2]),]
        x=groupinfo[,2]
        if(num==1)
        {
          y=data[groupinfo[,1]]
          z=lm(y~x)
          plot(x,y,main=name,ylab="beta",sub=region)
          lines(x,fitted(z))
        }else
        {
          for (i in 1:num)
          {
            y=data[name[i],][groupinfo[,1]]
            z=lm(y~x)
            plot(x,y,main=name[i],ylab="beta",sub=region)
            lines(x,fitted(z))
          }
        }
      }
      else
      {
        sampGroups = unique(groupinfo[,2])
        c1 = groupinfo[, 2] %in% gcase
        c2 = groupinfo[, 2] %in% gcontrol
        if(num==1)
        {
          con=data[c1]
          trt=data[c2]
          eval(parse(text=paste("eset=list(",gcase,"=con,",gcontrol,"=trt)",sep="")))
          boxplot(eset ,ylab="beta",main=name,sub=region,boxwex=0.2,col=c("chocolate2","cyan3"))
          stripchart(eset, method = "jitter", add = TRUE,vertical=T,pch=20) 
        }else
        {
          for (i in 1:num)
          {
            con=data[name[i],][c1]
            trt=data[name[i],][c2]
            eval(parse(text=paste("eset=list(",gcase,"=con,",gcontrol,"=trt)",sep="")))
            boxplot(eset,ylab="beta",main=name[i],sub=region,boxwex=0.2,col=c("chocolate2","cyan3"))
            stripchart(eset, method = "jitter", add = TRUE,vertical=T,pch=20)
          }
        }      
      }
    }
    else
    {
      cat("there are no difference methylation to plot\n")
    }
  }

  eset = data@bmatrix
  groupinfo = data@groupinfo
  
  ####define global variable
  dmrAnno <- 0
  samples <- 0
  
  annotation_path <- system.file("extdata", "annotation",package="LncDM")
  load(paste(annotation_path, '/dmrAnno.Rdata', sep=''))
  if (classes=="gene")
    classes="protein_coding" 
  annotation = data@annot
  transAnno=dmrAnno[dmrAnno[,"transcriptType"]==classes,]
  cgList=unique(unlist(strsplit(as.character(transAnno[,4]), ",")))
  cgList=intersect(cgList,rownames(annotation))
  eset=eset[cgList,]
  #test all sites
  if (contin=="ON")
  {
    dms=linear(eset=eset ,Padj =Padj , groupinfo= groupinfo)
  }else
  {
    dms=test(eset=eset,testmethod =testmethod ,Padj =Padj ,groupinfo= groupinfo,gcase = gcase,gcontrol = gcontrol,paired = paired)
    eset=eset[rownames(dms),]
  }
  cgList=rownames(dms)

  siteINregion=strsplit(as.character(transAnno[,4]), ",")
  siteINregion=lapply(siteINregion,FUN=function(x){x[x %in% cgList]})
  siteNum=unlist(lapply(siteINregion,FUN=function(x){length(x)}))
  siteINregion=siteINregion[siteNum>1]
  transAnno=transAnno[siteNum>1,]
  transAnno[,4]=unlist(lapply(siteINregion,FUN=function(x){paste(x,collapse=",")}))
  transAnno=as.matrix(transAnno)

  region1<-function(x)
  {
    name=unlist(strsplit(x[4],","))
    name=name[order(annotation[name,"Coordinate_38"])]
    sub_test=dms[name,]

    index=rep(0,length(sub_test[,1]))
    index[sub_test[,"Adjust Pval"]>rawpcut]=0
    index[sub_test[,"Adjust Pval"]<=rawpcut & sub_test[,3]<0]=-1
    index[sub_test[,"Adjust Pval"]<=rawpcut & sub_test[,3]>0]=1
    lengths=rle(index)$lengths
    if(all(lengths<2))
    {
      region_beta=matrix(,0,0)
    }else
    {
      lengths=cumsum(rle(index)$lengths)
      values=rle(index)$values
      temp=intersect(which(values!=0),which(rle(index)$lengths>1))
      regionSites=lapply(temp,FUN = function(x){if(x==1){name[1:lengths[x]]} else{name[(lengths[x-1]+1):lengths[x]]}})
      regionName=lapply(regionSites,FUN = function(y){
                 #start sites
                 y=unlist(y)
                 if (min(which(name %in% y))==1)
                 {
                   star=max(as.numeric(x[2]),as.numeric(annotation[y[1],"Coordinate_38"])-2000)
                 }else
                 {
                   if((as.numeric(annotation[y[1],"Coordinate_38"])-2000)<as.numeric(annotation[name[min(which(name %in% y))-1],"Coordinate_38"]))
                   {
                     star=floor(mean(c(as.numeric(annotation[y[1],"Coordinate_38"]),as.numeric(annotation[name[min(which(name %in% y))-1],"Coordinate_38"]))))
                   }else
                   {
                     star=as.numeric(annotation[y[1],"Coordinate_38"])-2000
                   }
                 }
                 #end sites
                 if (max(which(name %in% y))==length(name))
                 {
                   sto=min(as.numeric(x[3]),as.numeric(annotation[y[length(y)],"Coordinate_38"])+2000)
                 }else
                 {
                   if((as.numeric(annotation[y[length(y)],"Coordinate_38"])+2000)>as.numeric(annotation[name[max(which(name %in% y))+1],"Coordinate_38"]))
                   {
                     sto=floor(mean(c(as.numeric(annotation[y[length(y)],"Coordinate_38"]),as.numeric(annotation[name[max(which(name %in% y))+1],"Coordinate_38"]))))
                   }else
                   {
                     sto=as.numeric(annotation[y[length(y)],"Coordinate_38"])+2000
                   }
                 }
                 paste(paste(x[10],x[6],sep="_"),paste(paste(x[1],star,sep=":"),sto,sep="-"),sep=":")
      })
      names(regionSites)=regionName
      region_beta=lapply(regionSites,FUN = function(x){colMeans(eset[x,])})
      region_beta=as.matrix(t(as.data.frame(region_beta)))
      if (length(region_beta)!=0)
      {
        temp=matrix(rep(x,nrow(region_beta)),nrow=nrow(region_beta),byrow=TRUE)
        region_beta=cbind(temp,region_beta)
      }
      rownames(region_beta)=regionName
    }
    region_beta
  }
  regionList=apply(transAnno,1,FUN=region1)
  
  regionListNum=lapply(regionList,FUN=function(x){dim(x)[1]})
  regionList=regionList[unlist(regionListNum)>0]
  
  if(length(regionList)>0)
  {
  regionList=lapply(regionList,FUN=function(x){t(x)})
  regionEset=t(as.matrix(as.data.frame(regionList)))
  regionListName=unlist(lapply(regionList,FUN=function(x){colnames(x)}))
  rownames(regionEset)=regionListName
  
  transAnno=regionEset[,1:ncol(transAnno)]
  regionEset=regionEset[,(ncol(transAnno)+1):ncol(regionEset)]
  temp=matrix(as.numeric(regionEset),nrow=nrow(regionEset))
  colnames(temp)=colnames(eset)
  rownames(temp)=rownames(regionEset)
  colnames(transAnno)=colnames(dmrAnno)
  regionEset=temp
  save(regionEset,transAnno,file="regionEset.Rdata")
  if(contin=="ON")
  {
    dmr=linear(eset=regionEset,Padj =Padj , groupinfo= groupinfo)
    dmr=dmr[dmr[,1]<rawpcut & dmr[,2]<adjustpcut,]
  }else
  {
    dmr=test(eset=regionEset,testmethod =testmethod ,Padj =Padj , groupinfo= groupinfo,gcase = gcase,gcontrol = gcontrol,paired = paired)
    dmr=pvalue_filter(dmr,rawpcut = rawpcut, adjustpcut = adjustpcut, betadiffcut = betadiffcut)
  }

  if (sole & length(rownames(dmr))!=0)
  {
    name=rownames(dmr)
    loca=unlist(lapply(name,FUN=function(x){paste(strsplit(x,":")[[1]][2],strsplit(x,":")[[1]][3],sep=":")}))
    temp=unique(loca)
    temp1=table(loca)
    temp=temp[temp1[temp]==1]
    name=name[which(loca %in% temp)]
    dmr=dmr[name,]    
    OUT=dmr
  }else
  {
    OUT=dmr
  }
  
  if(length(dim(OUT))!=0)
  {
  cat("plot boxplot of difference methylation\n")
  pdf("dmr_boxplot.pdf")
    BoxOfdm(data=regionEset,DM=dmr,contin=contin,num=num,groupinfo=groupinfo,region="DMR",gcase = gcase, gcontrol = gcontrol)
  dev.off()

  pdf("heatmap of DMR.pdf")
  #require(gplots)
  cat("plot heatmap of DMR\n")
  if (nrow(dmr)>0)
  {
    heatmapData=regionEset[rownames(dmr),]
    col=c(rgb(0,191,0,maxColorValue =255),rgb(35,245,26,maxColorValue =255),rgb(145,249,12,maxColorValue =255),rgb(239,239,0,maxColorValue =255),rgb(249,237,167,maxColorValue =255),rgb(177,177,255,maxColorValue =255),rgb(141,141,255,maxColorValue =255),rgb(104,104,255,maxColorValue =255),rgb(66,66,255,maxColorValue =255),rgb(0,0,255,maxColorValue =255))
 
    if (contin=="OFF")
    {
      labcol=groupinfo[colnames(heatmapData),2]
     
      ColSideColors=rep("a",ncol(heatmapData))
      ColSideColors[labcol==gcase]= "chocolate2"
      ColSideColors[labcol==gcontrol]= "cyan3"

      hv <- heatmap.2(heatmapData,dendrogram="both",col = col, trace="none", margins = c(5,5),xlab = "Samples", main = "heatmap of DMR",labCol=FALSE,ColSideColors=ColSideColors,keysize=1.5,density.info="density",key.title = "color key",key.xlab = "beta",cexRow=0.2,adjRow = c(0.4,0))
      legend("topright",legend=c(gcase,gcontrol),fill=c("chocolate2","cyan3"),bty="n")
    }else
    {
      ColSideColors=rep("white",ncol(heatmapData))
      Col.col=rainbow(length(unique(samples)))
      sampGroups=as.factor(samples)
      for(i in 1:length(ColSideColors))
      {
        ColSideColors[i]=Col.col[sampGroups[i]]
      }
      hv <- heatmap.2(heatmapData,dendrogram="both",col = col, trace="none", margins = c(5,5),xlab = "Samples", main = "heatmap of DMR",labCol=FALSE,ColSideColors=ColSideColors,keysize=1.5,density.info="density",key.title = "color key",key.xlab = "beta",cexRow=0.2,adjRow = c(0.4,0))
      legend("topright",legend=levels(sampGroups),fill=Col.col,bty="n") 
    }
  }else
  {
    cat("There are no dmr.\n")
  }
  dev.off()
  OUT=cbind(OUT,transAnno[rownames(dmr),])
  write.table(OUT,file="dmr.txt",col.names=T,row.names=T,sep="\t",quote=F)
  write.table(regionEset[rownames(dmr),] ,file="dmr_betaMatrix.txt",col.names=T,row.names=T,sep="\t",quote=F)
  }else
  {
    cat("There are no unique DMR.")
  }
  }else
  {
    cat("There are no DMR.")
  }
}
