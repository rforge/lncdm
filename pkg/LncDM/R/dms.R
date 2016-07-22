dms <-
function(data,contin=c("ON","OFF"),classes=c("lincRNA","gene","processed_transcript","pseudogene"),testmethod = c("wilcox","limma", "t.test", "satterthwaite"), Padj = c("holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr", "none"), gcase = "case", gcontrol = "control", paired = FALSE,rawpcut = 0.05, adjustpcut = 0.05, betadiffcut = 0.3,XY=c(FALSE,"X","Y",c("X","Y")),tlog=FALSE,num)
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

   mhtplot<-function(data,lab.col="brown",axis.col="black",bg.col="white",col=c("black","red","green","turquoise","tan3","hotpink"),tlog=FALSE,XY=c(FALSE,"X","Y",c("X","Y")),...)
  {
    if(class(data) != "mhtObject")
      stop("It is not a S4 class \n")

    if(tlog == TRUE)
    {
      tmp<-log10(data@pvalue)
      data@pvalue<--tmp
    }

    f<-levels(factor(data@chr))
    f=order(f)

    if (XY[1])
    {
      if (length(XY)>2)
      {
      }else
      {
        if(XY[2] == "X")
          f[23]="Y"
        if(XY[2] == "Y")
          f[23]="X"
      }
    }else
    {
      f[c(23,24)]=c("X","Y")
    }

    colv<-col
    col<-rep(colv,ceiling(length(f)/length(colv)))
    maxt<-length(f)*10

    par(bg=bg.col,col.lab=lab.col,col.axis=axis.col,xaxs="r")
    par(lab=c(5,5,7))

    if(tlog==TRUE)
    {
      plot(rep(maxt,length(data@pvalue)),data@pvalue,type="n",xlim=c(0,maxt/10),xlab="Base Position",ylab="Observed -log10(adjustP-value)",axes=FALSE,ylim=c(0,ceiling(max(data@pvalue,na.rm=TRUE))))
      axis(1,seq(.5,length(f)-0.5,1),labels=f,tick=FALSE)
      axis(1,seq(0,length(f),1),labels=FALSE)
      axis(2,at=c(0,ceiling(max(data@pvalue,na.rm=TRUE))),lwd=0.5)

      abline(h=2,lty=3,col="black")
      abline(h=1.3,col="black")
    }else
    {
      plot(rep(maxt,length(data@pvalue)),data@pvalue,type="n",xlim=c(0,maxt/10),xlab="Base Position",ylab="P-value",axes=FALSE,ylim=c(0,ceiling(max(data@pvalue,na.rm=TRUE))))
      axis(1,seq(.5,length(f)-0.5,1),labels=f,tick=FALSE)
      axis(1,seq(0,length(f),1),labels=FALSE)
      abline(h=max(data@pvalue,na.rm=TRUE)*0.95,lty=3,col="black")
      axis(2,at=c(0,1),lwd=0.5)
      abline(h=0.01,col="black",lty=3)
      abline(h=0.05,col="black")
    }

    for ( i in 1:length(f) )
    {
      k<-data@coor[data@chr==f[i]]
      k<-k[!is.na(k)]
	  if(length(k)!=0){
      n=1
      xt<-0
      for( j in 1:length(k) )
      {
        xt[n]<-(k[j]-min(k))*.8/(max(k)-min(k))+(i-1)+0.1
        n<-n+1
      }
      y=data@pvalue[data@chr==f[i]]
      y<-y[!is.na(y)]
      points(xt,y,col=col[i],type="p",pch=20,...)
	  }
    }
  }

   mhtObject<-function(data,annotation)
  {
     rownames(annotation)=annotation[,1]
     chr=annotation[rownames(data),"CHR_38"]
     coor=as.numeric(annotation[rownames(data),"Coordinate_38"])
     pvalue=data[,1]
     setClass("mhtObject", representation(chr="character",coor="numeric",pvalue="numeric"))
     out=new("mhtObject",chr=chr,coor=coor,pvalue=pvalue)
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
  annotation = data@annot
  lincRNA =  data@linc
  processed_transcript =  data@processed
  protein_coding =  data@protein
  pseudogene =  data@pseu

  if (classes=="gene")
    cpgID=unique(unlist(strsplit(as.character(protein_coding[,4]), ",")))
  else
    eval(parse(text=paste("cpgID=unique(unlist(strsplit(as.character(",classes,"[,4]), \",\")))",sep="")))
  eset=eset[cpgID,]
  
  pdf("Manhattan.pdf",width=12)
  if (contin=="ON")
  {
    site_test=linear(eset,Padj,groupinfo)
    #plot Manhattan plots
	rownames(annotation)=annotation[,1]
    chr=annotation[rownames(site_test),"CHR_38"]
    coor=as.numeric(annotation[rownames(site_test),"Coordinate_38"])
    pvalue=site_test[,1]
    setClass("mhtObject", representation(chr="character",coor="numeric",pvalue="numeric"),where = topenv(parent.frame()))
    mht=new("mhtObject",chr=chr,coor=coor,pvalue=pvalue)
    #mht=mhtObject(data=site_test,annotation=annotation)
    mhtplot(data=mht,lab.col="brown",axis.col="black",bg.col="white",col=c("black","red","green","turquoise","tan3","hotpink"),tlog=tlog,XY=XY)
    site_test=site_test[site_test[,1]<rawpcut & site_test[,2]<adjustpcut,]    
    DMS_beta=eset[rownames(site_test),]
  }else
  {
    site_test=as.data.frame(test(eset,testmethod=testmethod,Padj=Padj,groupinfo = groupinfo,paired = paired,gcase = gcase,gcontrol = gcontrol))
    cat("plot Manhattan\n")
	rownames(annotation)=annotation[,1]
    chr=annotation[rownames(site_test),"CHR_38"]
    coor=as.numeric(annotation[rownames(site_test),"Coordinate_38"])
    pvalue=site_test[,1]
    setClass("mhtObject", representation(chr="character",coor="numeric",pvalue="numeric"),where = topenv(parent.frame()))
    mht=new("mhtObject",chr=chr,coor=coor,pvalue=pvalue)
    #mht=mhtObject(data=site_test,annotation=annotation)
    mhtplot(data=mht,lab.col="brown",axis.col="black",bg.col="white",col=c("black","red","green","turquoise","tan3","hotpink"),tlog=tlog,XY=XY)
    site_test=pvalue_filter(site_test,rawpcut = rawpcut, adjustpcut = adjustpcut, betadiffcut = betadiffcut)
    DMS_beta=eset[rownames(site_test),]
  }    
  dev.off()


  if (nrow(site_test)>1)
  {
    pdf("heatmap of DMS.pdf")
    #require(gplots)
    cat("plot heatmap of DMS\n")
    samples=groupinfo[colnames(eset),2]
    heatmapData=eset[rownames(site_test),]
    col=c(rgb(0,191,0,maxColorValue =255),rgb(35,245,26,maxColorValue =255),rgb(145,249,12,maxColorValue =255),rgb(239,239,0,maxColorValue =255),rgb(249,237,167,maxColorValue =255),rgb(177,177,255,maxColorValue =255),rgb(141,141,255,maxColorValue =255),rgb(104,104,255,maxColorValue =255),rgb(66,66,255,maxColorValue =255),rgb(0,0,255,maxColorValue =255))
    if (contin=="OFF")
    {
      labCol=c(samples[which(samples==gcase)],samples[which(samples==gcontrol)])
      heatmapData=cbind(heatmapData[,which(samples==gcase)],heatmapData[,which(samples==gcontrol)])
      ColSideColors=rep("a",ncol(heatmapData))
      ColSideColors[labCol==gcase]= "chocolate2"
      ColSideColors[labCol==gcontrol]= "cyan3"
      hv <- heatmap.2(heatmapData,dendrogram="none",col = col, trace="none", margins = c(5,5),xlab = "Samples",  main = "heatmap of DMS",labCol=FALSE,ColSideColors=ColSideColors,keysize=1.5,density.info="density",key.title = "color key",key.xlab = "beta",cexRow=0.2)
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
      hv <- heatmap.2(heatmapData, dendrogram="none",col = col, trace="none", margins = c(5,5),xlab = "Samples",  main = "heatmap of DMS",labCol=FALSE,ColSideColors=ColSideColors,density.info="density",key.title = "color key",key.xlab = "beta",cexRow=0.2) 
      legend("topright",legend=levels(sampGroups),fill=Col.col,bty="n") 
    }
    dev.off()
  }else
  {
    cat("There are not enough dms to plot heatmap.\n")
  }
  
  pdf("dms_boxplot.pdf")
  cat("plot boxplot of difference methylation\n")
  BoxOfdm(data=eset,DM=site_test,contin=contin,num=num,groupinfo=groupinfo,region="site",gcase = gcase, gcontrol = gcontrol)
  dev.off()

  write.table(site_test,file="site_test.txt",col.names=T,row.names=T,sep="\t",quote=F)
  write.table(DMS_beta,file="DMS_betaMatrix.txt",col.names=T,row.names=T,sep="\t",quote=F)  
  return(site_test)
}
