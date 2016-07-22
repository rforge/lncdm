dme <-
function(data,classes=c("lincRNA","gene","processed_transcript","pseudogene"),contin=c("ON","OFF"),testmethod = c("wilcox","limma", "t.test", "satterthwaite"), Padj = c("holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr", "none"), gcase = "case", gcontrol = "control", paired = FALSE,rawpcut = 0.05, adjustpcut = 0.05, betadiffcut = 0.14,num)
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

  groupinfo = data@groupinfo
  annotation = data@annotation
  transAnno = data@transAnno
  transBeta = data@transBeta
  rownames(transAnno)=paste(transAnno[,10],paste(transAnno[,1],transAnno[,2],transAnno[,3],sep="-"),sep=":")
  rownames(transBeta)=paste(transAnno[,10],paste(transAnno[,1],transAnno[,2],transAnno[,3],sep="-"),sep=":")
  
  ###define global variable
  DM1 <- 0
  p <- 0
  heatmaData <- 0
  if(classes == "gene" )
    classes = "protein_coding"

  RegionName=unique(transAnno[,6])

  #plot boxplots.
  eval(parse(text=paste("pdf(\"",classes,"_plot.pdf\")",sep="")))
  if(contin=="ON")
  {
    #require(gplots)
    zero=c(length(RegionName)+1)
    for(i in 1:length(RegionName))
    {
      eset=transBeta[transAnno[,6]==RegionName[i],]
      
      titl=RegionName[i]
      if(titl=="5'UTR")
        titl="5UTR"
      if(titl=="3'UTR")
        titl="3UTR"

      if (sum(transAnno[,6]==RegionName[i])>1)
      {
        eval(parse(text = paste("a",titl,"_test = as.data.frame(linear(eset,Padj=Padj,groupinfo=groupinfo))",sep="")))
        eval(parse(text=paste("DM1=a",titl,"_test",sep="")))
        if (length(DM1)<2)
        {
          zero=c(zero,i)
        }
        eval(parse(text = paste("BoxOfdm(data=eset,DM=","a", titl, "_test,contin=contin,num=num,groupinfo=groupinfo,region=\"",RegionName[i],"\")",sep = "")))
        eval(parse(text = paste("a",titl,"_test=","a",titl,"_test[","a",titl,"_test[,1]<rawpcut & ","a",titl,"_test[,2]<adjustpcut,]",sep="")))
        eval(parse(text=paste("DME_",titl,"_beta=as.data.frame(eset[rownames(","a",titl,"_test),])",sep="")))
        eval(parse(text=paste("DME_",titl,"_anno=as.data.frame(transAnno[rownames(","a",titl,"_test),])",sep="")))     
        eval(parse(text = paste("if (nrow(","a",titl,"_test)>1){p=1}else{p=0}",sep="")))
        if (p==1)
        {
          samples=groupinfo[, 2]
          eval(parse(text=paste("heatmaData=eset[rownames(","a",titl,"_test[order(","a",titl,"_test[,\"Adjust Pval\"])[1:min(1000,nrow(","a",titl,"_test))],]),]",sep="")))
          col=c(rgb(0,191,0,maxColorValue =255),rgb(35,245,26,maxColorValue =255),rgb(145,249,12,maxColorValue =255),rgb(239,239,0,maxColorValue =255),rgb(249,237,167,maxColorValue =255),rgb(177,177,255,maxColorValue =255),rgb(141,141,255,maxColorValue =255),rgb(104,104,255,maxColorValue =255),rgb(66,66,255,maxColorValue =255),rgb(0,0,255,maxColorValue =255))
         
          ColSideColors=rep("white",ncol(heatmapData))
          Col.col=rainbow(length(unique(samples)))
          sampGroups=as.factor(samples)
          for(i in 1:length(ColSideColors))
          {
            ColSideColors[i]=Col.col[sampGroups[i]]
          }  
          main=paste("heatmap of ",RegionName[i],sep="")
          hv <- heatmap.2(heatmapData,dendrogram="both",col = col, trace="none", margins = c(5,5),xlab = "Samples", main = main,labCol=FALSE,ColSideColors=ColSideColors,keysize=1.5,density.info="density",key.title = "color key",key.xlab = "beta",cexRow=0.2,adjRow = c(0.4,0)) 
          legend("topright",legend=levels(sampGroups),fill=Col.col,bty="n") 
        }
      }
      else
      {
        zero=c(zero,i)
      }
    }
    RegionName=RegionName[-zero]
  }

  else
  {   
    #require(gplots)
    zero=c(length(RegionName)+1)
    for(i  in  1:length(RegionName))
    {
      cat("calculating", RegionName[i], "\n")
      eset=transBeta[transAnno[,6]==RegionName[i],]
      titl=RegionName[i]
      if(titl=="5'UTR")
        titl="5UTR"
      if(titl=="3'UTR")
        titl="3UTR"
        
      if(sum(transAnno[,6]==RegionName[i])>1)
      {
        eval(parse(text = paste("a",titl,"_test=as.data.frame(test(eset,testmethod=testmethod,Padj=Padj,groupinfo = groupinfo,paired = paired,gcase = gcase,gcontrol = gcontrol))", sep = "")))
        eval(parse(text = paste("a",titl,"_test=pvalue_filter(","a",titl,"_test,rawpcut = rawpcut, adjustpcut = adjustpcut, betadiffcut = betadiffcut)", sep = "")))
        eval(parse(text=paste("DME_",titl,"_beta=as.data.frame(eset[rownames(","a",titl,"_test),])",sep="")))    
        eval(parse(text=paste("DME_",titl,"_anno=as.data.frame(transAnno[rownames(","a",titl,"_test),])",sep="")))       
        eval(parse(text=paste("DM1=a",titl,"_test",sep="")))
        if (length(rownames(DM1))<1)
        {
          zero=c(zero,i)
        }
        eval(parse(text = paste("BoxOfdm(data=eset,DM=", "a",titl, "_test,contin=contin,num=num,groupinfo=groupinfo,region=\"",RegionName[i],"\",gcase = gcase, gcontrol = gcontrol)",sep = "")))
        eval(parse(text = paste("if (nrow(","a",titl,"_test)>1){p=1}else{p=0}",sep="")))
        if (p==1)
        {
          samples=groupinfo[, 2]
          eval(parse(text=paste("heatmaData=eset[rownames(","a",titl,"_test[order(","a",titl,"_test[,\"Adjust Pval\"])[1:min(1000,nrow(","a",titl,"_test))],]),]",sep="")))
          col=c(rgb(0,191,0,maxColorValue =255),rgb(35,245,26,maxColorValue =255),rgb(145,249,12,maxColorValue =255),rgb(239,239,0,maxColorValue =255),rgb(249,237,167,maxColorValue =255),rgb(177,177,255,maxColorValue =255),rgb(141,141,255,maxColorValue =255),rgb(104,104,255,maxColorValue =255),rgb(66,66,255,maxColorValue =255),rgb(0,0,255,maxColorValue =255))
         
          labcol=c(samples[which(samples==gcase)],samples[which(samples==gcontrol)])
          heatmapData=cbind(heatmaData[,which(samples==gcase)],heatmaData[,which(samples==gcontrol)])
          ColSideColors=rep("a",ncol(heatmapData))
          ColSideColors[labcol==gcase]= "chocolate2"
          ColSideColors[labcol==gcontrol]= "cyan3"

          main=paste("heatmap of ",RegionName[i],sep="")
          hv <- heatmap.2(heatmapData,dendrogram="both",col = col, trace="none", margins = c(5,5),xlab = "Samples", main = main,labCol=FALSE,ColSideColors=ColSideColors,keysize=1.5,density.info="density",key.title = "color key",key.xlab = "beta",cexRow=0.2,adjRow = c(0.4,0))
          legend("topright",legend=c(gcase,gcontrol),fill=c("chocolate2","cyan3"),bty="n")
        }
      }else
      {
        zero=c(zero,i)
      }
    }
    RegionName=RegionName[-zero]
  }
  dev.off()
#####################################return result###################################
  if (length(RegionName)==0)
  {
    cat("There are no DME.")
  }else
  {
  for(i in 1:length(RegionName))
  {
    titl=RegionName[i]
    if(titl=="5'UTR")
      titl="5UTR"
    if(titl=="3'UTR")
      titl="3UTR"
    eval(parse(text = paste("if (nrow(","a",titl,"_test)>1){p=1}else{p=0}",sep="")))
    if(p==0)
    {
      eval(parse(text = paste("a",titl,"_out=cbind(t(DME_",titl,"_anno),","a",titl,"_test)",sep = "")))
      eval(parse(text = paste("a",titl,"_out[6]=RegionName[i]",sep = ""))) 
    }else
    {
      eval(parse(text = paste("a",titl,"_out=cbind(DME_",titl,"_anno,","a",titl,"_test)",sep = "")))
      eval(parse(text = paste("a",titl,"_out[,6]=RegionName[i]",sep = "")))
    }
       
  }
  #require(WriteXLS)
  RegionName[RegionName=="5'UTR"]="5UTR"
  RegionName[RegionName=="3'UTR"]="3UTR"
  eval(parse(text = paste("WriteXLS(paste(\"a\",RegionName, \"_out\", sep = \"\"), ExcelFileName = \"",classes,"_DME.xls\", SheetNames = RegionName, row.names = TRUE)",sep = "")))
  
  dir.create("diffPvalue")
  for(i in 1:length(RegionName))
  {
    titl=RegionName[i]
    if(titl=="5'UTR")
      titl="5UTR"
    if(titl=="3'UTR")
      titl="3UTR"
    eval(parse(text = paste("write.table(a",titl,"_out,file=\"diffPvalue/DME_",titl,"_diffPvalue.txt\",col.names=T,row.names=T,sep=\"\t\",quote=F)",sep = "")))
  }
    
  dir.create("betaMatrix")
  for(i in 1:length(RegionName))
  {
    titl=RegionName[i]
    if(titl=="5'UTR")
      titl="5UTR"
    if(titl=="3'UTR")
      titl="3UTR"
    eval(parse(text = paste("write.table(DME_",titl,"_beta,file=\"betaMatrix/DME_",titl,"_betaMatrix.txt\",col.names=T,row.names=T,sep=\"\t\",quote=F)",sep = "")))
  }
  }
}
