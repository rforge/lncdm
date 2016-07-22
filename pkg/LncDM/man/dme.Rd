\name{dme}
\alias{dme}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Find the differential methylation elements of lincRNA, protein coding gene ,processed transcript and pseudogene
}
\description{
Find the differential methylation elements or the elements that related with the phenotype. The elements are belong to lincRNA, protein coding gene, processed transcript and pseudogene.
}
\usage{
dme(data,classes=c("lincRNA","gene","processed_transcript","pseudogene"),contin=c(
"ON","OFF"),testmethod = c("wilcox","limma", "t.test", "satterthwaite"), Padj = c(
"holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr", "none"), gcase = 
"case", gcontrol = "control", paired = FALSE,rawpcut = 0.05, adjustpcut = 0.05,
 betadiffcut = 0.14,num)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{ The objects of class LincMethy450 which return from \code{\link{loaddata}}. The beta matrixs of sites. A site per row and A sample per column. }
  \item{classes}{ Whose element will be calculated. }
  \item{contin}{ If phenotype is continuous,contin is 'ON',use linear regression to find the elements that related with the phenotype. }
  \item{testmethod}{ The method to do the test to find dme while contin is 'OFF' which means phenotype is discontinuous.  }
  \item{Padj}{ The method of multiple testing adjustment to adjust P value. }
  \item{gcase}{ The name of case group while contin is 'OFF'. }
  \item{gcontrol}{ The name of case group while contin is 'OFF'. }
  \item{paired}{ Whether compare in pairs while do t.test. }
  \item{rawpcut}{ It is the threshold for cutting raw P value. }
  \item{adjustpcut}{ It is the threshold for cutting adjust P value. }
  \item{betadiffcut}{ The minimum differential between two groups' means while contin is 'OFF'.   }
  \item{num}{ The number which is the parameter of of elements to plot. }
}
\details{
dme is designed to find differential methylation elements or the transcripts' elements that related with the continuous phenotype. If contin is 'ON', it means the phenotype is continuous, and linear regression will be used. If the phenotype isn't continuous, test such as t test will be used.
}
\value{
dme will return two excel files that one contains the transcripts' elements whose P value less than rawpcut, adjust P less than adjustpcut and the differ of the means of two groups more than betadiffer, while another is the beta matrix of these significant elements. There are box plot for most significative elements and heat map all significative elements.
}
\author{
Hui Zhi\email{zhihui013201@gmail.com},Yanxun Su\email{hmu_yanxunsu@163.com},Xin Li\email{lixin920126@163.com}
}

\seealso{
 See Also  \code{\link{dms}} and  \code{\link{dmr}}
}
\examples{
  Dir <- system.file("extdata/localdata",package="LncDM")
  dir.create(paste(Dir,"/dme",sep=""))
  ###user can set the dir of their own
  setwd(paste(Dir,"/dme",sep=""))
  ###load the result of loaddata()
  load(paste(Dir,"/loadData.Rdata",sep=""))
  ###dme is based on the result of the regionLevel()
  Region <- regionLevel(data=loadData,indexmethod = "mean",classes="lincRNA")
  dme(data=Region,classes="lincRNA",contin="OFF",testmethod = "t.test", Padj = "fdr", 
  gcase = "case", gcontrol = "control", paired = FALSE,rawpcut = 0.05, adjustpcut = 0.05, 
  betadiffcut = 0.3,num=1)
}


