\name{dms}
\alias{dms}
\alias{mhtObject}
\title{
Find the differential methylation CpG sites in the lincRNA, protein coding gene ,processed transcript and pseudogene
}
\description{
Find the differential methylation sites or the sites that related with the phenotype who are in the lincRNA etc.
}
\usage{
dms(data,contin=c("ON","OFF"),classes=c("lincRNA","gene","processed_transcript",
"pseudogene"),testmethod = c("wilcox","limma", "t.test", "satterthwaite"), Padj =
c("holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr", "none"), gcase =
"case", gcontrol = "control", paired = FALSE,rawpcut = 0.05, adjustpcut = 0.05,
 betadiffcut = 0.3,XY=c(FALSE,"X","Y",c("X","Y")),tlog=FALSE,num)
}
\arguments{
  \item{data}{ The objects of class LincMethy450 which return from \code{\link{loaddata}}. The beta matrixs of sites. A site per row and A sample per column. }
  \item{classes}{ Whose CoG sites will be calculated }
  \item{contin}{ If phenotype is continuous,contin is 'ON',and use linear regression to find the sites that related with the phenotype. }
  \item{testmethod}{ The method to do the test to find dms while contin is 'OFF' which means phenotype is discontinuous.  }
  \item{Padj}{ The method of multiple testing adjustment to adjust P value. }
  \item{gcase}{ The name of case group while contin is 'OFF'. }
  \item{gcontrol}{ The name of case group while contin is 'OFF'. }
  \item{paired}{ Whether compare in pairs while do t.test. }
  \item{rawpcut}{ The threshold for cutting raw P value. }
  \item{adjustpcut}{ The threshold for cutting adjust P value. }
  \item{betadiffcut}{ The minimum differential between two groups' means while contin is 'OFF'. }
  \item{num}{ The number of sites,elements or regions to plot. }
  \item{XY}{ The chromosomes that where CpG sites were filtered. }
  \item{tlog}{ Whether do the log transform for the P values. }
}
\details{
dms is designed to find differential methylation sites or the CpG sites that related with the continuous phenotype. If contin is 'ON', which means the phenotype is continuous, and linear regression will be used. If the phenotype isn't continuous, test such as t test will be used.
}
\value{
dms will return two txt files that one contains the sites whose P value less than rawpcut, adjust P less than adjustpcut and the differ of the means of two groups more than betadiffer while another is the beta matrix of dms. There are box plot for most significative sites,heat map and Manhattan of all significative sites.
}
\author{
Hui Zhi\email{zhihui013201@gmail.com},Yanxun Su\email{hmu_yanxunsu@163.com},Xin Li\email{lixin920126@163.com}
}

\seealso{
 See Also  \code{\link{dme}} and \code{\link{dmr}}
}
\examples{
   Dir <- system.file("extdata/localdata",package="LncDM")
   dir.create(paste(Dir,"/dms",sep=""))
   ###user can set the dir of their own
   setwd(paste(Dir,"/dms",sep=""))
   ###load the result of loaddata()
   load(paste(Dir,"/loadData.Rdata",sep=""))
   dms(data=loadData,contin="OFF",classes="lincRNA",testmethod = "t.test", Padj = "fdr", 
   gcase = "case", gcontrol = "control", paired = FALSE,rawpcut = 0.05, adjustpcut = 0.05, 
   betadiffcut = 0.3,XY=c(FALSE,"X","Y"),tlog=FALSE,num=1)
}


