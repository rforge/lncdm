\name{dmr}
\alias{dmr}
\title{
Find the differential methylation regions in the elements of lincRNA, protein coding gene ,processed transcript and pseudogene
}
\description{
Find the differential methylation regions or the regions that related with the phenotype who are in the elements of lincRNA etc.
}
\usage{
dmr(data,contin=c("ON","OFF"),classes=c("lincRNA","gene","processed_transcript",
"pseudogene"),testmethod = c("wilcox","limma", "t.test", "satterthwaite"), Padj =
c("holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr", "none"), gcase =
"case", gcontrol = "control", paired = FALSE,rawpcut = 0.05, adjustpcut = 0.05,
 betadiffcut = 0.3,num,sole=FALSE)
}
\arguments{
  \item{data}{ The objects of class LincMethy450 which return from \code{\link{loaddata}}. The beta matrixs of sites. A site per row and A sample per column. }
  \item{classes}{ Whose regions will be calculated. }
  \item{contin}{ If phenotype is continuous,contin is 'ON',and use linear regression to find the regions that related with the phenotype. }
  \item{testmethod}{ The method to do the test to find dmr while contin is 'OFF' which means phenotype is discontinuous.  }
  \item{Padj}{ The method of multiple testing adjustment to adjust P value. }
  \item{gcase}{ The name of case group while contin is 'OFF'. }
  \item{gcontrol}{ The name of case group while contin is 'OFF'. }
  \item{paired}{ Whether compare in pairs while do t.test. }
  \item{rawpcut}{ It is the threshold for cutting raw P value. }
  \item{adjustpcut}{ It is the threshold for cutting adjust P value. }
  \item{betadiffcut}{ The minimum differential between two groups' means while contin is 'OFF'.   }
  \item{num}{ The number which is the box plot of different methylation regions to plot. }
  \item{sole}{ Whether keep the no duplicate records. }
}
\details{
dmr is designed to find differential methylation regions or the regions in the elements of lincRNA etc. that related with the continuous phenotype. If contin is 'ON', it means the phenotype is continuous, and linear regression will be used. If the phenotype isn't continuous, test such as t test will be used.
}
\value{
dmr will return two txt files that one contains the transcripts' elements and the region of chromosome whose P value less than rawpcut, adjust P less than adjustpcut and the differ of the means of two groups more than betadiffer, while another is the beta matrix of these significant regions. There are box plot for most significative regions and heat map all significative regions.
}
\author{
Hui Zhi\email{zhihui013201@gmail.com},Yanxun Su\email{hmu_yanxunsu@163.com},Xin Li\email{lixin920126@163.com}
}

\seealso{
 See Also  \code{\link{dms}} and \code{\link{dme}}
}
\examples{
   Dir <- system.file("extdata/localdata",package="LncDM")
   dir.create(paste(Dir,"/dmr",sep=""))
   ###user can set the dir of their own
   setwd(paste(Dir,"/dmr",sep=""))
   ###load the result of loaddata()
   load(paste(Dir,"/loadData.Rdata",sep=""))
   dmr(data=loadData,contin="OFF",classes="lincRNA",testmethod = "t.test", Padj = "fdr", 
   gcase = "case", gcontrol = "control", paired = FALSE,rawpcut = 0.05, adjustpcut = 0.05, 
   betadiffcut = 0.3,num=1,sole=FALSE)
}


