\name{regionLevel}
\alias{regionLevel}
\alias{RegionMethy450}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Calculate the beta values of elements in lincRNA etc
}
\description{
Calculate the beta values of elements in lincRNA etc. based on reannotation information.

LincRNA and processed transcript's elements are TSS1500, TSS200, 1_exon, genebody and intron. Protein coding gene and pseudogene transcript's elements are 5'UTR, 3'UTR, TSS1500, TSS200, 1_exon, genebody and intron.
}
\usage{
regionLevel(data,indexmethod = c("mean", "median"),classes=c("gene","lincRNA",
"processed_transcript","pseudogene"))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{ The objects of class LincMethy450 which return from \code{\link{loaddata}}. The beta matrix of sites. A site per row and A sample per column. }
  \item{indexmethod}{ The method to calculate the beta if there are more than two CpG sites in a element. }
  \item{classes}{ Whose CpG sites will be calculated. }
}
\details{
The function, regionLevel, is designed to calculate the beta value of gene, lincRNA, processed_transcript and pseudogenes' elements. If indexmethod is "mean" ,and the mean of CpG sites' beta value will be regarded as the element's beta. If classes is "gene", only the elements of protein coding genes' transcript will be calculated. If want to calculate more than one class, please do it respectively.
}
\value{
This methods will return a object of class RegionMethy450. It contains some matrix of elements' beta values. A row per transcript, a col per sample.
}
\author{
Hui Zhi\email{zhihui013201@gmail.com},Yanxun Su\email{hmu_yanxunsu@163.com},Xin Li\email{lixin920126@163.com}
}

\seealso{
 See Also  \code{\link{dms}}, \code{\link{dme}} and \code{\link{dmr}}
}
\examples{
\dontrun{
  Dir <- system.file("extdata/localdata",package="LncDM")
  setwd(Dir)
  ###load the result of loaddata()
  #load(paste(Dir,"/loadData.Rdata",sep=""))
  Region <- regionLevel(data=loadData,indexmethod = "mean",classes="lincRNA")
  ###save the region data in order to caculate dme
  save(Region,file="Region.Rdata")
}
}


