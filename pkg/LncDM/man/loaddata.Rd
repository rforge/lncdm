\name{loaddata}
\alias{loaddata}
\alias{LincMethy450}
\title{
Load and preprocessed raw data
}
\description{
Read the file of signal_intensities, calculate the beta value, filter the unqualified samples and sites. Plot the heat map, box plot, density plot and density bean plot of CpG sites, and bar plot for detect P-value of samples.
}
\usage{
loaddata(fileDir,is_beta=FALSE,beta_method=c("M/(M+U)","M/(M+U+100)"),groupfile
,samplefilter = FALSE,contin=c("ON","OFF"),samplefilterperc = 0.75, XYchrom = 
c(FALSE, "X","Y", c("X", "Y")),sitefilter = FALSE, sitefilterperc = 0.75,
 filterDecetP =0.05, normalization  = FALSE,transfm = c(FALSE, "arcsinsqr", "logit")
,snpfilter=c(FALSE,"within_10","prob_snp"),gcase="case",gcontrol ="control",skip=0
,imputation=c("mean","min","knn"),knn.k=10)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{fileDir}{ The folder name of samples' signal_intensities files. }
  \item{is_beta}{ Logical. The signal_intensities is beta value or not. }
  \item{beta_method}{ The method for calculating the beta. }
  \item{groupfile}{ The name of phenotype file. }
  \item{samplefilter}{ Logical. Filter the samples whose most detection P values aren't significative or not. }
  \item{contin}{ 'ON' means the phenotype is continuous,just like age etc. 'OFF' means it is discontinued. }
  \item{samplefilterperc}{ A number in [0,1]. The samples whose percent of the significative detection P values less than this number will be filtered. }
  \item{XYchrom}{ The CpG sites in X or Y chromosome should or shuoldn't be filtered. }
  \item{sitefilter}{ Logical. Filter the sites whose most detection P values aren't significative or not. }
  \item{sitefilterperc}{ A number in [0,1]. The sites whose percent of the significative detection P values less than this number will be filtered. }
  \item{filterDecetP}{ Threshold: value of significative detection P. Always 0.05 or 0.01. }
  \item{normalization}{ Logical. Normalization for the different chips or not. }
  \item{transfm}{ Data transformation for beta or not. Contains 'arcsinsqr' and 'logit'. }
  \item{snpfilter}{ The CpG sites that contain SNP sites with 10bp or 50bp shuold or shouldn't be filetered. }
  \item{gcase}{ The name of case group while contin is 'OFF'. }
  \item{gcontrol}{ The name of case group while contin is 'OFF'. }
  \item{skip}{ integer: the number of lines of the data file to skip before beginning to read signal_intensities data, the first row must be signal values. }
  \item{imputation}{ The method to fill the NA.Contains 'mean', 'min' and 'knn'. }
  \item{knn.k}{ The K number if imputation is 'knn'. }
}
\details{ 
Loaddata is designed to load and process the methylated data for the package. It provides two methods to calculate the beta value,which means the ratio of methylation,M/(M+U) and M/(M+U+100),M means the intensity of methylation and U means the intensity of unmethylation.
For the methylated data, a file per sample. In the signal_intensities file, there are four columns, CpG ID, Methylated_Intensity, Unmethylated_Intensity and Detection_P_value.
The groupfile that explain the phenotype of samples. Distinguish the case or control. The samples that at the same group have the same label. The sample IDs are same as the names of corresponding signal_intensities file (without File Suffixes).   
Loaddata also call the other function to plot the heat map,box plot, density plot and density bean plot of CpG sites, and bar plot for detect P-value of samples.  
}
\value{
Loaddata will return an object of class LincMethy450. And return some plots to describe the information of data.
}
\author{
Hui Zhi\email{zhihui013201@gmail.com},Yanxun Su\email{hmu_yanxunsu@163.com},Xin Li\email{lixin920126@163.com}
}

\seealso{
 See Also  \code{\link{dms}}, \code{\link{dme}} and \code{\link{dmr}}
}
\examples{
\dontrun{
 ##the directory of phenotype and 450k methylation's sample data
 Dir <- system.file("extdata/localdata",package="LncDM")
 setwd(Dir)
 ###phenotype file's name
 groupfile <- "BRCA_pheno.txt"
 ###our methylation data in the subdirectory "Level_2" is just example data, when you 
 ###run this function, please prepare complete sample files, and change default directory 
 ###to yourself
 loadData <- loaddata(fileDir="Level_2",is_beta=FALSE,beta_method="M/(M+U)",groupfile=groupfile,
 samplefilter = TRUE,contin="OFF",samplefilterperc = 0.75,XYchrom = c(FALSE, "X","Y"),sitefilter = TRUE, 
 sitefilterperc = 0.75,filterDecetP=0.05,normalization  = FALSE,transfm = FALSE,snpfilter=c(FALSE,"prob_snp"),
 gcase="case",gcontrol="control",skip=2,imputation="knn",knn.k=10)
 ###save the loadData in order to caculate dms,dmr and dme
 save(loadData,file="loadData.Rdata",compress="xz")
}
}
