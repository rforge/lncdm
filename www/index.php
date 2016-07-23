
<!-- This is the project specific website template -->
<!-- It can be changed as liked or replaced by other content -->

<?php

$domain=ereg_replace('[^\.]*\.(.*)$','\1',$_SERVER['HTTP_HOST']);
$group_name=ereg_replace('([^\.]*)\..*$','\1',$_SERVER['HTTP_HOST']);
$themeroot='r-forge.r-project.org/themes/rforge/';

echo '<?xml version="1.0" encoding="UTF-8"?>';
?>
<!DOCTYPE html
	PUBLIC "-//W3C//DTD XHTML 1.0 Transitional//EN"
	"http://www.w3.org/TR/xhtml1/DTD/xhtml1-transitional.dtd">
<html xmlns="http://www.w3.org/1999/xhtml" xml:lang="en" lang="en   ">

  <head>
	<meta http-equiv="Content-Type" content="text/html; charset=UTF-8" />
	<meta name="Keywords" content="LncDM R package">
	<title><?php echo $group_name; ?></title>
	<link href="http://<?php echo $themeroot; ?>styles/estilo1.css" rel="stylesheet" type="text/css" />
  </head>

<style>
   #title {
     margin:0 auto;
     height: 40px;
	 text-align:center;
	 box-shadow:-2px 0 6px #a9a9aa; /*◊Û±ﬂ“ı”∞*/
	 box-shadow:2px 0 6px #a9a9aa; /*”“±ﬂ“ı”∞*/
   }
   #title h2 {
     line-height:40px;
   }
   #contents {
     margin-top:10px;
     text-align:left;
   }
   #contents #content a{
     font-family:"TimesNewRoman",Georgia,Serif;
     font-size:18px;
     text-decoration: none;
	 color: #5e5f5f;
	 font-weight:bold;
   }
   #contents #content a:hover{
     text-decoration:underline;
   }
   h3 {
     text-align:left;
   }
   #No1 {
     margin-left:10px;
	 margin-right:10px;
   }
   #No2 {
     margin-left:10px;
	 margin-right:10px;
   }
   #No3 {
     margin-left:10px;
	 margin-right:10px;
   }
   #No3 a{
     color:#2822a4;
   }
   #tail {
     margin-left:10px;
	 margin-right:10px;
   }
   #tail a{
     color:#2822a4;
   }
</style>

<body>

<!-- R-Forge Logo -->
<table border="0" width="100%" cellspacing="0" cellpadding="0">
<tr><td>
<a href="http://r-forge.r-project.org/"><img src="http://<?php echo $themeroot; ?>/imagesrf/logo.png" border="0" alt="R-Forge Logo" /> </a> </td> </tr>
</table>


<!-- get project title  -->
<!-- own website starts here, the following may be changed as you like -->
<div id="title">
  <h2>LncDM:An R package for performing DNA methylation pattern of lncRNAs based on a reannotation method</h2>
</div>

<div id="No1">
  <h3>Introduction</h3>
  <p style="text-align:justify; text-justify:inter-ideograph;">
     LncDM is a novel computational method for identifying different methylation site,different methylation region and different methylation element in specific disease.Here we provide multiple options for user to analyze protein coding gene,lincRNA and other long non-coding RNA,pseudogene.LncDM is based on reannotation method to use Illumina HumanMethylation450 BeadChip and Gencode's all annotation information,and get different transcript's CpG value.We consider gene function region to caculate DM pattern.LncDM provide several test options and threshold value for users.Users can get data visualization in many aspects.  
  </p>
</div>

<div id="No2">
  <h3>Installation</h3>
  <p>1.The LncDM package requires R version >= 3.30 .</p>
  <p>2.The LncDM package requires the following packages to be installed: beanplot,reshape,gplots,WriteXLS,MASS,impute,limma,preprocessCore. If your system does not have them installed, the easiest way to install them is to issue the following command at the R prompt:
  source("http://bioconductor.org/biocLite.R"); biocLite(c("impute","limma","preprocessCore")); install.packages(c("beanplot","reshape","gplots","WriteXLS","MASS"),repos="http://cran.r-project.org");</p>
  <p>3.Install LncDM package:
  install.packages("LncDM", repos="http://R-Forge.R-project.org");</p>
</div>

<div id="No3">
  <h3>Tutorial</h3>
  <p>A vignentte that illustrate various aspects of LncDM is available <a href="download/LncDM.pdf"><strong>here</strong></a></p>
  <p>The user manual of LncDM package could be found <a href="download/LncDM-manual.pdf"><strong>here</strong></a></p>
</div>

<!-- end of project description -->


<div id="tail">
 <p> The <strong>project summary page</strong> you can find <a href="http://<?php echo $domain; ?>/projects/<?php echo $group_name; ?>/"><strong>here</strong></a>. </p>
 <p>Author: Hui Zhi&lt;zhihui013201@gmail.com&gt;</p>
</div>

</body>
</html>
