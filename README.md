# DiffCoEx
![](https://github.com/siddsabata/DiffCoEx/blob/main/data/diffcoex_app.png)

This was a group project done for 02-601, Programming for Scientists, taught at Carnegie Mellon University by Dr. Phillip Compeau during the Fall of 2024. 

Authors: Jason Hyun, Darrick Lo, Siddharth Sabata, and Katherine Wang

### Background
In this project, we created a streamlined Go-R data pipeline to speed up data processing and significance testing of differential gene coexpression results from microarray data. A user friendly interface was created using R-Shiny to simplify viewing the results. 

We implemented two different algorithms: [DiffCoEx](https://bmcbioinformatics.biomedcentral.com/articles/10.1186/1471-2105-11-497), and [CoXpress](https://bmcbioinformatics.biomedcentral.com/articles/10.1186/1471-2105-7-509), and compared clustering results between the two. Significance testing was then performed on DiffCoEx results, due to the significantly smaller number of gene modules. 

Note: Generative AI used for this project

### Package dependencies
Make sure you have RStudio installed with the latest version of R. 

WGCNA

Shiny

CoXpress

RColorBrewer

preprocessCore

flashClust

ggplot2

tidyverse

Also make sure to have Bioconductor installed. 

Installing CoXpress: 
From https://sourceforge.net/projects/coxpress/files/coxpress/R%203.0/ , download coXpress_1.5.tar.gz (the tar file, not the zip file). 
RStudio > Tools > Install Packages... > select "Install from: Package Archive File (.zip; .tar.gz) > select the tar file.

To run this application, open `app.R` and press run. 
### Testing 
When running `go test`, make sure to go into the subdirectories to test all of the supporting executables. The root directory contains the testing for the data processing done in Go, and everything else is organized in their respective directories. 
### Usage 
Your computer will not give the executables from Go permission to run originally. You have to manually give each executable permission to run after pressing each button (except for clustering). You should get a warning that pop up regarding permissions to run an executable. On a Mac you can manually give permission to each executable after it gets blocked by pressing the allow button. This is under the Privacy tab. Once you have given permission to the executable, press the button again and it will work. 

Use the demonstration video as reference in using the app. 

NOTE: the rat data will take a very long time to run. Golub is much faster. 

NOTE: there is no indication that the significance testing process has begun. Simply press the button once and wait. Do not over load the app with instructions. 
