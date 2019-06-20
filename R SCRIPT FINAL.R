####### R SCRIPT FOR ANALYSIS OF DATA IN COOPER ET AL (2019) #######

##  Reference:  Cooper, K.M., Bolam, S.G., Downie, A.-L., Barry J. (2019) Biological-based habitat  classification approaches promote cost-efficient monitoring: an example using seabed assemblages. Journal of Applied Ecology

## Required Files: 
# 1. R SCRIPT FINAL.R (i.e. this script. Available from: https://doi.org/10.14466/CefasDataHub.56). 
# 2. C5922DATASET13022017REDACTED.csv (Available from https://doi.org/10.14466/CefasDataHub.34)
# 2. C5922DATASETFAM13022017REDACTED.csv (used to generate file: C5922DATASETFAM13022017REDACTED.csv Also available from https://doi.org/10.14466/CefasDataHub.56). 
# 3. PARTBAGG05022017.csv (Available from https://doi.org/10.14466/CefasDataHub.34)
# 4. EUROPE.shp (European Coastline) (Available from https://doi.org/10.14466/CefasDataHub.34)
# 5. EuropeLiteScoWal.shp (European Coastline with UK boundaries) (Available from https://doi.org/10.14466/CefasDataHub.34)
# 6. UKSeaMap2016_SedimentsDissClip.shp ((Available from:https://doi.org/10.14466/CefasDataHub.56)
# 7. StudyArea.shp (Available from:https://doi.org/10.14466/CefasDataHub.56)
# 8. AvCur_Final.tif (Not provided - see originator of data, details in table 1 of the paper)
# 9. Chla_Final.tif (Not provided - see originator of data, details in table 1 of the paper)
# 10. Depth_Final.tif (Not provided - see originator of data, details in table 1 of the paper)
# 11. Gravel_Final.tif (Not provided - see originator of data, details in table 1 of the paper)
# 12. Mud_Final.tif (Not provided - see originator of data, details in table 1 of the paper)
# 13. Sal_Final.tif (Not provided - see originator of data, details in table 1 of the paper)
# 14. Sand_Final.tif (Not provided - see originator of data, details in table 1 of the paper)
# 15. SPM_Final.tif (Not provided - see originator of data, details in table 1 of the paper)
# 16. Stress_Final.tif (Not provided - see originator of data, details in table 1 of the paper)
# 17. Temp_Final.tif (Not provided - see originator of data, details in table 1 of the paper)
# 18. WOV_Final.tif (Not provided - see originator of data, details in table 1 of the paper)
# 19. FaunalCluster.tif (Available from:https://doi.org/10.14466/CefasDataHub.56)
# 20. PhysicalCluster.tif (Available from:https://doi.org/10.14466/CefasDataHub.56)
# 21. FaunalClusterClip.tif (Available from:https://doi.org/10.14466/CefasDataHub.56)
# 22. PhysicalClusterClip.tif (Available from:https://doi.org/10.14466/CefasDataHub.56)

#NB (.tifs created in: CDP:\C7331_Marine Aggregate Reg\Working_Area\SPECIES DIST MOD\R\SDM13.R)

## Required folder structure:
# C:\FINAL\R 
#               \DATA
#                   \FINAL_MODEL
#                       \FaunalCluster
#                       \PhysicalCluster
#                   \FINAL_RASTERS
#                   \UKSeaMap2016
#               \OUTPUTS

## Folder Contents: 
# C:\FINAL\R (file 1) 
# C:\FINAL\DATA (files 2:5)
# C:\FINAL\DATA\UKSeaMap2016\ (file: 6)
# C:\FINAL\DATA\FINAL_RASTERS (files 8:18)
# C:\FINAL\DATA\FINAL_MODEL (file: 7)
# C:\FINAL\DATA\FINAL_MODEL\FaunalCluster (file 19)
# C:\FINAL\DATA\FINAL_MODEL\FaunalCluster (file 20)
# C:\FINAL\OUTPUTS (.png and .csv files resulting from script)

## Notes:  # The script will generate the figures and content of tables in Cooper et al (2018). The script is divided into 57 individual steps. The dataset employed in this study is initially built,using steps 2-7, from the raw data file C5922DATASET13022017REDACTED.csv .As the resulting family level dataset is made available it is possible to omit these steps.Note also that where raster data (required files 8-18) are unavailable then omit steps 17-26 and 28-35. The outputs from these steps (FaunalCluster.tif,PhysicalCluster.tif,FaunalClusterClip.tif,PhysicalClusterClip.tif) are provided for use with this script.

## Set working directory
setwd('C:/Users/kmc00/OneDrive - CEFAS/C6264 FINAL')


#### 1. PREPARE MAPPING LAYERS ####

## Load packages
library(ggplot2)
library(rgdal)
library(maptools)
library(plyr)
library(rgeos)
library(sp)
library(rgdal)
library(raster)

## Produce a high definition European coast map
# Load shapefile
eu <- readOGR("DATA","EUROPE")
eu@data$id <- rownames(eu@data)

## Create a data.frame from our spatial object
euDF <- fortify(eu, region = "id")

## Merge the "fortified" data with the data from our spatial object. Object euDFis a dataframe of the polygon coordinates.  This is combined with the other attribute data from object eu
euDF2 <- merge(euDF, eu@data, by = "id")

## Create layer for coastline including borders between UK countries
eubound = readOGR("DATA","EuropeLiteScoWal")
eubound@data$id = rownames(eubound@data)
eubound.points = fortify(eubound, region="id")
eubound.df = join(eubound.points, eubound@data, by="id")

## Create layer for EUNIS. Note that 'readOGR' doesn't work with this shapefile so have to use 'readShapePoly'instead
EUNIS=readShapePoly("DATA/UKSeaMap2016/UKSeaMap2016_SedimentsDissClip.shp")
EUNIS@data$id = rownames(EUNIS@data)
EUNIS.points = fortify(EUNIS, region="id")
EUNIS.df = join(EUNIS.points, EUNIS@data, by="id")

## Bring in study area polygon
StudyArea = readOGR("DATA/FINAL_MODEL/StudyArea.shp")
StudyArea@data$id = rownames(StudyArea@data)
StudyArea.points = fortify(StudyArea, region="id")
StudyArea.df = join(StudyArea.points, StudyArea@data, by="id")
#View(StudyArea.df)


#### 2. BUILD FAMILY LEVEL DATASET: BRING IN RAW DATA ####

## Load raw data (https://doi.org/10.14466/CefasDataHub.34). If using C5922DATASETFAM13022017REDACTED.csv skip to step'LOAD DATA' below.
dataall2=read.csv("DATA/C5922DATASET13022017REDACTED.csv", header=T,na.strings=c("NA", "-","?","<null>"),
                  stringsAsFactors=F,check.names=FALSE)

## View dataset
#View(dataall2)

## Remove 1st column (repeat of row labels)
dataall2 = dataall2[,-1]

## View dataset
#View(dataall2)# check 1st column removed

## Check dataset dimensions
dim(dataall2)# 33198 13587

## Check column order
options(max.print=1000000)# to see all column names
names(dataall2) # Structure: Faunal data (A & B) followed by other variables  

## Select only variables of interest
# 1:13450 Faunal data
# 13451 Order
# 13452 SampleCode
# 13453 SurveyName
# 13455 Latitude_WGS84
# 13456 Longitude_WGS84
# 13458 Gear
# 13460 Sieve
# 13461 Year
# 13462 Month
# 13465:13563 Sed data
# 13564 Sector
# 13565 Source
# 13569 PSASubSample
# 13572 Treatment
# 13575 Sal
# 13576 Temp
# 13577 Chla
# 13578 SPM
# 13579 Depth
# 13580 WOV
# 13581 AvCur
# 13582 Stress
# 13583 IncCol
# 13584 TopSieveEmpty
# 13585 Data
# 13586 Programme
# 13587 WentworthSuitability
dataall3=dataall2[,c(1:13450,13465:13563,13451:13453,13455:13456,13458,13460:13462,13564:13565,
                     13569,13572,13575:13587)]

## Check column names
names(dataall3)

## Check dataset dimensions
dim(dataall3) # 33198 (samples) 13575 (variables)


#### 3. BUILD FAMILY LEVEL DATASET: SPLIT RAW DATA INTO SECTIONS A AND B ####

## The dataset is too large to manipulate as a whole, so it must be split into2 chunks(A & B)

## Create a dataframe object 'dataparta' for Part A
dataparta=dataall3[1:21101,]# samples 1:21101

## Create a dataframe object 'datapartb' for Part B
datapartb=dataall3[21102:33198,]# samples 21102:33198

## Save parts A and B as .csv files
write.csv(dataparta,file = "DATA/dataparta.csv",row.names=TRUE)
write.csv(datapartb,file = "DATA/datapartb.csv",row.names=TRUE)

## Remove objects that are no longer required to free up memory
rm(dataall2,dataall3)
gc()
rm(dataparta,datapartb)
gc()


#### 4. BUILD FAMILY LEVEL DATASET: GENERATE TAXON MATRIX (FAMILY LEVEL) FOR PART A ####

##  If you experience memory problems close R and other programmes and start from this Step (having set the working directory)           

## Load the data for Part A
dataparta=read.csv("DATA/dataparta.csv", header=T,na.strings=c("NA", "-", "?","<null>"),
                   stringsAsFactors=F)

## View df 'dataparta'
#View(dataparta)

## Remove 1st column
dataparta$X=NULL

## Check dimensions of Part A - should be 21101 13575
dim(dataparta)# it is

## View names for df 'dataparta'
options(max.print=1000000)# to see all column names
names(dataparta)

## Select only the faunal data. Note that this includes columns for additional taxa from part B
# (all with zero abundance)
dataparta2=dataparta[,1:13450]

## View df 'dataparta2
View(dataparta2)

## Change NAs to 0. NAs occur for taxa relating to part B
dataparta2[is.na(dataparta2)] <- 0

## Transpose the faunal data (keep column names)
t_dataparta = setNames(data.frame(t(dataparta2[,-1])), dataparta2[,1])

## Bring in aggregation data. Note that Part A includes full taxon list, hence use of the Part B
# Aggregation file
agg=read.csv("DATA/PARTBAGG05022017.csv", header=T,
             na.strings=c("NA", "-","?","<null>"),stringsAsFactors=F,check.names=FALSE)

## View aggreation file
#View(agg)

## Join the aggregation file to data
alla=cbind(agg,t_dataparta)

## Remove unwanted taxa
alla2 = subset(alla, Include=='Y')

## Identify dimensions of df 'alla2'
dim(alla2)# note the total number of columns

## View df 'alla2'
#View(alla2)# identify the col nos for Family and samples

## Take only the Family and sample columns
alla2fam=alla2[,c(19,31:21131)]

## Remove objects that are no longer required to free up memory
rm(alla)
gc()
rm(alla2)
gc()
rm(t_dataparta)
gc()
rm(dataparta)
gc()
rm(dataparta2)
gc()

## Aggregate data (sum) by family for multiple samples
famabundPARTA=aggregate(. ~ Family, alla2fam, sum)

## Remove objects that are no longer required to free up memory
rm(alla2fam)
gc()
rm(agg)
gc()

## Write file and close
write.csv(famabundPARTA,file = "DATA/famabundPARTA.csv",row.names=TRUE)


#### 5. BUILD FAMILY LEVEL DATASET: GENERATE TAXON MATRIX (FAMILY LEVEL) FOR PART B ####

## Load the data for Part B
datapartb=read.csv("DATA/datapartb.csv", header=T,na.strings=c("NA", "-", "?","<null>"),
                   stringsAsFactors=F)

## Remove 1st column
datapartb$X=NULL

## View names for df 'datapartb'
options(max.print=1000000)# to see all column names
names(datapartb)

## Select only the faunal data
datapartb2=datapartb[,1:13450]

## Transpose the data (keep column names)
t_datapartb = setNames(data.frame(t(datapartb2[,-1])), datapartb2[,1])

## Bring in aggregation data
agg=read.csv("DATA/PARTBAGG05022017.csv", header=T,na.strings=c("NA", "-", "?","<null>"),
             stringsAsFactors=F)

## Join aggregation file to data
allb=cbind(agg,t_datapartb)

## Remove unwanted taxa
allb2 = subset(allb, Include=='Y')
names(allb2)

## Take only the Family col and the sample cols (NB Impossible to aggregate other factor cols)
allb2fam=allb2[,c(19,31:12127)]
dim(allb2fam)#12913 12098

## Change NAs to 0. 
allb2fam[is.na(allb2fam)] <- 0

## Remove objects that are no longer required to free up memory
rm(agg,allb,allb2,datapartb,datapartb2,t_datapartb)
gc()

## Aggregate data (sum) by family for multiple samples
famabundPARTB=aggregate(. ~ Family, allb2fam, sum)
names(famabundPARTB)

## Remove objects that are no longer required to free up memory
rm(allb2fam)
gc()

## Write file
write.csv(famabundPARTB,file = "DATA/famabundPARTB.csv",row.names=TRUE)

## Remove objects that are no longer required to free up memory
rm(famabundPARTA,famabundPARTB)
gc()


#### 6. BUILD FAMILY LEVEL DATASET: BUILD FAMILY ABUND MATRIX FOR PARTS A AND B ####

## Load data from Part A
famabundPARTA=read.csv("DATA/famabundPARTA.csv", header=T,na.strings=c("NA", "-","?","<null>"),
                       stringsAsFactors=F,check.names=FALSE)

## View df 'famabundPARTA'
#View(famabundPARTA)

## Remove 1st column
famabundPARTA[1]=NULL

## Check dimensions of df 'famabundPARTA'
dim(famabundPARTA)#774 21102

## Load data from Part B
famabundPARTB=read.csv("DATA/famabundPARTB.csv", header=T,na.strings=c("NA", "-","?","<null>"),
                       stringsAsFactors=F,check.names=FALSE)

## Remove 1st column
famabundPARTB[1]=NULL

## View df 'famabundPARTB'
#View(famabundPARTB)

## Check dimensions of df 'famabundPARTB'
dim(famabundPARTB)#774 12098

## Merge the dataframes 'famabundPARTA' and 'famabundPARTB'
total <- merge(famabundPARTA,famabundPARTB,by="Family")

##Check dimensions of df 'total'
dim(total)#774 33199

## View df 'total'. Note samples as columns
#View(total)

## Transpose df 'total' so samples as rows and variables as columns
t_total = setNames(data.frame(t(total[,-1])), total[,1])

## Check orientation of df 't_total'
#View(t_total)

## Check names of df 't_total'
names(t_total)

## Change row names (sample codes) to a column
t_total <- cbind(Sample = rownames(t_total), t_total)
rownames(t_total) <- NULL

## Check new column (Sample) added and old row names (Sample codes)removed
#View(t_total)# all correct

## Write file
write.csv(t_total,file = "DATA/t_total.csv",row.names=TRUE)


#### 7. BUILD FAMILY LEVEL DATASET: ADD OTHER VARIABLES TO FAM ABUND MATRIX ####

## Load required dfs
dataparta=read.csv("DATA/dataparta.csv", header=T,na.strings=c("NA", "-", "?","<null>"),
                   stringsAsFactors=F)
datapartb=read.csv("DATA/datapartb.csv", header=T,na.strings=c("NA", "-", "?","<null>"),
                   stringsAsFactors=F)

## Remove 1st columns of dfs 'dataparta' and 'datapartb'
dataparta$X=NULL
datapartb$X=NULL

## View names of df 'dataparta'
names(dataparta)
names(datapartb)

## Select only the non-faunal data (i.e. meta and sediment data)
dataparta2other=dataparta[,c(1,13451:13575)]
datapartb2other=datapartb[,c(1,13451:13575)]

## Join together df 'dataparta2other' and 'datapartb2other' (other variables -  parts A and B)
datapartab2other=rbind(dataparta2other,datapartb2other)

## Merge faunal (family) matrix with other variables
data <- merge(t_total,datapartab2other,by="Sample")

## Check dimensions of df 'data'
dim(data)# 33198 900

## Check names of df 'data'
names(data)

## View df 'data'
#View(data)

## Save df 'data' (Sample/variable matrix, Family level)
write.csv(data,file = "DATA/C5922DATASETFAM13022017REDACTED.csv",row.names=TRUE)

## Remove objects that are no longer required to free up memory
rm(data,dataparta,dataparta2other,datapartab2other,datapartb,datapartb2other,famabundPARTA,
   famabundPARTB,t_total,total)
gc()


#### 8. LOAD DATA ####

## Read in data from csv file. This family level dataset was created using steps in 'BUILD FAMILY LEVEL DATASET' - see above
data=read.csv("DATA/C5922DATASETFAM13022017REDACTED.csv", header=T,na.strings=c("NA", "-","?","<null>"),stringsAsFactors=F,check.names=FALSE)

## Remove 1st col
data[1] <- NULL 

## Check dataset dimensions
dim(data)# 33198 900
names(data)

## Identify cols which should be type 'numeric'
names(data) #Faunal data (2:775) sed data (776:874) coordinates (878:879) explanatory variable(888:895)
## Vector for cols which need to be numeric
ix <- c(2:775,776:874,878:879,888:895)
data[ix] <- lapply(data[ix], as.numeric) 
str(data[888:895])

## Identify cols which should be type 'character'
names(data)# sample (1),sorder (875)  samplecode(876), surveyname (877), gear (880), sieve (881), sector (884), source (885), psasubsample (886), incol (896), topsieveempty (897),data (898),programme (899),wentworthsuitability (900)
ix2 <- c(1,875,876,877,880,881,884,885,886,896,897,898,899,900)
data[ix2] <- lapply(data[ix2], as.character) 

## Identify cols which should be type 'factor'
names(data)# year (882), month (883), treatment (887) 
ix3 <- c(882,883,887)
data[ix3] <- lapply(data[ix3], as.factor)


#### 9. PREPARE DATA  ####

## First remove known 'impact' samples
dim(data)# 33198   900
data1=subset(data, Treatment=="R"|Treatment=="NA")
dim(data1)# 29202   900

## Change NAs in sieve data cols (776:874) to zero so you can sum
names(data1)
data1[, 776:874][is.na(data1[, 776:874])] <- 0

## Add in sed summary variables for Gravel, Sand and Mud
data1$Gravel=rowSums(data1[,776:808])# 125mm - 2mm
data1$Sand=rowSums(data1[,809:835])# 1.7mm - 0.0625mm
data1$Mud=rowSums(data1[,836:874])# 0.045mm - Pan

## Add column for TotalPercent
data1["TotalPercent"]=NA
names(data1)# Check col added - it is

## Populate col TotalPercent by summing across the sieve cols: 125mm (776) to Pan (874)
data1$TotalPercent=rowSums(data1[,776:874])

## Subset data for TotalPercent is ~100
data1.1 <- subset(data1, TotalPercent > 99 & TotalPercent < 101)
dim(data1.1) # 23601   904

## Remove rows where there is a value of NA in one of the phy variable columns: Sal, Temp, Chla, SPM, Depth, WOV, AvCur,Stress, Gravel, Sand, Mud.
names(data1.1)
data1.2=data1.1[complete.cases(data1.1[,c(888:895,901:903)]),]
dim(data1.2)# 22011   904

## Create a df 'test.points' which are the coordinates (lon, lat) from df data1.2
names(data1.2)
test.points=data1.2[,879:878]

## Define coordinates for df 'test.points'
coordinates(test.points) <- ~Longitude_WGS84+Latitude_WGS84

## Define CRS for df 'test.points'
proj4string(test.points) <- CRS(proj4string(StudyArea))

## Id samples within study area. NAs signify points outside study area
ss=over(test.points, StudyArea)
names(ss)

## Combine df data1 with results of in/out assessment
data1ss=cbind(data1.2, ss)

## Change name of col 'id' to 'WithinStudyArea'
names(data1ss)
colnames(data1ss)[905]="WithinStudyArea"

#### Remove points which are out of shapefile (or raster) extent
#https://gis.stackexchange.com/questions/181586/how-to-remove-points-which-are-out-of-shapefile-or-raster-extent 
dim(data1ss)# 22011   905
data1.5=data1ss[!is.na(data1ss$WithinStudyArea),]
dim(data1.5)# 21939   905

## Drop col 'WithinStudyArea'
data1.5$WithinStudyArea<- NULL
names(data1.5)

## Subset data by gear (all 0.1m2 grabs)
data2 = subset(data1.5, Gear=="MHN" | Gear=="DG" | Gear=="VV" | Gear=="SM")
dim(data2)# 21314   904

## Subset by sieve size (1mm only)
data3 = subset(data2, Sieve=='1')
dim(data3)# 19409   904

## Remove samples from the faunal data (df data3) where no fauna present
names(data3) # 2:775
data4 = data3[ rowSums(data3[,2:775])!=0, ]
dim(data4)# 19314   904


#### 10. REMOVE 'REPLICATES' TO ADDRESS SPATIAL AUTOCORRELATION ####

## Set coordinates
coordinates(data4) <- c("Longitude_WGS84", "Latitude_WGS84")


## Work out 50m distanc in decimal degrees. 1 degree of latitude =111,000m
50/111000# degrees for 50m #0.0004504505

## Set distance within which to remove replicates (zero)
zd <- zerodist(data4,zero = 0.0004504505)

## Drop replicates
data4norep <- data4[-zd[,2], ]
dim(data4norep)# 13742   902

## Change class to df
data4.2=data.frame(data4norep)
class(data4.2)
names(data4.2)

## Drop col 'optional'
data4.2=data4.2[,1:904]

#data4.2=data4 # just to get SA 1st round
#names(data4)


#### 11. PREPARE DATA CONT'D ####

## Create matrix for faunal data
dat=as.matrix(data4.2[,2:775])

## Identify columns with no abundance records (i.e. all zeros)
i<- (colSums(dat, na.rm=T) != 0)

## Remove columns with all zeros (i.e. taxon with no abundances)
matnonzero <- dat[, i]

## Change df 'matnonzero' back to df
matnonzero=as.data.frame(matnonzero)
 
## Stitch data back together
data4.5=cbind(data4.2$Sample,matnonzero,data4.2[,776:904])

## Change name of 1st column to 'sample
colnames(data4.5)[1]="sample"

## Check dimensions of df 'data4.5'
dim(data4.5)# 13742   765

## Show names of df 'data4.5'
names(data4.5)  


#### 12. SPLIT DATA INTO TRAIN AND TEST SUBSETS ####

## 70% of the sample size
smp_size <- floor(0.7 * nrow(data4.5))

## Randomly sample dataset
set.seed(123)#set the seed to make your partition reproducible
train_ind <- sample(seq_len(nrow(data4.5)), size = smp_size)

## Create a train and a test set
traindat <- data4.5[train_ind, ]
testdat <- data4.5[-train_ind, ]

## Check dimensions of train and test sets
dim(traindat)#  9619  765
dim(testdat)#  4123  765
names(traindat)

## Plot trains and test sites to make sure similar spatial coverage
names(traindat) # note cols for lat long
coordtraindat=traindat[,740:739]
coordtestdat=testdat[,740:739]

## Plot study area and train and test sets
plot(StudyArea)
points(coordtraindat, col = "blue", cex = 0.01)
points(coordtestdat, col = "red", cex = 0.01)

## Create a faunal data subset for train set 
names(traindat)
data5=traindat[,2:636]
dim(data5)# 9619  635

## Remove any faunal variable with no abundance records
names(data5)
dattrain=as.matrix(data5[,1:635])# 1st create a matrix

## Identify cols with no data - results to object i2
i2<- (colSums(dattrain, na.rm=T) != 0)

## Remove empty cols
matnonzerotrain <- dattrain[, i2]# drop zero cols

## Change df 'matnonzerotrain' from matrix to df
matnonzerotrain=as.data.frame(matnonzerotrain)
dim(matnonzerotrain)# 9619  613

## Change name of object 'matnonzerotrain' to 'data5'
data5.1=matnonzerotrain
dim(data5.1)
names(data5.1)

## Change class of df data5 to a matrix
data6=data.matrix(data5.1)

## Create a df 'pos' for Sample, Latitude_WGS84 and Longitude_WGS84 
names(traindat)
pos=traindat[,c(1,739:740)]
names(pos)


#### 13. FAUNAL CLUSTERING ####

## Transform the data (fourth-root transformation)
datat=data6^(0.25)

## Perform k-means clusterinig of data. Results (cluster group) to the object 'results'
set.seed(1234)
results=kmeans(datat,12,algorithm="MacQueen",iter.max=100,nstart=25)

## Number of samples belonging to each cluster group
results$size#  387  317  832  410  643 1413  337 1266 2719  748  134  413

## Add cluster group from kmeans results file to df 'pos' which includes 'Sample',
# 'Latitude_WGS84' and 'Longitude_WGS84'
faunal.cluster=cbind(pos,results$cluster)

## Change name of col 'results$cluster' to 'ClusterNum'
names(faunal.cluster)[4]<-paste("ClusterNum")

## Add a new empty col 'FaunalCluster' to df 'faunal.cluster
faunal.cluster["FaunalCluster"]=NA
dim(faunal.cluster)# 9619    5
names(faunal.cluster)

## Populate FaunalCluster col with new names (see dendrogram Fig 3a in Cooper & Barry 2017)
faunal.cluster$FaunalCluster[faunal.cluster$ClusterNum == 11] <- "A1"
faunal.cluster$FaunalCluster[faunal.cluster$ClusterNum == 2]<- "A2a"
faunal.cluster$FaunalCluster[faunal.cluster$ClusterNum == 3] <- "D2b"
faunal.cluster$FaunalCluster[faunal.cluster$ClusterNum == 9]<- "D2c"
faunal.cluster$FaunalCluster[faunal.cluster$ClusterNum == 7] <- "B1b"
faunal.cluster$FaunalCluster[faunal.cluster$ClusterNum == 8] <- "D2d"
faunal.cluster$FaunalCluster[faunal.cluster$ClusterNum == 4] <- "A2b"
faunal.cluster$FaunalCluster[faunal.cluster$ClusterNum == 1] <- "B1a"
faunal.cluster$FaunalCluster[faunal.cluster$ClusterNum == 12] <- "D1"
faunal.cluster$FaunalCluster[faunal.cluster$ClusterNum == 5] <- "C1b"
faunal.cluster$FaunalCluster[faunal.cluster$ClusterNum == 10] <- "C1a"
faunal.cluster$FaunalCluster[faunal.cluster$ClusterNum == 6]<- "D2a"

## RS Map of faunal cluster distribution with all train samples
p2= ggplot()+
  geom_point(data=faunal.cluster,aes(Longitude_WGS84,Latitude_WGS84,col=FaunalCluster),size=0.45,show.legend = TRUE)+#
  geom_polygon(data = euDF2, aes(x=long, y=lat, group = group),fill="white",colour="black",size=0.1)+
  scale_colour_manual(values = c("blue2","cyan1","#05aae1","plum2","darkorchid3","green3","palegreen1","#b40202","red1","darkorange","yellow","#b4b404"),name="Cluster")+
  guides(colour = guide_legend(override.aes = list(size=3)))+ # Change size of legend dots
  coord_map(xlim = c(-10.7, 4),ylim = c(48, 62))+ #set x,y limits of plot
  theme_bw(base_size = 24)+ 
  labs(x="Longitude",y="Latitude")

fig4a=p2+theme(legend.key.size = unit(1, "cm"))+
  guides(colour = guide_legend(override.aes = list(size=6)))

fig4a


#### 14. PHYSICAL CLUSTERING ####

## Create a df 'phy.data' which is the subset of chosen variables (see above) from df 'traindat'. Variables of interest are: Sample (1), Latitude_wgs84 (739), Longitude_wgs84 (740),Sal (749), Temp (750), Chla (751), SPM (752), Depth (753), WOV (754), AvCur (755), Stress (756), Gravel (762), Sand (763), Mud (764)
names(traindat)
phy.data=traindat[,c(1,739,740,749:756,762:764)]
names(phy.data)# Check selected columns present
str(phy.data)

## Change bathy values to positive
phy.data$Depth=abs(phy.data$Depth)

## Remove variables from df phy.data that won't be used for clustering (i.e. Sample, Latitude_WGS84, and Longitude_WGS84)
phy.data.final=phy.data[,4:14] 
names(phy.data.final)# Check unwanted varaibales removed
dim(phy.data.final)#9619   11

## Pairs plot
#View(phy.data.final)
pairs(phy.data.final)
res <- cor(phy.data.final)
round(res, 2)# Note correlations> 0.75 for: Chla/SPM and Gravel/Sand. Therefore remove Chla and Sand
phy.data.final2=phy.data.final[,c(1,2,4,5,6,7,8,9,11)]
str(phy.data.final2)
#View(phy.data.final2)

## Identify any skewness in the data. Intuitively, the skewness is a measure of symmetry. As a rule, negative skewness indicates that the mean of the data values is less than the median, and the data distribution is left-skewed. Positive skewness would indicates that the mean of the data values is larger than the median, and the data distribution is right-skewed.
summary(phy.data.final2)# RS (Mean>median) observed for SPM (3), Depth (4), Stress (7), Gravel (8), Mud (9)

## Create a copy of df phy.data.final
phy.data.finalT=phy.data.final2

## Transform relevant columns
phy.data.finalT[,c(3,4,7,8,9)]=log(phy.data.finalT[,c(3,4,7,8,9)]+0.1)
#View(phy.data.finalT)

## Normalisation of phy variable data (subtract by pop mean and sd)
df <- scale(phy.data.finalT)# scale the data 

## Dimensions of df 'df'
dim(df)# 9619   9
#View(df)
class(df)

## Remove Sand from the model as it covaries with gravel
#df2=as.data.frame(df)
#df3=df2[,c(1:9,11)]
#names(df3)

## kmeans clutering of physical variables
set.seed(1234)
phy.results=kmeans(df,12,algorithm="MacQueen",iter.max=100,nstart=25)

## Create a df 'phy.cluster' for the cluster groups
phy.cluster=as.data.frame(phy.results$cluster)

## Create an object 'phy.output' with sample coordinates and cluster group.
phy.output=cbind(phy.data$sample,phy.data$Latitude_WGS84,
                 phy.data$Longitude_WGS84,phy.cluster)

## Reinstate the appropriate column names
colnames(phy.output) <- c("Sample", "Latitude_WGS84","Longitude_WGS84","PhyCluster")

## Change variable PhyCluster from an integer to a factor
phy.output$PhyCluster=as.factor(phy.output$PhyCluster)

## Check names of df 'phy.output'
names(phy.output)

## Number of samples belonging to each cluster group
phy.results$size#  1405  923  633  589  910  284  336  618 1181 1320  359 1061

## Produce a map of physical cluster group distributuion
PhyClusMap= ggplot()+
  geom_polygon(data = euDF2, aes(x=long, y=lat, group = group),fill="white",colour="black",size=0.05)+
  geom_point(data=phy.output,aes(Longitude_WGS84,Latitude_WGS84,col=PhyCluster), size=0.45,show.legend = TRUE)+
  coord_map(xlim = c(-10.7, 4),ylim = c(48, 62)) + #set x,y limits of plot
  theme_bw(base_size=24)+
  guides(colour = guide_legend(override.aes = list(size=3)))+ # Change size of legend dots
  labs(x="Longitude",y="Latitude")#

fig8a=PhyClusMap+theme(legend.key.size = unit(1, "cm"))+
  guides(colour = guide_legend(override.aes = list(size=8))) 

fig8a


#### 15. FIGURE S1: PHY CLUSTER VARIABLE BOX PLOTS  ####

## Create a df with required data: Sal,Temp,Chl a, SPM, Depth, WOV, AvCur, Stress
# (from df 'phy.data.final') and Sample, latitude_wgs84, longitude_wgs84 and PhyCluster (phy cluster group) (from df 'phy.output).
bp=cbind(phy.data.final,phy.output)
names(bp)
dim(bp)# 9619   15

## Change depths back to negative values
#View(bp)
bp$Depth=bp$Depth*-1

## Create a df 'bp2' with only the phy variables and PhyCluster group
bp2=bp[,c(1:11,15)]

## Check required variable present
names(bp2)

## Change colnames for bp
colnames(bp2)= c("Sal","Temp","Chl a","SPM","Depth","WOV","AvCur"     
                 ,"Stress","Gravel","Sand","Mud","PhyCluster")

## Call library 'reshape2' - required for melting data
library(reshape2)

## Melt data into suitable form for facetting 
dat <- melt(bp2,"PhyCluster") # Produces a df of 3 cols (PhyCluster, variable, value)
#View(dat)

## Plotting
phyboxplots=ggplot(dat, aes(PhyCluster,value,fill=PhyCluster)) + 
  geom_boxplot(outlier.shape=NA) + # don't display outliers
  theme_bw(base_size = 24)+
  facet_wrap(~variable,scales="free_y")+
  theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.4)) +
  labs(y = "Value")

figS1=phyboxplots+theme(legend.position="none")

## Save plot to an image file (png or tiff)
png("OUTPUTS/Figure_S1.png", width = 1200, height = 933,units = "px", pointsize = 24)
#tiff("OUTPUTS/Figure_S1.tif", width = 1200, height = 933,units = "px", pointsize = 24,compression = "lzw")

figS1
dev.off()


#### 16. FIGURE S2: PHY VARIABLE BOX PLOTS FOR FAUNAL CLUSTERS ####

## Find required variables
names(traindat) # Sample (1),Latitude_WGS84 (739),Longitude_WGS84 (740), Sal (749),Temp (750), Chla (751), SPM (752), Depth (753), WOV (754), AvCur (755), Stress (756), Gravel (762), Sand (763), Mud (764)
names(faunal.cluster)# ClusterNum (4), FaunalCluster (5)
names(phy.output)# PhyCluster (4)

## Check objects have same number of rows
dim(traindat) # 9619  765
dim(phy.output)# 9619    4
dim(faunal.cluster)# 9619    5

## Bring data together in one df
allres=cbind(traindat[,c(1,739,740,749:756,762:764)],faunal.cluster[,c(4,5)],phy.output[4])
#View(allres)
names(allres)

## Select cols for FaunalCluster(16), Sal(4), Temp(5),Chla(6), SPM(7), Depth(8), WOV(9), AvCur(10), Stress (11), Gravel (12), Sand (13), Mud (14)
data4fbp=allres[,c(16,4:11,12:14)]
names(data4fbp)
str(data4fbp)

## Change FaunalCluster from chr to factor
data4fbp$FaunalCluster=as.factor(data4fbp$FaunalCluster)
str(data4fbp)
dim(data4fbp)# 9619   12

## Update column names for df 'data4fbp'
colnames(data4fbp)=c("FaunalCluster","Sal","Temp","Chl a","SPM","Depth","WOV","AvCur","Stress","Gravel","Sand", "Mud")

## Call library 'reshape2' - required for melting data
library(reshape2)
library(ggplot2)

## Melt data into suitable form for facetting 
dat4 <- melt(data4fbp,"FaunalCluster") # Produces a df of 3 cols (PhyCluster, variable, value)
#View(dat4)

## Plotting
plot12=ggplot(dat4, aes(FaunalCluster,value,fill=FaunalCluster)) + 
  scale_fill_manual(values=c("blue2","cyan1","#05aae1","plum2","darkorchid3","green3","palegreen1","#b40202","red1","darkorange","yellow","#b4b404"))+
  geom_boxplot(outlier.shape=NA) + # don't display outliers
  theme_bw(base_size = 24)+
  theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.4)) +#
  facet_wrap(~variable,scales="free_y")+# was"free"
  labs(y = "Value")

plot13=plot12+theme(legend.position="none")

## Save plot to an image file (png or tiff)
png("OUTPUTS/Figure_S2.png", width = 1200, height = 933,units = "px", pointsize = 24)
#tiff("OUTPUTS/Figure_S2.tif", width = 1200, height = 933,units = "px", pointsize = 24,compression = "lzw")
plot13
dev.off()


#### 17. RANDOM FOREST (BIO): PREPARE DATA ####

## Start with clusteringresults df
names(allres)

## Select cols for: Sample (1), Longitude_WGS84 (3),Latitude_WGS84 (2), FaunalCluster(16), PhyCluster (17)
SMDdata=allres[,c(1,3,2,16,17)]
#View(SMDdata)#
names(SMDdata)

## In case you need to save/reload data
#write.csv(SMDdata, file = "DATA/SMDdata.csv",row.names=FALSE)
#SMDdata=read.csv("DATA/SMDdata.csv", header=T,na.strings=c("NA", "-","?","<null>"),stringsAsFactors=F,check.names=FALSE)

## Get columns in correct order
FaunalCluster=SMDdata[,c(1:4)]
#View(FaunalCluster)

## Change names of cols
colnames(FaunalCluster)=c("Sample","lon","lat","cluster")
head(FaunalCluster)

## Swap cluster labels for numbers
FaunalCluster$cluster[FaunalCluster$cluster == "A1"] <- "1"
FaunalCluster$cluster[FaunalCluster$cluster == "A2a"] <- "2"
FaunalCluster$cluster[FaunalCluster$cluster == "A2b"] <- "3"
FaunalCluster$cluster[FaunalCluster$cluster == "B1a"] <- "4"
FaunalCluster$cluster[FaunalCluster$cluster == "B1b"] <- "5"
FaunalCluster$cluster[FaunalCluster$cluster == "C1a"] <- "6"
FaunalCluster$cluster[FaunalCluster$cluster == "C1b"] <- "7"
FaunalCluster$cluster[FaunalCluster$cluster == "D1"] <- "8"
FaunalCluster$cluster[FaunalCluster$cluster == "D2a"] <- "9"
FaunalCluster$cluster[FaunalCluster$cluster == "D2b"] <- "10"
FaunalCluster$cluster[FaunalCluster$cluster == "D2c"] <- "11"
FaunalCluster$cluster[FaunalCluster$cluster == "D2d"] <- "12"
FaunalCluster$cluster=as.numeric(FaunalCluster$cluster)

dim(FaunalCluster)# 9619    4

## Take only the coordinates
FaunalCluster2 <- FaunalCluster[,2:3]
head(FaunalCluster2)


#### 18. RANDOM FOREST (BIO): CREATE RASTER STACK FOR ENV PREDICTOR VARIABLES ####

## Load libraries
library(ggplot2)
library(rgdal)
library(maptools)
library(plyr)
library(dismo)
library(raster)

## Identify raster files
list <- list.files(path='DATA/FINAL_RASTERS',pattern=".tif", full.names=TRUE)
list

## Create raster stack
predictors <- stack(list)
predictors

## Simple names for predictor variables
names(predictors)=c("AvCur","Depth","Gravel","Mud","Sal","SPM","Stress","Temp","WOV")#"Sand",
names(predictors)

## Plot raster stack
plot(predictors)


#### 19. RANDOM FOREST (BIO): OVERLAY FAUNAL CLUSTER RECORDS ON ENV RASTERS ####

## Produce map with samples overlaid

## Call library
library(maptools)
data(wrld_simpl)

#And now plot:
## Plot first layer of the RasterStack
plot(predictors, 1)

## Add basic coastline note the "add=TRUE" argument with plot
plot(wrld_simpl, add=TRUE)

## Add sample records using the points function
points(FaunalCluster2, col='blue')


#### 20. RANDOM FOREST (BIO): EXTRACT PREDICTOR VARIABLES FROM RASTER STACK ####

head(FaunalCluster)

## Create a df for predictor variables
sdata <- extract(predictors, FaunalCluster[,2:3])
#View(sdata)
class(sdata)

##Change from matrix to df
sdata=as.data.frame(sdata)

## Combine response and predictor variables into one dataframe
sdata2=cbind(FaunalCluster$Sample,FaunalCluster$cluster,sdata)
colnames(sdata2)[1] <- "Sample"
colnames(sdata2)[2] <- "Cluster"

## Change cols to appropriate type
sdata2$Cluster=as.factor(sdata2$Cluster)
sdata2$Sample=as.character(sdata2$Sample)
sdata2$AvCur=as.numeric(sdata2$AvCur)
#sdata2$Chla=as.numeric(sdata2$Chla)
sdata2$Depth=as.numeric(sdata2$Depth)
sdata2$Gravel=as.numeric(sdata2$Gravel)
sdata2$Mud=as.numeric(sdata2$Mud)
sdata2$Sal=as.numeric(sdata2$Sal)
#sdata2$Sand=as.numeric(sdata2$Sand)
sdata2$SPM=as.numeric(sdata2$SPM)
sdata2$Stress=as.numeric(sdata2$Stress)
sdata2$Temp=as.numeric(sdata2$Temp)
sdata2$WOV=as.numeric(sdata2$WOV)

## First check cols of correct type
str(sdata2)
dim(sdata2)# 9619   11


#### 21. RANDOM FOREST (BIO): MAKE TRAINING AND TESTING SET ####

## Equal splitting code, needs the caTools library;  sdata is the dataframe with cols for response and predictor variables

## Call library
#install.packages("caTools")
library(caTools)

## Vector for response variable
Y <- sdata2[,2]

## Take a random 90% of samples in proportion across the 12 groups
set.seed(2)
msk= sample.split( Y, SplitRatio = 9/10, group = NULL )

## The training set
train = sdata2[ msk,] 
dim(train)#8657   11
#View(train)

## Remove station labels for train
train2 =train[,2:11]
#View(train2)

## The test set
test  = sdata2[ !msk,]
dim(test)#962  11
#View(test)

## Remove station labels for test
test2 =test[,2:11]
#View(test2)
str(test2)

## Check number of observations for train (TRUE) and test (FALSE) sets
print(table(Y, msk)) 

## Check number of samples in train and test sets is equal to total
dim(sdata) # 9619   9
dim(train)+dim(test)# 9619   22


#### 22. RANDOM FOREST (BIO): DO MODELLING ####

## Call library
#install.packages("randomForest")
library(randomForest)

## Model
model <- factor(Cluster) ~ AvCur +  Depth + Gravel+ Mud + Sal + SPM + Stress +Temp + WOV#+ Sand,Chla +

## Run model
rf2 <- randomForest(model, data=train2,na.action=na.exclude)


#### 23. RANDOM FOREST (BIO): EXAMINE HOW VARIABLES ARE INFLUENCING THE MODEL ####

## Produce plots 
varImpPlot(rf2)
preds <- names(rf2$forest$xlevels)

for (i in 1:length(preds)){
  partialPlot(rf2, train2, preds[i], which.class ='1')
  next
}


#### 24. RANDOM FOREST (BIO): EVALUATE THE MODEL PERFORMANCE ####

## Predict cluster group for test set
pred <- predict(rf2, newdata = test2)
table(pred, test2$Cluster)

## We can test the accuracy as follows:
(10+17+18+16+20+21+24+19+59+65+200+102)/ nrow(test2)# 59%

## Matches for nearest cluster neighbour(s), based on dendrogram in Cooper & Barry (2017) fig 3a.
(10+0+0+17+0+18+1+16+6+14+20+21+7+9+24+19+3+2+3+4+59+0+23+5+65+7+2+200+13+13+102)/ nrow(test2)# 71%


#### 25. RANDOM FOREST (BIO): PRODUCE FULL COVERAGE RASTER FOR BIO ####

##Use model to predict cluster group for each raster cell
pr <- predict(predictors, rf2)

## Basic plot
plot(pr, main='Random Forest, regression')#
plot(wrld_simpl, add=TRUE, border='dark grey')

## Update attribute table
pr2 <- ratify(pr)
rat <- levels(pr2)[[1]]
rat$Pixel_Values <- c(1,2,3,4,5,6,7,8,9,10,11,12)
rat$Class_Names <- c("A1","A2a", "A2b", "B1a", "B1b", "C1a", "C1b", "D1", "D2a","D2b", "D2c","D2d")
levels(pr2) <- rat
plot(pr2)


#### 26. RANDOM FOREST (BIO): OUTPUT RASTER AS TIFF ####

writeRaster(pr2,'data/final_model/FaunalCluster/FaunalCluster.tif',overwrite=TRUE,format = "GTiff")


#### 27. FIGURE 2d: BIO CLUSTER CLASSIFICATION PLOT ####

## Clip raster in QGIS
# Bring in vector for StudyArea and raster. Raster>Extraction>Clipper. Save as FaunalClusterClip.tif

library(raster)
library(ggplot2)
library(scales)

## Load raster
fraster = raster('data/final_model/FaunalCluster/FaunalClusterClip.tif')

## convert the raster to points for plotting
fraster.p <- rasterToPoints(fraster)

## Make the points a dataframe for ggplot
fdf <- data.frame(fraster.p)

## Make appropriate column headings
colnames(fdf) <- c("Longitude", "Latitude", "FCluster")
#View(fdf)
str(fdf)

## Change numbers to codes
fdf$FCluster[fdf$FCluster == "1"] <- "A1"
fdf$FCluster[fdf$FCluster == "2"] <- "A2a"
fdf$FCluster[fdf$FCluster == "3"] <- "A2b"
fdf$FCluster[fdf$FCluster == "4"] <- "B1a"
fdf$FCluster[fdf$FCluster == "5"] <- "B1b"
fdf$FCluster[fdf$FCluster == "6"] <- "C1a"
fdf$FCluster[fdf$FCluster == "7"] <- "C1b"
fdf$FCluster[fdf$FCluster == "8"] <- "D1b"
fdf$FCluster[fdf$FCluster == "9"] <- "D2a"
fdf$FCluster[fdf$FCluster == "10"] <- "D2b"
fdf$FCluster[fdf$FCluster == "11"] <- "D2c"
fdf$FCluster[fdf$FCluster == "12"] <- "D2d"

## Make Cluster a factor
fdf$FCluster=as.factor(fdf$FCluster)

#Now make the map
modfc=ggplot(data=fdf, aes(y=Latitude, x=Longitude)) +
  geom_raster(aes(fill=FCluster)) +
  scale_fill_manual(values = c("blue2","cyan1","#05aae1","plum2","darkorchid3","green3","palegreen1","#b40202","red1","darkorange","yellow","#b4b404"),name="Cluster")+
  geom_polygon(data = eubound, aes(x=long, y=lat, group = group),fill="black",colour="black",size=0.15)+
  labs(x="Longitude",y="Latitude")+
  theme_bw(base_size=12)+
  theme(legend.position="bottom")+
  theme(legend.title=element_blank())+
  theme(axis.title.y=element_blank(),
        axis.text.y=element_blank(),axis.title.x=element_blank(),
        axis.text.x=element_blank())+
  coord_quickmap(xlim = c(-10.1, 3.2),ylim=c(48.5, 58.2))+
  theme(legend.text=element_text(size=18))

modfc  

#ggsave("OUTPUTS/modfc.png")


#### 28. RANDOM FOREST (PHY): PREPARE DATA ####

## Start with clusteringresults df
names(allres)

## Select cols for: Sample (1), Longitude_WGS84 (3),Latitude_WGS84 (2), FaunalCluster(16), PhyCluster (17)
SMDdata=allres[,c(1,3,2,16,17)]
#View(SMDdata)#
names(SMDdata)

## Get columns in correct order
PhysicalCluster=SMDdata[,c(1,2,3,5)]
#View(PhysicalCluster)

## Change names of cols
colnames(PhysicalCluster)=c("Sample","lon","lat","cluster")
head(PhysicalCluster)

## Check PhysicalCluster$cluster is factor
str(PhysicalCluster$cluster)#it is
#class(PhysicalCluster)
#View(PhysicalCluster)


#### 29. RANDOM FOREST (PHY): EXTRACT PREDICTOR VARIABLES FROM RASTER STACK ####

## Create a df for predictor variables. Raster stack is object'predictors'
physdata <- extract(predictors, PhysicalCluster[,2:3])
#View(physdata)

##Change from matrix to df
physdata=as.data.frame(physdata)

## Combine response and predictor variables into one dataframe
physdata2=cbind(PhysicalCluster$Sample,PhysicalCluster$cluster,physdata)
colnames(physdata2)[1] <- "Sample"
colnames(physdata2)[2] <- "Cluster"
str(physdata2)

## Change cols to appropriate type
physdata2$Cluster=as.factor(physdata2$Cluster)
physdata2$Sample=as.character(physdata2$Sample)
physdata2$AvCur=as.numeric(physdata2$AvCur)
#physdata2$Chla=as.numeric(physdata2$Chla)
physdata2$Depth=as.numeric(physdata2$Depth)
physdata2$Gravel=as.numeric(physdata2$Gravel)
physdata2$Mud=as.numeric(physdata2$Mud)
physdata2$Sal=as.numeric(physdata2$Sal)
#physdata2$Sand=as.numeric(physdata2$Sand)
physdata2$SPM=as.numeric(physdata2$SPM)
physdata2$Stress=as.numeric(physdata2$Stress)
physdata2$Temp=as.numeric(physdata2$Temp)
physdata2$WOV=as.numeric(physdata2$WOV)

## First check cols of correct type
str(physdata2)


#### 30. RANDOM FOREST (PHY): MAKE TRAINING AND TESTING SET ####

## Equal splitting code, needs the caTools library;  sdata is the dataframe with cols for response and predictor variables

## Call library
library(caTools)

## Vector for response variable
Yphy <- physdata2[,1]

## Take a random 90% of samples in proportion across the 12 groups
set.seed(1)
msk= sample.split( Yphy, SplitRatio = 9/10, group = NULL )

## The training set
trainPhy = physdata2[ msk,] 
#View(trainPhy)

## The test set
testPhy  = physdata2[ !msk,]

## Check number of observations for train (TRUE) and test (FALSE) sets
print(table(Yphy, msk)) 

## Check number of samples in train and test sets is equal to total
dim(physdata2) # 9619   11
dim(trainPhy)+dim(testPhy)# 9619   22


#### 31. RANDOM FOREST (PHY): DO MODELLING ####

## Call library
library(randomForest)

## Model
model <- factor(Cluster) ~ AvCur  + Depth + Gravel+ Mud + Sal + SPM + Stress +Temp + WOV#+ Sand,+ Chla

## Run model
set.seed(1)
rf2phy <- randomForest(model, data=trainPhy,na.action=na.exclude)


#### 32. RANDOM FOREST (PHY): EXAMINE HOW VARIABLES ARE INFLUENCING THE MODEL ####

## Anna'S plots 
varImpPlot(rf2phy)
preds <- names(rf2phy$forest$xlevels)

for (i in 1:length(preds)){
  partialPlot(rf2, train, preds[i], which.class ='1')
  next
}


#### 33. RANDOM FOREST (PHY): EVALUATE THE MODEL PERFORMANCE ####

## Predict cluster group for test set
pred <- predict(rf2phy, newdata = testPhy)
table(pred, testPhy$Cluster)

## We can test the accuracy as follows:
(114+89+49+65+96+22+33+59+126+108+30+64)/ nrow(testPhy)# 89%


#### 34. RANDOM FOREST (PHY): PRODUCE FULL COVERAGE RASTER FOR PHY ####

##Use model to predict cluster group for each raster cell
prPhy <- predict(predictors, rf2phy)

## Basic plot
plot(prPhy, main='Random Forest, regression')#
plot(wrld_simpl, add=TRUE, border='dark grey')

## Update attribute table
pr2phy <- ratify(prPhy)
rat <- levels(pr2phy)[[1]]
rat$Pixel_Values <- c(1,2,3,4,5,6,7,8,9,10,11,12)
rat$Class_Names <- c("1","2", "3", "4", "5", "6", "7", "8", "9","10", "11","12")
levels(pr2phy) <- rat
plot(pr2phy)


#### 35. RANDOM FOREST (PHY): OUTPUT RASTER AS TIFF ####

writeRaster(pr2phy,'data/final_model/PhysicalCluster/PhysicalCluster2.tif',overwrite=TRUE,format = "GTiff")


#### 36. FIGURE 2c: PHY CLUSTER CLASSIFICATION PLOT ####

## Clip raster in QGIS
# Bring in vector for StudyArea and raster. Raster>Extraction>Clipper. Save as FaunalClusterClip.tif

## Load libraries
library(raster)
library(ggplot2)
library(scales)

## Clip raster in QGIS
#Bring in raster and vector for clip extent. Raster>Extraction>Clipper

## Load raster
praster = raster('data/final_model/physicalCluster/PhysicalClusterClip.tif')

#convert the raster to points for plotting
praster.p <- rasterToPoints(praster)

#Make the points a dataframe for ggplot
pdf <- data.frame(praster.p)

#Make appropriate column headings
colnames(pdf) <- c("Longitude", "Latitude", "PCluster")
#View(fdf)
str(pdf)

## Make Cluster a factor
pdf$PCluster=as.factor(pdf$PCluster)

#Now make the map
modpc=ggplot(data=pdf, aes(y=Latitude, x=Longitude)) +
  geom_raster(aes(fill=PCluster)) +
  scale_fill_manual(values = c("#F8766D","#DE8C00","#B79F00","#7CAE00","#00BA38","#00C08B","#00BFC4","#00B4F0","#619CFF","#C77CFF","#F564E3","#FF64B0"),name="Cluster")+#remove this row to accept automatic colours+
  geom_polygon(data = eubound, aes(x=long, y=lat, group = group),fill="black",colour="black",size=0.15)+
  labs(x="Longitude",y="Latitude")+
  theme_bw(base_size=12)+
  theme(legend.position="bottom")+
  theme(legend.title=element_blank())+
  theme(axis.title.y=element_blank(),
        axis.text.y=element_blank(),axis.title.x=element_blank(),
        axis.text.x=element_blank())+
  coord_quickmap(xlim = c(-10.1, 3.2),ylim=c(48.5, 58.2))+
  theme(legend.text=element_text(size=18))

modpc  

#ggsave("OUTPUTS/modpc.png")

## To see colours
#scales::show_col(scales::hue_pal()(12)) 


#### 37. FIGURE 2b: EUNIS 4 PLOT ####

## Basic plot
plot(EUNIS)

## Include only the EUNIS sediment classes
EUNISL4.df=subset(EUNIS.df,EUNIS.df$EUNIScomb=="A5.13"|
                    EUNIS.df$EUNIScomb=="A5.14"|
                    EUNIS.df$EUNIScomb=="A5.15"|
                    EUNIS.df$EUNIScomb=="A5.23 or A5.24"|
                    EUNIS.df$EUNIScomb=="A5.25 or A5.26"|
                    EUNIS.df$EUNIScomb=="A5.27"|
                    EUNIS.df$EUNIScomb=="A5.33"|
                    EUNIS.df$EUNIScomb=="A5.33 or A5.34"|
                    EUNIS.df$EUNIScomb=="A5.34"|
                    EUNIS.df$EUNIScomb=="A5.35"|
                    EUNIS.df$EUNIScomb=="A5.35 or A5.36"|
                    EUNIS.df$EUNIScomb=="A5.36"|
                    EUNIS.df$EUNIScomb=="A5.37"|
                    EUNIS.df$EUNIScomb=="A5.43"|
                    EUNIS.df$EUNIScomb=="A5.44"|
                    EUNIS.df$EUNIScomb=="A5.45")


##Drop unused levels
EUNISL4.df$EUNIScomb=factor(EUNISL4.df$EUNIScomb)

## Check levels for EUNIScomb
levels(EUNISL4.df$EUNIScomb)

## Plot map
pL4=ggplot()+
  geom_polygon(data=EUNISL4.df, aes(x=long, y=lat, group=group,fill=EUNIScomb), size=0.15)+
  geom_polygon(data = eubound, aes(x=long, y=lat, group = group),fill="black",colour="black",size=0.15)+
  labs(x="Longitude",y="Latitude")+ 
  coord_quickmap(xlim = c(-10.1, 3.2),ylim=c(48.5, 58.2))+
  theme_bw(base_size=12)+
  theme(legend.position="bottom")+
  theme(legend.title=element_blank())+
  theme(axis.title.y=element_blank(),axis.text.y=element_blank(),axis.title.x=element_blank(),axis.text.x=element_blank())+
  theme(legend.text=element_text(size=18))

pL4


#### 38. FIGURE 2a: EUNIS 3 PLOT ####

## Make a copy of df 'EUNISL4.df'
EUNISL3.df=EUNISL4.df
#str(EUNISL3.df)

## Collapse level 4 to level 3
EUNISL3.df$EUNIScomb=as.character(EUNISL3.df$EUNIScomb)
EUNISL3.df$EUNIScomb[EUNISL3.df$EUNIScomb == "A5.13"] <- "A5.1"
EUNISL3.df$EUNIScomb[EUNISL3.df$EUNIScomb == "A5.14"] <- "A5.1"
EUNISL3.df$EUNIScomb[EUNISL3.df$EUNIScomb == "A5.15"] <- "A5.1"
EUNISL3.df$EUNIScomb[EUNISL3.df$EUNIScomb == "A5.23 or A5.24"] <- "A5.2"
EUNISL3.df$EUNIScomb[EUNISL3.df$EUNIScomb == "A5.25 or A5.26"] <- "A5.2"
EUNISL3.df$EUNIScomb[EUNISL3.df$EUNIScomb == "A5.27"] <- "A5.2"
EUNISL3.df$EUNIScomb[EUNISL3.df$EUNIScomb == "A5.33"] <- "A5.3"
EUNISL3.df$EUNIScomb[EUNISL3.df$EUNIScomb == "A5.33 or A5.34"] <- "A5.3"
EUNISL3.df$EUNIScomb[EUNISL3.df$EUNIScomb == "A5.34"] <- "A5.3"
EUNISL3.df$EUNIScomb[EUNISL3.df$EUNIScomb == "A5.35"] <- "A5.3"
EUNISL3.df$EUNIScomb[EUNISL3.df$EUNIScomb == "A5.35 or A5.36"] <- "A5.3"
EUNISL3.df$EUNIScomb[EUNISL3.df$EUNIScomb == "A5.36"] <- "A5.3"
EUNISL3.df$EUNIScomb[EUNISL3.df$EUNIScomb == "A5.37"] <- "A5.3"
EUNISL3.df$EUNIScomb[EUNISL3.df$EUNIScomb == "A5.43"] <- "A5.4"
EUNISL3.df$EUNIScomb[EUNISL3.df$EUNIScomb == "A5.44"] <- "A5.4"
EUNISL3.df$EUNIScomb[EUNISL3.df$EUNIScomb == "A5.45"] <- "A5.4"

## Drop unwanted levels
EUNISL3.df$EUNIScomb=factor(EUNISL3.df$EUNIScomb)

## Plot map
pL3=ggplot()+
  geom_polygon(data=EUNISL3.df, aes(x=long, y=lat, group=group,fill=EUNIScomb), size=0.15)+
  geom_polygon(data = eubound, aes(x=long, y=lat, group = group),fill="black",colour="black",size=0.15)+
  labs(x="Longitude",y="Latitude")+ 
  coord_quickmap(xlim = c(-10.1, 3.2),ylim=c(48.5, 58.2))+
  theme_bw(base_size=12)+
  theme(legend.position="bottom")+
  theme(legend.title=element_blank())+
  theme(axis.title.x=element_blank(),axis.text.x=element_blank(),axis.title.y=element_blank(),axis.text.y=element_blank())+
  theme(legend.text=element_text(size=18))

pL3


#### 39. FIGURE 2: COMBINE PLOTS 2a-2d ####

#install.packages("cowplot")
require(cowplot)

hcR=plot_grid(pL3,NULL,pL4,NULL,modpc,NULL,modfc, labels = c(" a)                                  EUNIS 3","","b)                                EUNIS 4",""," c)                                 PHY","    d)                                  BIO"),nrow = 2,rel_widths = c(1,0.05,1,0.05,1,0.05,1),label_size = 26,hjust=0.01,vjust=1,align = 'h')#scale = 0.95

png("OUTPUTS/Figure_2.png",width=45, height=63, units="cm", res=400)#res=600
#tiff("OUTPUTS/Figure_2.tif",width=45, height=63, units="cm", res=250,compression = "lzw")

hcR 
dev.off()


#### 40. READY TEST DATA FOR POWER ANALYSIS ETC ####

## Start with the test data (fauna and other variables)
dim(testdat)# 4123  765
names(testdat)

## Select only the faunal data
data5test=testdat[,2:636]

## Check number of samples
dim(data5test) # 4123  635

## Check df 'data5test' is just the faunal data
names(data5test)# it is

## Create a df 'postest' for Sample, Latitude_WGS84 and Longitude_WGS84 
names(testdat)
postest=testdat[,c(1,739:740)]

## Check names of df 'postest'
names(postest)
#View(postest)

## Calculate Richness and abundance for test samples
library(vegan)
Richness_test = specnumber(data5test) # Species Richness(S)
Abundance_test=rowSums(data5test) # Abundance

## Join richness and abund results to df univar_test
univar_test=rbind(Richness_test,Abundance_test)
#View(univar_test)

## Transpose data
univar_test_t=t(univar_test)
#View(univar_test_t)
dim(univar_test_t)# 4123    2

## Add transposed data to positions etc
univar_testfinal=cbind(postest,univar_test_t)
#View(univar_testfinal)
#dim(postest)#4123    3

## Get positions of test samples. These will be used too extract cluster groups from rasters
postest2=postest[,3:2]
dim(postest2)# 4123    2

## Convert to a spatial data frame
coordinates(postest2) <- ~Longitude_WGS84+Latitude_WGS84

## Define CRS
proj4string(postest2) <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"

## Plot positions
plot(postest2)

## Bring in raster for modelled Faunal cluster
fcraster=raster('data/final_model/FaunalCluster/FaunalCluster.tif')
#plot(fcraster)

## Bring in raster for modelled Physical cluster
pcraster=raster('data/final_model/PhysicalCluster/PhysicalCluster.tif')
#plot(pcraster)

## Compile rasters into raster stack
modelled<-stack(fcraster,pcraster)
plot(modelled)

## Extract modelled Faunal and Physical cluster groups
testmodgrp <- extract(modelled, postest2)

## Extracted data is matrix - change to df
class(testmodgrp)
testmodgrp=as.data.frame(testmodgrp)
#View(testmodgrp)

## Join extracted cluster groups to test data univariate results
testdat6=cbind(univar_testfinal,testmodgrp)
#View(testdat6)

## Change modelled faunal cluster numbers to codes
testdat6$FaunalCluster[testdat6$FaunalCluster == "1"] <- "A1"
testdat6$FaunalCluster[testdat6$FaunalCluster == "2"] <- "A2a"
testdat6$FaunalCluster[testdat6$FaunalCluster == "3"] <- "A2b"
testdat6$FaunalCluster[testdat6$FaunalCluster == "4"] <- "B1a"
testdat6$FaunalCluster[testdat6$FaunalCluster == "5"] <- "B1b"
testdat6$FaunalCluster[testdat6$FaunalCluster == "6"] <- "C1a"
testdat6$FaunalCluster[testdat6$FaunalCluster == "7"] <- "C1b"
testdat6$FaunalCluster[testdat6$FaunalCluster == "8"] <- "D1"
testdat6$FaunalCluster[testdat6$FaunalCluster == "9"] <- "D2a"
testdat6$FaunalCluster[testdat6$FaunalCluster == "10"] <- "D2b"
testdat6$FaunalCluster[testdat6$FaunalCluster == "11"] <- "D2c"
testdat6$FaunalCluster[testdat6$FaunalCluster == "12"] <- "D2d"
#View(testdat6)

## Now get Enis classes from EUNIS Seamap 2016 layer
EUNIST=readOGR("DATA/UKSeaMap2016/UKSeaMap2016_SedimentsDissClip.shp")# Load shapefile

## Check number of samples in test dataset
#test=as.data.frame(postest2)
#View(test)# 4123-this is correct

## Extract Eunis classes
uksm <- over(postest2[1:4123,],EUNIST)
str(uksm)

## Add emply cols for receiving EUNIS classes
testdat6$UKSM2016_Level3=NA
testdat6$UKSM2016_Level4=NA

## Now paste extracted data into relevant columns
testdat6$UKSM2016_Level3=uksm$EUNIScomb
testdat6$UKSM2016_Level4=uksm$EUNIScomb
#View(testdat6)

## Check data type for EUNIS codes
str(testdat6)# factors

## Make a EUNIS level 3 code from EUNIS level 4 code
testdat6$UKSM2016_Level3=as.character(testdat6$UKSM2016_Level3)# 1st change col to type character

## Change codes from level 4 to level 3
testdat6$UKSM2016_Level3[testdat6$UKSM2016_Level3 == "A5.13"] <- "A5.1"
testdat6$UKSM2016_Level3[testdat6$UKSM2016_Level3 == "A5.14"] <- "A5.1"
testdat6$UKSM2016_Level3[testdat6$UKSM2016_Level3 == "A5.15"] <- "A5.1"
testdat6$UKSM2016_Level3[testdat6$UKSM2016_Level3 == "A5.23 or A5.24"] <- "A5.2"
testdat6$UKSM2016_Level3[testdat6$UKSM2016_Level3 == "A5.25 or A5.26"] <- "A5.2"
testdat6$UKSM2016_Level3[testdat6$UKSM2016_Level3 == "A5.27"] <- "A5.2"
testdat6$UKSM2016_Level3[testdat6$UKSM2016_Level3 == "A5.33"] <- "A5.3"
testdat6$UKSM2016_Level3[testdat6$UKSM2016_Level3 == "A5.33 or A5.34"] <- "A5.3"
testdat6$UKSM2016_Level3[testdat6$UKSM2016_Level3 == "A5.34"] <- "A5.3"
testdat6$UKSM2016_Level3[testdat6$UKSM2016_Level3 == "A5.35"] <- "A5.3"
testdat6$UKSM2016_Level3[testdat6$UKSM2016_Level3 == "A5.35 or A5.36"] <- "A5.3"
testdat6$UKSM2016_Level3[testdat6$UKSM2016_Level3 == "A5.36"] <- "A5.3"
testdat6$UKSM2016_Level3[testdat6$UKSM2016_Level3 == "A5.37"] <- "A5.3"
testdat6$UKSM2016_Level3[testdat6$UKSM2016_Level3 == "A5.43"] <- "A5.4"
testdat6$UKSM2016_Level3[testdat6$UKSM2016_Level3 == "A5.44"] <- "A5.4"
testdat6$UKSM2016_Level3[testdat6$UKSM2016_Level3 == "A5.45"] <- "A5.4"

## Remove unwanted records for Eunis level 3 codes (non sediment)
testdat6$UKSM2016_Level3=as.character(testdat6$UKSM2016_Level3)
testdat6$UKSM2016_Level3[testdat6$UKSM2016_Level3 == "A3.1"] <- NA
testdat6$UKSM2016_Level3[testdat6$UKSM2016_Level3 == "A3.2"] <- NA
testdat6$UKSM2016_Level3[testdat6$UKSM2016_Level3 == "A3.3"] <- NA
testdat6$UKSM2016_Level3[testdat6$UKSM2016_Level3 == "A4.1"] <- NA
testdat6$UKSM2016_Level3[testdat6$UKSM2016_Level3 == "A4.12"] <- NA
testdat6$UKSM2016_Level3[testdat6$UKSM2016_Level3 == "A4.2"] <- NA
testdat6$UKSM2016_Level3[testdat6$UKSM2016_Level3 == "A4.27"] <- NA
testdat6$UKSM2016_Level3[testdat6$UKSM2016_Level3 == "A4.3"] <- NA
testdat6$UKSM2016_Level3[testdat6$UKSM2016_Level3 == "A4.33"] <- NA
testdat6$UKSM2016_Level3[testdat6$UKSM2016_Level3 == "A6"] <- NA
testdat6$UKSM2016_Level3[testdat6$UKSM2016_Level3 == "Na"] <- NA

## Change col back to type:factor
testdat6$UKSM2016_Level3=factor(testdat6$UKSM2016_Level3)
levels(testdat6$UKSM2016_Level3)# check levels
#View(testdat6)

## Use predict function to id acutual faunal group of test samples. Train and test datasets must have same families
## First create a dataframe for family names from the train faunal data (df data5)
template=colnames(data6)# this is a vector
template=as.data.frame(template) # this is a df
#View(template)

## Change name of col 1 from 'template' to 'Family'
colnames(template)[1] <- "Family"
dim(template)# 613   1

## Create a version of the test data that only includes families from the train data. Merge two data frames by ID- Merge keeps only the common columns

## Start with sample code and test faunal data 
names(testdat)
data5testnew=testdat[,1:636]
dim(data5testnew)# 4123  636

## Change any nas to zero
data5testnew[is.na(data5testnew)]=0

## Transpose, keeping names
data5testnew_t=setNames(data.frame(t(data5testnew[,-1])),data5testnew[,1])
str(data5testnew_t)

## Make row names 1st col 'Family'
data5testnew_t <- data.frame(Family=rownames(data5testnew_t),data5testnew_t,row.names=NULL)

## Do Merge
testfamall <- merge(template,data5testnew_t,by="Family")
#View(testfamall)
dim(testfamall)# should be 613 (taxa) and how ever many samples you have (in this case 4124)

## Transpose df 'testfamall' so rows are samples and cols are variables. Note this changes object to a matrix.
testfamallt=t(testfamall)
#View(testfamallt)

## Taxa names are in 1st row - need to make these the column headers
colnames(testfamallt) = testfamallt[1, ]
testfamallt = testfamallt[-1, ]# remove 1st row (taxon names)

## Row names in df 'testfamallt' are the Sample codes - these need to be removed 
row.names(testfamallt) <- NULL
#View(testfamallt)

##Change class of object 'testfamallt' from matrix to dataframe
testfamallt2=as.data.frame(testfamallt)
#class(testfamallt2)
#str(testfamallt2)

## Note that data (taxon counts) are factors. Need to convert to numeric
testfamallt2[] <- lapply(testfamallt2, function(x) as.numeric(as.character(x)))
#str(testfamallt2)
names(testfamallt2)

## Get column order right
testfamallt3=testfamallt2[,c(611,1:610,612,613)]
#testfamallt3=testfamallt2

## Transform the data
testfamallt4=testfamallt3^(0.25)

## Identify faunal cluster groups for test samples (predict::Flexclust)
library("flexclust")

## Convert kmeans output to kcca type (this is necessary when using predict function in flexclust
# package).
resultsA=as.kcca(results,datat)

## Check predicted cluster groups are the same as those from clustering of training data set
pred_train <- predict(resultsA)
pred_train # they are

## Confirm sizes of train abd test datasets
dim(datat)# 9619  613
dim(testfamallt2)# 4123  613
names(testfamallt2)
class(testfamallt2)

## Now use predict function to predict cluster groups for test data.
pred_test <- predict(resultsA, newdata=testfamallt4)
pred_test

## Add cluster group from predict to testdata file
testdat7=cbind(testdat6,pred_test)
names(testdat7)

## Change name of col 'testdat7$pred_test' to 'ActualFC'
names(testdat7)[10]<-paste("ActualFC")
str(testdat7)

## Change cluster number to appropriate code
testdat7$ActualFC[testdat7$ActualFC == 11] <- "A1"
testdat7$ActualFC[testdat7$ActualFC == 2]<- "A2a"
testdat7$ActualFC[testdat7$ActualFC == 3] <- "D2b"
testdat7$ActualFC[testdat7$ActualFC == 9]<- "D2c"
testdat7$ActualFC[testdat7$ActualFC == 7] <- "B1b"
testdat7$ActualFC[testdat7$ActualFC == 8] <- "D2d"
testdat7$ActualFC[testdat7$ActualFC == 4] <- "A2b"
testdat7$ActualFC[testdat7$ActualFC == 1] <- "B1a"
testdat7$ActualFC[testdat7$ActualFC == 12] <- "D1"
testdat7$ActualFC[testdat7$ActualFC == 5] <- "C1b"
testdat7$ActualFC[testdat7$ActualFC == 10] <- "C1a"
testdat7$ActualFC[testdat7$ActualFC == 6]<- "D2a"

## Add column for simplified predicted faunal cluster
testdat7$FaunalClusterPred4=NA
testdat7$FaunalClusterPred4=testdat7$FaunalCluster
#View(testdat7)
#str(testdat7)


##Check faunal cluster groups of test samples match broadscale patterns
names(testdat7)

## RS Map of faunal cluster distribution with all test samples
p2t= ggplot()+
  geom_point(data=testdat7,aes(Longitude_WGS84,Latitude_WGS84,col=FaunalCluster),size=0.45,show.legend = TRUE)+#
  geom_polygon(data = euDF2, aes(x=long, y=lat, group = group),fill="white",colour="black",size=0.1)+
  scale_colour_manual(values = c("blue2","cyan1","#05aae1","plum2","darkorchid3","green3","palegreen1","#b40202","red1","darkorange","yellow","#b4b404"),name="Cluster")+
  guides(colour = guide_legend(override.aes = list(size=3)))+ # Change size of legend dots
  coord_map(xlim = c(-10.7, 4),ylim = c(48, 62))+ #set x,y limits of plot
  theme_bw(base_size = 24)+ 
  labs(x="Longitude",y="Latitude")

fig4at=p2t+theme(legend.key.size = unit(1, "cm"))+
  guides(colour = guide_legend(override.aes = list(size=6)))

fig4at

## Change faunal cluster codes to simplified version A,B,C and D
testdat7$FaunalClusterPred4=as.character(testdat7$FaunalClusterPred4)
testdat7$FaunalClusterPred4[testdat7$FaunalClusterPred4 == "D2d"] <- "D"
testdat7$FaunalClusterPred4[testdat7$FaunalClusterPred4 == "B1a"]<- "B"
testdat7$FaunalClusterPred4[testdat7$FaunalClusterPred4 == "D2c"] <- "D"
testdat7$FaunalClusterPred4[testdat7$FaunalClusterPred4 == "D2b"]<- "D"
testdat7$FaunalClusterPred4[testdat7$FaunalClusterPred4 == "A2b"] <- "A"
testdat7$FaunalClusterPred4[testdat7$FaunalClusterPred4 == "C1a"] <- "C"
testdat7$FaunalClusterPred4[testdat7$FaunalClusterPred4 == "A1"] <- "A"
testdat7$FaunalClusterPred4[testdat7$FaunalClusterPred4 == "D2a"] <- "D"
testdat7$FaunalClusterPred4[testdat7$FaunalClusterPred4 == "D1"] <- "D"
testdat7$FaunalClusterPred4[testdat7$FaunalClusterPred4 == "A2a"] <- "A"
testdat7$FaunalClusterPred4[testdat7$FaunalClusterPred4 == "C1b"] <- "C"
testdat7$FaunalClusterPred4[testdat7$FaunalClusterPred4 == "B1b"]<- "B"

testdat7$FaunalClusterPred4=as.factor(testdat7$FaunalClusterPred4)

## Change col names
str(testdat7)
colnames(testdat7)[4] <- "Richness"
colnames(testdat7)[5] <- "Abundance"
colnames(testdat7)[6] <- "FaunalClusterPred12"
colnames(testdat7)[7] <- "PhyClusterPred"
colnames(testdat7)[10] <- "FaunalCluster"

## Change data types as necessary
str(testdat7)
testdat7$sample=as.character(testdat7$sample)
testdat7$FaunalCluster=as.factor(testdat7$FaunalCluster)
testdat7$FaunalClusterPred12=as.factor(testdat7$FaunalClusterPred12)
testdat7$PhyClusterPred=as.factor(testdat7$PhyClusterPred)
str(testdat7)
#View(testdat7)
#names(testdat7)

names(testdat7)
## Get rid of any rows with no data (due to EUNIS classes)
testdat7 = testdat7[complete.cases(testdat7), ]

names(testdat7)
dim(testdat7)# 4000   11
str(testdat7)

testdat7.1<-testdat7[!(testdat7$UKSM2016_Level4=="A5.34"),]
dim(testdat7.1)#3999   11

testdat8<-testdat7.1[!(testdat7.1$UKSM2016_Level4=="A5.43"),]
dim(testdat8)# 3989   11

## Get rid of any unused levels
testdat8 <- droplevels(testdat8)
str(testdat8)
#View(testdat8)

#write.csv(testdat8,"OUTPUTS/testdat8.csv",row.names = F)

testdat9=testdat8[,c(4:9,11)]
#View(testdat9)
str(testdat9)

#write.csv(testdat9,file="outputs/testdat9.csv",row.names = F)

## Change data to long format
library(reshape)
testdat10 <- melt(testdat9, id=c("Richness","Abundance"))
#View(testdat10)
str(testdat10)

## Get summary data for Richness
library(dplyr)
sumres=testdat10%>%group_by(variable,value)%>%summarise(n=n(),mean=mean(Richness,na.rm=TRUE),sd=sd(Richness,na.rm=TRUE),cv= sd(Richness,na.rm=TRUE)/mean(Richness,na.rm=TRUE),delta=0.2*mean(Richness,na.rm=TRUE))
View(sumres)

write.csv(sumres,file="outputs/sumres.csv",row.names = F)

## Get summary data for Abundance
sumresABUN=testdat10%>%group_by(variable,value)%>%summarise(n=n(),mean=mean(Abundance,na.rm=TRUE),sd=sd(Abundance,na.rm=TRUE),cv= sd(Abundance,na.rm=TRUE)/mean(Abundance,na.rm=TRUE),delta=0.2*mean(Abundance,na.rm=TRUE))
View(sumresABUN)

write.csv(sumresABUN,file="outputs/sumresABUN.csv",row.names = F)


#### 41. AUTOCORRELATION: PREPARE DATA ####

## Bring in data. This need to be done with the full dataset (i.e. without reps removed within 50m)
sadata=read.csv("OUTPUTS/tesdat8.csv",header=T,stringsAsFactors=T,check.names=FALSE)
str(sadata)

## Source of surveynames
names(data)

## Change sadata$sample to chr
sadata$sample=as.character(sadata$sample)
str(sadata)

## Take just cols for sample (to be used for matching) and surveyname
names(data)
surveyname=data[,c(1,877)]
names(surveyname)
str(surveyname)
colnames(surveyname)=c("sample","surveyname")

## Join survey names to sadata
library(dplyr)
sadata2=left_join(sadata,surveyname,by="sample")
names(sadata2)

## Take relevant columns (sample	lat	lon	richness	abundance	EUNIS_3	BIO_4	EUNIS_4	PHY	BIO_12	survey)
sadata3=sadata2[,c(1,2,3,4,5,8,11,9,7,6,12)]
names(sadata3)

## Update col names
colnames(sadata3)=c("sample","lat","lon","richness","abundance","EUNIS_3","BIO_4","EUNIS_4","PHY","BIO_12","survey")
View(sadata3)
str(sadata3)

## Change to factors
sadata3$PHY=as.factor(sadata3$PHY)
sadata3$survey=as.factor(sadata3$survey)

## Save data
#write.csv(sadata3,file="OUTPUTS/spacordatav2.csv",row.names = F)


#### 42. AUTOCORRELATION: RICHNESS PLOTS ####
library(cowplot)

## Load data. Note the data used here include all samples (i.e. before exclusion of sample reps within 50m - see step REMOVE 'REPLICATES' TO ADDRESS SPATIAL AUTOCORRELATION)
data=read.csv("OUTPUTS/spacordatav2.csv", header=T,na.strings=c("NA", "-","?","<null>"),check.names=FALSE)
str(data)
#View(data)

## Change PHY to factor
data$PHY=factor(data$PHY,levels = as.character(1:12))

require(ncf)

# Create fata frame for plot output
ncf.cor <- data.frame(Split=character(0),x=numeric(0),y=numeric(0),stringsAsFactors = FALSE)
ncfvals <- data.frame(Factor=character(0), Level=character(0), x=numeric(0),e=numeric(0),y=numeric(0),stringsAsFactors = FALSE)

# Loop through the factor variables
for (i in 6:10){
  
  # Extract factor levels
  splits <- as.character(levels(data[,i]))
  
  # Loop through factor levels
  for (j in 1:length(splits)){
    
    # Extraxt data frame with x,y and richness for each habitat class level in turn 
    td1 <- data[data[,i]==splits[j], c('lon','lat','richness')]
    # Spline correlogram
    splc <- spline.correlog(td1[[1]], 
                            td1[[2]], 
                            td1[[3]],
                            resamp=10,
                            npoints = 500,
                            latlon = TRUE)
    # Add factor, level, and correlogram plot data to output dataframe
    ncf.cor <- rbind(ncf.cor,
                     data.frame(Factor= names(data)[i],
                                Level=splits[j],
                                x=splc$real$predicted$x[1,],
                                y=splc$real$predicted$y[1,],
                                stringsAsFactors = FALSE))
    ncfvals <- rbind(ncfvals,
                     data.frame(Factor= names(data)[i],
                                Level=splits[j],
                                summary(splc)$estimate,
                                stringsAsFactors = FALSE))
  }
  
}


## Spatial autocorrelation plot
require(ggplot2)

plotlist <- vector("list", length = 5)
names(plotlist) <- names(data)[6:10]

for (i in 1:length(plotlist)) {
  
  plotlist[[i]] <-  ggplot(ncf.cor[ncf.cor$Factor==names(plotlist)[i],],aes(x,y)) +
    geom_line(size=0.5) +
    geom_hline(yintercept=0, size=0.5, linetype="solid", color = "grey10") +
    geom_vline(xintercept=50, size=0.5, linetype="dashed", color = "grey10") +
    ylim(c(-1,1)) +
    xlim(c(0,500)) +
    facet_wrap(~Level,ncol = 4) +
    ggtitle(names(plotlist)[i]) +
    #labs(x = "Distance (m)")+
    #labs(y = "Correlation")+
    theme_bw(base_size=14)+ theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                                  panel.background = element_blank(), axis.line = element_line(colour = "black"),axis.title.x=element_blank(),axis.title.y=element_blank())+
    theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.4))

}

## Table with main results. AC distance is the distance at which autocorrelation goes to 0, Max AC is the maximum correlation coefficient at closest distance
#install.packages("tables")
require(tables)
actab <- tabular(Factor(Factor,'Classification')
        *Factor(Level,'Class')*AllObs(ncfvals)~ 
          ('AC distance (m)'= round(x,2)) 
        + ('Max AC'= round(y,2)),data=ncfvals)

actab

## Variogram plots: 4 grp comparison
saplotrich4grp<-plot_grid(plotlist$EUNIS_3,plotlist$BIO_4,NULL,labels = "a)",nrow = 3,label_size = 14, hjust=0.04,rel_heights=c(1,1,0.1))#rel_widths = c(1, 1.15),  
saplotrich4grp

## Variogram plots:12 grp comparison
saplotrich12grp<-plot_grid(plotlist$EUNIS_4,plotlist$PHY,plotlist$BIO_12,NULL, labels = "a)",nrow = 4,label_size = 14, hjust=0.04,rel_heights=c(1,1,1,0.05))#rel_widths = c(1, 1.15),
saplotrich12grp


#### 43. AUTOCORRELATION: ABUNDANCE PLOTS ####

require(ncf)
str(data)
# Change PHY into a factor
data$PHY <- factor(data$PHY,levels=as.character(1:12))

# Create fata frame for plot output
ncf.cor <- data.frame(Split=character(0),x=numeric(0),y=numeric(0),stringsAsFactors = FALSE)
ncfvals <- data.frame(Factor=character(0), Level=character(0), x=numeric(0),e=numeric(0),y=numeric(0),stringsAsFactors = FALSE)

# Loop through the factor variables
for (i in 6:10){
  
  # Extract factor levels
  splits <- as.character(levels(data[,i]))
  
  # Loop through factor levels
  for (j in 1:length(splits)){
    
    # Extraxt data frame with x,y and richness for each habitat class level in turn 
    td1 <- data[data[,i]==splits[j], c('lon','lat','abundance')]
    # Spline correlogram
    splc <- spline.correlog(td1[[1]], 
                            td1[[2]], 
                            td1[[3]],
                            resamp=10,
                            npoints = 500,
                            latlon = TRUE)
    # Add factor, level, and correlogram plot data to output dataframe
    ncf.cor <- rbind(ncf.cor,
                     data.frame(Factor= names(data)[i],
                                Level=splits[j],
                                x=splc$real$predicted$x[1,],
                                y=splc$real$predicted$y[1,],
                                #splc$real$predicted,
                                stringsAsFactors = FALSE))
    ncfvals <- rbind(ncfvals,
                     data.frame(Factor= names(data)[i],
                                Level=splits[j],
                                summary(splc)$estimate,
                                stringsAsFactors = FALSE))
  }
  
}


## Spatial autocorrelation plot
require(ggplot2)

plotlist <- vector("list", length = 5)
names(plotlist) <- names(data)[6:10]

for (i in 1:length(plotlist)) {
  
  plotlist[[i]] <-  ggplot(ncf.cor[ncf.cor$Factor==names(plotlist)[i],],aes(x,y)) +
    geom_line(size=0.5) +
    geom_hline(yintercept=0, size=0.5, linetype="dashed", color = "grey10") +
    geom_vline(xintercept=50, size=0.5, linetype="dashed", color = "grey10") +
    ylim(c(-1,1)) +
    xlim(c(0,500)) +
    facet_wrap(~Level,ncol = 4) +
    ggtitle(names(plotlist)[i]) +
    #labs(x = "Distance (m)")+
    #labs(y = "Correlation")+
    theme_bw(base_size=14)+ theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                                  panel.background = element_blank(), axis.line = element_line(colour = "black"),axis.title.x=element_blank(),axis.title.y=element_blank())+
    theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.4))+
    theme(axis.text.y = element_text(colour ="white"))#,axis.ticks.y=element_blank()
  
}

## Table with main results. AC distance is the distance at which autocorrelation goes to 0, Max AC is the maximum correlation coefficient at closest distance
require(tables)
actab <- tabular(Factor(Factor,'Classification')
                 *Factor(Level,'Class')*AllObs(ncfvals)~ 
                   ('AC distance (m)'= round(x,2)) 
                 + ('Max AC'= round(y,2)),data=ncfvals)
actab

## Variogram plots: 4 grp comparison
saplotabund4grp<-plot_grid(plotlist$EUNIS_3,plotlist$BIO_4,NULL, labels = "b)",nrow = 3,label_size = 14, hjust=0.04,rel_heights=c(1,1,0.1))
saplotabund4grp

## Variogram plots: 12 grp comparison
saplotabund12grp<-plot_grid(plotlist$EUNIS_4,plotlist$PHY,plotlist$BIO_12, labels = "b)",nrow = 4,label_size = 14, hjust=0.04,rel_heights=c(1,1,1,0.05))
saplotabund12grp


#### 44. AUTOCORRELATION: VARIOGRAM PLOTS FOR RICHNESS AND ABUNDANCE TOGETHER ####

## 4 grp comparison
png("OUTPUTS/Figure S3.png",width=30, height=11, units="cm", res=1000)
vario4=plot_grid(NULL,saplotrich4grp,saplotabund4grp, nrow = 1,label_size = 14,rel_widths=c(0.02,1,1))+
  draw_label("Distance (m)", x=0.53, y=  0, vjust=-0.5, angle= 0)+
  draw_label("Correlation", x=  0, y=0.5, vjust= 1.5, angle=90)
vario4
dev.off()

## 12 grp comparison
png("OUTPUTS/Figure S4.png",width=30, height=39, units="cm", res=1000)
vario12=plot_grid(NULL,saplotrich12grp,saplotabund12grp, nrow = 1,label_size = 14,rel_widths=c(0.02,1,1))+
draw_label("Distance (m)", x=0.53, y=  0, vjust=-0.5, angle= 0)+
  draw_label("Correlation", x=  0, y=0.5, vjust= 1.5, angle=90)
vario12
dev.off()


#### 45. FIGURE 1: MAP OF SAMPLE LOCATIONS USED FOR MODELLING/TESTING ####

#View(coordtraindat)
#View(testdat7)

library(shadowtext)
dim(testdat7)
## Plotting
PSam2=ggplot()+
  geom_polygon(data=StudyArea.df, aes(x=long, y=lat, group=group), fill="grey",alpha = 0.3,colour="grey",size=0.35)+#fill=Level4
  geom_polygon(data = eubound, aes(x=long, y=lat, group = group),fill="black",colour="gray 92",size=0.2)+
  geom_point(data=coordtraindat,aes(Longitude_WGS84,Latitude_WGS84),col="blue",size=0.16,show.legend = FALSE)+#train samples
  geom_point(data=testdat7,aes(Longitude_WGS84,Latitude_WGS84),col="red",size=0.16,show.legend = FALSE)+#test samples
  labs(x="Longitude",y="Latitude")+ 
  coord_quickmap(xlim = c(-10.1, 3.2),ylim=c(48.5, 58.2))+
  annotate("text",x=c(-1.4),y=c(53),label=c("UK"),color="white", size=10)+
  annotate("text",x=c(-1,-3.9,-3.6,-6.8),y=c(52.2,56.8,52.55,54.6),label=c("England", "Scotland","Wales", "Northern \n Ireland"),color="dark grey",size=7, angle=c(0,0,90,0))+# Country labels
  geom_shadowtext(aes(x =2.33 , y =54.9 ,label="Dogger \n Bank"), color="black",bg.color='gray 92', size=7)+
  geom_shadowtext(aes(x =-6 , y =51.3 ,label="Celtic \nSea"), color="black",bg.color='gray 92', size=7)+
  geom_shadowtext(aes(x =-1.1, y =50.19 ,label="English Channel"), color="black",bg.color='gray 92', size=7)+
  geom_shadowtext(aes(x =-4.75, y =53.9 ,label="Irish Sea"), color="black",bg.color='gray 92', size=7)+
  geom_shadowtext(aes(x =1, y =56.5 ,label="North Sea"), color="black",bg.color='gray 92', size=7)+
  geom_shadowtext(aes(x =-8.5 , y =50 ,label="Celtic \nShelf"), color="black",bg.color='gray 92', size=7)+
  geom_shadowtext(aes(x =-5.5, y = 59.25,label="Hebridean \n Shelf"), color="black",bg.color='gray 92', size=7)+
  theme_bw(base_size=24)

png("OUTPUTS/Figure_1.png",width = 29.7,height = 36,units = "cm", res = 600,pointsize = 48)
#tiff("OUTPUTS/Figure_1.tif",width=29.7, height=36, units="cm", res=250,compression = "lzw")
PSam2
dev.off()


#### 46. TABLE 4: POWER ANALYSIS (RICHNESS) ####

## Bring in data (if starting from this point)
sumres=read.csv("OUTPUTS/sumres.csv",header=T,stringsAsFactors=T,check.names=FALSE)
str(sumres)
View(sumres)

## See results for each group
sumres[25:28,c(1,2,5,7)]#Eunis 3
sumres[41:44,c(1,2,5,7)]#FC4
sumres[29:40,c(1,2,5,7)]#Eunis 4
sumres[13:24,c(1,2,5,7)]#PC
sumres[1:12,c(1,2,5,7)]#FC

## Now do power test changing values for delta and sd. Results to Table 4.
print(power.t.test(power=0.9,             # Power
                  delta=2.120725,              # Change
                  sd= 11.62357  ,                # Population Standard Deviation
                  sig.level=0.05,         # Significance level
                  type="two.sample",      # Type of study
                  alternative="two.sided"))


#### 47. TABLE 4: POWER ANALYSIS (ABUNDANCE) ####

## Bring in data (if starting from this point)
testdat9=read.csv("outputs/testdat9.csv",header=T,stringsAsFactors=T,check.names=FALSE)
str(testdat9)

## Make PhyClusterPred a factor
testdat9$PhyClusterPred=as.factor(testdat9$PhyClusterPred)

##Start with df 'testdat9'. This is a df with 6010 rows and cols for Abundance, Richness, Cluster groupings, Eunis classes
#View(testdat9)

library(emon)

## Start with
testdat11=testdat9

## Add new col 'logAbund'to df 'testdat11'
testdat11$logAbund=log(testdat11$Abundance)
#View(testdat11)

## Change data to long format
library(reshape)
testdat12 <- melt(testdat11, id=c("Richness","Abundance","logAbund"))
#View(testdat12)
str(testdat12)

## Get some summary data
library(dplyr)

sumresABUN2=testdat12%>%group_by(variable,value)%>%summarise(n=n(),mean=mean(logAbund,na.rm=TRUE),sd=sd(logAbund,na.rm=TRUE))
#View(sumresABUN2)

## See results for each group
sumresABUN2[25:28,c(1,2,4,5)]#Eunis 3
sumresABUN2[41:44,c(1,2,4,5)]#FC4
sumresABUN2[29:40,c(1,2,4,5)]#Eunis 4
sumresABUN2[13:24,c(1,2,4,5)]#PC
sumresABUN2[1:12,c(1,2,4,5)]#FC

## Now do power using the logged data.pars1 is the mean and sd. pars2 is the sd. Change n1 and n1 until result = desired power (i.e. 0.9)  
power.groups(change=20, change.type="M", n1=c(1035), n2=c(1035), pars1=c(2.88, 1.28), pars2=1.28, test='P', distribution="Lognormal")


#### 48. TABLES 2 & 3: RSS/MSS ASSESSMENT (RICHNESS) ####

## Start with df testdat8
names(testdat9)
str(testdat9)

## Subset required data
random=testdat9[,c(1,2,5,7,6,4,3)]
names(random)
dim(random)# 3989    7
colnames(random)= c("Richness","Abundance","EUNISL3.4","Faunal.4","EUNISL4.12","Physical.12","Faunal.12")
#View(random)

## Create vectors for variables
rich = random$Richness
abun = random$Abundance

## Create vectors for factors
eunisl3.4 = random$EUNISL3.4
faunal.4 = random$Faunal.4
eunis12 = random$EUNISL4.12
phys12 = random$Physical.12
faunal12 = random$Faunal.12

faunal.4.weights = (table(faunal.4)-1) / sum(table(faunal.4))
eunisl3.4.weights = (table(eunisl3.4)-1) / sum(table(faunal.4))
eunis12.weights = (table(eunis12)-1) / 3989# 3989 is the number of samples
faunal12.weights = (table(faunal12)-1) / 3989
phys12.weights = (table(phys12)-1) / 3989

## Numbers of samples
table(eunisl3.4)# EUNIS 3
table(faunal.4)#  BIO 4
table(eunis12)#   EUNIS 4
table(phys12)#    PHY
table(faunal12)#  BIO 12

e4.vars = tapply(rich, eunisl3.4, var) 
e4.obs.weights = sum(e4.vars * eunisl3.4.weights)  

f4.vars = tapply(rich, faunal.4, var)
f4.obs.weights = sum(f4.vars * faunal.4.weights)   

e12.vars = tapply(rich, eunis12, var)
e12.obs.weights = sum(e12.vars * eunis12.weights)  

f12.vars = tapply(rich, faunal12, var)
f12.obs.weights = sum(f12.vars * faunal12.weights)  

p12.vars = tapply(rich, phys12, var)
p12.obs.weights = sum(p12.vars * phys12.weights) 

## Within group variances
e4.vars#  EUNIS 3
f4.vars#  BIO 4
e12.vars# EUNIS 4
p12.vars# PHY
f12.vars# BIO 12

## MSS values
e4.obs.weights#  EUNIS 3
f4.obs.weights#  BIO 4
e12.obs.weights# EUNIS 4
p12.obs.weights# PHY
f12.obs.weights# BIO 12

## Doing the randomisation MSS
nreps = 1000
f4.ran.weights = rep(0, nreps)
e4.ran.weights = rep(0, nreps)
e12.ran.weights = rep(0, nreps)
f12.ran.weights = rep(0, nreps)
p12.ran.weights = rep(0, nreps)

for (k in 1:nreps) {
  rfaunal.4 = sample(faunal.4)
  vars = tapply(rich, rfaunal.4, var)
  f4.ran.weights[k] = sum(vars * faunal.4.weights)
  
  reunis.4 = sample(eunisl3.4)
  vars = tapply(rich, reunis.4, var)
  e4.ran.weights[k] = sum(vars * eunisl3.4.weights)
  
  rand = sample(eunis12)
  vars = tapply(rich, rand, var)
  e12.ran.weights[k] = sum(vars * eunis12.weights)
  
  rand = sample(faunal12)
  vars = tapply(rich, rand, var)
  f12.ran.weights[k] = sum(vars * faunal12.weights)
  
  rand = sample(phys12)
  vars = tapply(rich, rand, var)
  p12.ran.weights[k] = sum(vars * phys12.weights)
}


#These will be the same unless the you contrained the number of samples in each cluster.
mean(e4.ran.weights) # EUNIS 3 
mean(f4.ran.weights) # BIO 4      
var(f4.ran.weights)
var(e4.ran.weights) 

## MSS (12 GROUP COMPARISON)
mean(e12.ran.weights)# RANDOM (EUNIS 4)
mean(p12.ran.weights)# RANDOM (PHY)
mean(f12.ran.weights)# RANDOM (BIO 12)

## Calculate randomisation p-value
biggerf = f4.obs.weights >= f4.ran.weights
(pf = (sum(biggerf) + 1)/(nreps + 1))

biggere = e4.obs.weights >= e4.ran.weights
(pe = (sum(biggere) + 1)/(nreps + 1))

f4.obs.weights                        
e4.obs.weights                             
e12.obs.weights                          
f12.obs.weights                           
p12.obs.weights                         

sum(f4.vars)
sum(e4.vars)

## OPTIMAL 4-cluster solution using k-means
nrows = length(rich)
optimal = kmeans(rich, centers=4, nstart=5)
optimal$tot.withinss / nrows

opt.sol = optimal$clus
opt.vars = tapply(rich, opt.sol, var)
opt.weights = (table(opt.sol)-1) / sum(table(opt.sol))

## Number of samples 
table(opt.sol)# OPTIMAL 4

## Within group variances
opt.vars # OPTIMAL 4

## MSS
sum(opt.weights*opt.vars) #OPTIMAL 4

## OPTIMAL 12-cluster solution using k-means
nrows = length(rich)
optimal12 = kmeans(rich, centers=12, nstart=5)
optimal12$tot.withinss / nrows

opt12.sol = optimal12$clus
opt12.vars = tapply(rich, opt12.sol, var)
opt12.weights = (table(opt12.sol)-1) / sum(table(opt12.sol))

## Within group variances
opt12.vars # OPTIMAL 12

## Number of samples 
table(opt12.sol) # OPTIMAL 12

## MSS
sum(opt12.weights*opt12.vars) # OPTIMAL 12


#### 49. TABLES 2 & 3: RSS/MSS ASSESSMENT (ABUNDANCE) ####

e4.vars = tapply(abun, eunisl3.4, var) 
e4.obs.weights = sum(e4.vars * eunisl3.4.weights)  

f4.vars = tapply(abun, faunal.4, var)
f4.obs.weights = sum(f4.vars * faunal.4.weights)   

e12.vars = tapply(abun, eunis12, var)
e12.obs.weights = sum(e12.vars * eunis12.weights)  

f12.vars = tapply(abun, faunal12, var)
f12.obs.weights = sum(f12.vars * faunal12.weights)  

p12.vars = tapply(abun, phys12, var)
p12.obs.weights = sum(p12.vars * phys12.weights) 

## Within group variances
e4.vars#  EUNIS 3
f4.vars#  BIO 4
e12.vars# EUNIS 4
p12.vars# PHY
f12.vars# BIO 12

## MSS values
e4.obs.weights#  EUNIS 3
f4.obs.weights#  BIO 4
e12.obs.weights# EUNIS 4
p12.obs.weights# PHY
f12.obs.weights# BIO 12

## Doing the randomisation MSS
nreps = 1000
f4.ran.weights = rep(0, nreps)
e4.ran.weights = rep(0, nreps)
e12.ran.weights = rep(0, nreps)
f12.ran.weights = rep(0, nreps)
p12.ran.weights = rep(0, nreps)

for (k in 1:nreps) {
  rfaunal.4 = sample(faunal.4)
  vars = tapply(abun, rfaunal.4, var)
  f4.ran.weights[k] = sum(vars * faunal.4.weights)
  
  reunis.4 = sample(eunisl3.4)
  vars = tapply(abun, reunis.4, var)
  e4.ran.weights[k] = sum(vars * eunisl3.4.weights)
  
  rand = sample(eunis12)
  vars = tapply(abun, rand, var)
  e12.ran.weights[k] = sum(vars * eunis12.weights)
  
  rand = sample(faunal12)
  vars = tapply(abun, rand, var)
  f12.ran.weights[k] = sum(vars * faunal12.weights)
  
  rand = sample(phys12)
  vars = tapply(abun, rand, var)
  p12.ran.weights[k] = sum(vars * phys12.weights)
}

## MSS (4 GROUP COMPARISON)
#These will be the same unless the you contrained the number of samples in each cluster.
mean(e4.ran.weights) # EUNIS 3 
mean(f4.ran.weights) # BIO 4      
var(f4.ran.weights)
var(e4.ran.weights) 

## MSS (12 GROUP COMPARISON)
mean(e12.ran.weights)# RANDOM (EUNIS 4)
mean(p12.ran.weights)# RANDOM (PHY)
mean(f12.ran.weights)# RANDOM (BIO 12)

## Calculate randomisation p-value
biggerf = f4.obs.weights >= f4.ran.weights
(pf = (sum(biggerf) + 1)/(nreps + 1))

biggere = e4.obs.weights >= e4.ran.weights
(pe = (sum(biggere) + 1)/(nreps + 1))

f4.obs.weights   
e4.obs.weights  
e12.obs.weights 
f12.obs.weights 
p12.obs.weights 

sum(f4.vars)
sum(e4.vars)

## OPTIMAL 4-cluster solution using k-means
nrows = length(abun)
optimal = kmeans(abun, centers=4, nstart=5)
optimal$tot.withinss / nrows

opt.sol = optimal$clus
opt.vars = tapply(abun, opt.sol, var)
opt.weights = (table(opt.sol)-1) / sum(table(opt.sol))

## Number of samples 
table(opt.sol)# OPTIMAL 4

## Within group variances
opt.vars # OPTIMAL 4

## MSS
sum(opt.weights*opt.vars) #OPTIMAL 4

## OPTIMAL 12-cluster solution using k-means
nrows = length(abun)
optimal12 = kmeans(abun, centers=12, nstart=5)
optimal12$tot.withinss / nrows

opt12.sol = optimal12$clus
opt12.vars = tapply(abun, opt12.sol, var)
opt12.weights = (table(opt12.sol)-1) / sum(table(opt12.sol))

## Within group variances
opt12.vars # OPTIMAL 12

## Number of samples 
table(opt12.sol) # OPTIMAL 12

## MSS 
sum(opt12.weights*opt12.vars) # OPTIMAL 12


#### 50. FIGURE 6a: BAR GRAPHS FOR FAUNAL CLUSTER BY EUNIS LEVEL 3 ####

## Data in df 'testdat8'
names(testdat8)

##  Create a subset of data for histograms (cols UKSM2016_Level3 (5) and FaunalCluster (7)
data4eunis3fchist=testdat8[,c(8,10)]
names(data4eunis3fchist)
str(data4eunis3fchist)

### Use this if you want to exclude any EUNIS classes 
## First get levels for EUNIS 4
data4eunis3fchist$UKSM2016_Level3=as.factor(data4eunis3fchist$UKSM2016_Level3)
levels(data4eunis3fchist$UKSM2016_Level3)

## Exclude non-sediment EUNIS classes
data4eunis3fchist1=subset(data4eunis3fchist,
                          data4eunis3fchist$UKSM2016_Level3=="A5.1"|
                           data4eunis3fchist$UKSM2016_Level3=="A5.2"|
                           data4eunis3fchist$UKSM2016_Level3=="A5.3"|
                           data4eunis3fchist$UKSM2016_Level3=="A5.4")


## Create a histogram for a column of data
bgE3= ggplot(na.omit(data4eunis3fchist1), aes(FaunalCluster,fill=FaunalCluster))+  
  geom_bar()+
  guides(fill=FALSE)+
  scale_fill_manual(values = c("blue2","cyan1","#05aae1","plum2","darkorchid3","green3","palegreen1","#b40202","red1","darkorange","yellow","#b4b404"))+
  theme_bw(base_size=18)+  
  theme(plot.title = element_text(size=18))+
  theme(axis.text.x = element_text(angle=90, hjust=1,vjust=0.5, size=18))+
  theme(axis.text.y = element_text(hjust=1, size=18))+
  theme(axis.text.x=element_blank())+
  labs(y = "Count")+
  scale_x_discrete("Faunal Cluster")+  
  facet_wrap(~UKSM2016_Level3,nrow=1,scales="free_y")+#
  theme(strip.text.x = element_text(size = 18))

bgE3

## Use cow plot to add a white stip to LHS of plot - to give more room for label 'a)'
library(cowplot)
bgE32=plot_grid(NULL,bgE3,nrow = 1,rel_widths = c(0.03,1))#,label_size = 30,hjust=0.04,rel_heights = c(0.6, 1.8))A5.1plot
bgE32


#### 51. FIGURE 6b: BAR GRAPHS FOR FAUNAL CLUSTER BY EUNIS LEVEL 4 ####

## Data in df 'data.new'
names(testdat8)

##  Create a subset of data for hitograms (cols UKSM2016_Level4 (6) and FaunalCluster (7)
data4eunisfchist=testdat8[,9:10]
names(data4eunisfchist)

### Use this if you want to exclude any EUNIS classes 
## First get levels for EUNIS 4
data4eunisfchist$UKSM2016_Level4=as.factor(data4eunisfchist$UKSM2016_Level4)
levels(data4eunisfchist$UKSM2016_Level4)

## Exclude non-sediment EUNIS classes
data4eunisfchist1=subset(data4eunisfchist,data4eunisfchist$UKSM2016_Level4=="A5.13"|
                           data4eunisfchist$UKSM2016_Level4=="A5.14"|
                           data4eunisfchist$UKSM2016_Level4=="A5.15"|
                           data4eunisfchist$UKSM2016_Level4=="A5.23 or A5.24"|
                           data4eunisfchist$UKSM2016_Level4=="A5.25 or A5.26"|
                           data4eunisfchist$UKSM2016_Level4=="A5.27"|
                           data4eunisfchist$UKSM2016_Level4=="A5.33"|
                           data4eunisfchist$UKSM2016_Level4=="A5.33 or A5.34"|
                           data4eunisfchist$UKSM2016_Level4=="A5.34"|
                           data4eunisfchist$UKSM2016_Level4=="A5.35"|
                           data4eunisfchist$UKSM2016_Level4=="A5.35 or A5.36"|
                           data4eunisfchist$UKSM2016_Level4=="A5.36"|
                           data4eunisfchist$UKSM2016_Level4=="A5.37"|
                           data4eunisfchist$UKSM2016_Level4=="A5.43"|
                           data4eunisfchist$UKSM2016_Level4=="A5.44"|
                           data4eunisfchist$UKSM2016_Level4=="A5.45")

## Create a histogram for a column of data
bgE4= ggplot(na.omit(data4eunisfchist1), aes(FaunalCluster,fill=FaunalCluster))+  
  geom_bar()+
  guides(fill=FALSE)+
  scale_fill_manual(values = c("blue2","cyan1","#05aae1","plum2","darkorchid3","green3","palegreen1","#b40202","red1","darkorange","yellow","#b4b404"))+
  theme_bw(base_size=18)+  
  theme(plot.title = element_text(size=18))+
  theme(axis.text.x = element_text(angle=90, hjust=1,vjust=0.5, size=18))+
  theme(axis.text.y = element_text(hjust=1, size=16))+
  labs(y = "Count")+
  scale_x_discrete("Faunal Cluster")+  
  facet_wrap(~UKSM2016_Level4,nrow=4,scales="free_y")+#
  theme(strip.text.x = element_text(size = 18))

bgE4

## Create seperate facet plots for the different Eunis L4 classes grouped by Eunis L3. Start by creating subsets of the data for different groups
A5.1=subset(data4eunisfchist1,data4eunisfchist1$UKSM2016_Level4=="A5.13"|data4eunisfchist1$UKSM2016_Level4=="A5.14"|data4eunisfchist1$UKSM2016_Level4=="A5.15")

A5.2=subset(data4eunisfchist1,data4eunisfchist1$UKSM2016_Level4=="A5.23 or A5.24"|data4eunisfchist1$UKSM2016_Level4=="A5.25 or A5.26"|data4eunisfchist1$UKSM2016_Level4=="A5.27")

A5.3=subset(data4eunisfchist1,data4eunisfchist1$UKSM2016_Level4=="A5.33"|data4eunisfchist1$UKSM2016_Level4=="A5.34"|data4eunisfchist1$UKSM2016_Level4=="A5.35"|data4eunisfchist1$UKSM2016_Level4=="A5.36"|data4eunisfchist1$UKSM2016_Level4=="A5.37")

A5.4=subset(data4eunisfchist1,data4eunisfchist1$UKSM2016_Level4=="A5.43"|data4eunisfchist1$UKSM2016_Level4=="A5.44"|data4eunisfchist1$UKSM2016_Level4=="A5.45")

## Fudge to make sure all Faunal classes are included the level 5.3 plots. Do this by introducing 1 record for each group not represented. These records will be whited out in the plot. Failure to do this will result in the L5.3 plots having missing x classes.
vecUKSM2016_Level4=c("A5.37","A5.37","A5.37","A5.37","A5.37")
vecFaunalCluster=c("A1","A2a","A2b","B1a","C1a")
fudge=data.frame(vecUKSM2016_Level4,vecFaunalCluster)
#View(fudge)
colnames(fudge)=c("UKSM2016_Level4","FaunalCluster")
A5.3v2 <- rbind(A5.3, fudge)

## Produce plots
A5.1plot= ggplot(na.omit(A5.1), aes(FaunalCluster,fill=FaunalCluster))+  
  geom_bar()+
  guides(fill=FALSE)+
  scale_fill_manual(values = c("blue2","cyan1","#05aae1","plum2","darkorchid3","green3","palegreen1","#b40202","red1","darkorange","yellow","#b4b404"))+
  theme_bw(base_size=18)+  
  theme(plot.title = element_text(size=18))+
  theme(axis.text.x = element_text(angle=90, hjust=1,vjust=0.5, size=18))+
  theme(axis.text.y = element_text(hjust=1, size=18))+
  theme(axis.title.x=element_blank())+
  labs(y = "Count")+
  scale_x_discrete("Faunal Cluster")+  
  facet_wrap(~UKSM2016_Level4,nrow=5,scales="free_y")+#
  theme(strip.text.x = element_text(size = 18))
A5.1plot

A5.2plot= ggplot(na.omit(A5.2), aes(FaunalCluster,fill=FaunalCluster))+  
  geom_bar()+
  guides(fill=FALSE)+
  scale_fill_manual(values = c("blue2","cyan1","#05aae1","plum2","darkorchid3","green3","palegreen1","#b40202","red1","darkorange","yellow","#b4b404"))+
  theme_bw(base_size=18)+  
  theme(plot.title = element_text(size=18))+
  theme(axis.text.x = element_text(angle=90, hjust=1,vjust=0.5, size=18))+
  theme(axis.text.y = element_text(hjust=1, size=18))+
  theme(axis.title.x=element_blank())+
  theme(axis.title.y=element_blank())+
  labs(y = "Count")+
  scale_x_discrete("Faunal Cluster")+  
  facet_wrap(~UKSM2016_Level4,nrow=5,scales="free_y")+#
  theme(strip.text.x = element_text(size = 18))
A5.2plot

A5.3plot= ggplot(na.omit(A5.3v2), aes(FaunalCluster,fill=FaunalCluster))+  
    geom_bar()+
  guides(fill=FALSE)+
  scale_fill_manual(values = c("white","white","white","white","darkorchid3","white","palegreen1","#b40202","red1","darkorange","yellow","#b4b404"))+
  theme_bw(base_size=18)+  
  theme(plot.title = element_text(size=18))+
  theme(axis.text.x = element_text(angle=90, hjust=1,vjust=0.5, size=18))+
  theme(axis.text.y = element_text(hjust=1, size=18))+
  theme(axis.title.x=element_blank())+
  theme(axis.title.y=element_blank())+
  labs(y = "Count")+
  scale_x_discrete("Faunal Cluster")+  
  facet_wrap(~UKSM2016_Level4,nrow=5,scales="free_y")+#
  theme(strip.text.x = element_text(size = 18))
A5.3plot

A5.4plot= ggplot(na.omit(A5.4), aes(FaunalCluster,fill=FaunalCluster))+  
  geom_bar()+
  guides(fill=FALSE)+
  scale_fill_manual(values = c("blue2","cyan1","#05aae1","plum2","darkorchid3","green3","palegreen1","#b40202","red1","darkorange","yellow","#b4b404"))+
  theme_bw(base_size=18)+  
  theme(plot.title = element_text(size=18))+
  theme(axis.text.x = element_text(angle=90, hjust=1,vjust=0.5, size=18))+
  theme(axis.text.y = element_text(hjust=1, size=18))+
  theme(axis.title.x=element_blank())+
  theme(axis.title.y=element_blank())+
  labs(y = "Count")+
  scale_x_discrete("Faunal Cluster")+  
  facet_wrap(~UKSM2016_Level4,nrow=5,scales="free_y")+#
  theme(strip.text.x = element_text(size = 18))
A5.4plot

## As plots for A5.1s, A5.2s have 1 less plot we'll need to add in some blnk space (so plots will line up in final plot). Also A5.4s has 2 less plots
final5.1=plot_grid(A5.1plot,NULL,nrow = 2,rel_heights = c(1, 0.31))

##First add in 2 blank plots at bottom of 5.1, 5.2 and 5.4
final5.2=plot_grid(A5.2plot,NULL,nrow = 2,rel_heights = c(1, 0.31))

##First add in 2 blank plots at bottom of 5.1, 5.2 and 5.4
final5.4=plot_grid(A5.4plot,NULL,nrow = 2,rel_heights = c(1, 0.905))

## Now stitch together subplots to create a final level 4 plot
require(cowplot)

EUNIS4Plots=plot_grid(final5.1,final5.2,A5.3plot,final5.4,nrow = 1,rel_heights = c(1,1,1,1),rel_widths = c(1.08,1,1,0.98))

## Add in a white space to LHS of plot 'EUNIS4Plots' to allow more space for label 'b)'
EUNIS4Plots2=plot_grid(NULL,EUNIS4Plots,nrow = 1,rel_widths = c(0.03,1))


#### 52. FIGURE 6c: BAR GRAPHS FOR FAUNAL CLUSTER BY PHY CLUSTER ####

## Data in df 'data.new'
names(testdat8)

##  Create a subset of data for hitograms (cols ModPhy (11) and FaunalCluster (7)
data4phycluschist=testdat8[,c(7,10)]
names(data4phycluschist)
str(data4phycluschist)

## Change order of factors
data4phycluschist$PhyClusterPred = factor(data4phycluschist$PhyClusterPred, levels=c('1','7','2','8','3','9','4','10','5','11','6','12'))

levels(data4phycluschist$PhyClusterPred)

## Fix the factor order
data4phycluschist$PhyClusterPred = factor(data4phycluschist$PhyClusterPred)

## Create a histogram for a column of data
bgPhyClus= ggplot(na.omit(data4phycluschist), aes(FaunalCluster,fill=FaunalCluster))+  
  geom_bar()+
  guides(fill=FALSE)+
  scale_fill_manual(values = c("blue2","cyan1","#05aae1","plum2","darkorchid3","green3","palegreen1","#b40202","red1","darkorange","yellow","#b4b404"))+
  theme_bw(base_size=18)+  
  theme(plot.title = element_text(size=16))+
  theme(axis.text.x = element_text(angle=90, hjust=1,vjust=0.5, size=18))+
  theme(axis.text.y = element_text(hjust=1, size=18))+
  labs(y = "Count")+
  scale_x_discrete("Faunal Cluster")+  
  facet_wrap(~PhyClusterPred,nrow=6,scales="free_y")+#
  theme(strip.text.x = element_text(size = 18))

bgPhyClus


#### 53. FIGURE 6d: BAR GRAPHS FOR FAUNAL CLUSTER BY MODELLED FAUNAL CLUSTER ####

## Data in df 'data.new'
names(testdat8)

##  Create a subset of data for hitograms (cols ModRSMP (10) and FaunalCluster (7)
FCbyModFCdata=testdat8[,c(6,10)]
names(FCbyModFCdata)
str(FCbyModFCdata)
#View(FCbyModFCdata)

## Change order of factors
FCbyModFCdata$FaunalClusterPred = factor(FCbyModFCdata$FaunalClusterPred, levels=c('A1','C1b','A2a','D1','A2b','D2a','B1a','D2b','B1b','D2c','C1a','D2d'))

levels(FCbyModFCdata$FaunalClusterPred)

## Fix the factor order
FCbyModFCdata$FaunalClusterPred = factor(FCbyModFCdata$FaunalClusterPred)

## Produce plot
library(ggplot2)
bgRSMPModFC= ggplot(na.omit(FCbyModFCdata), aes(FaunalCluster,fill=FaunalCluster))+  
  geom_bar()+
  guides(fill=FALSE)+
  scale_fill_manual(values = c("blue2","cyan1","#05aae1","plum2","darkorchid3","green3","palegreen1","#b40202","red1","darkorange","yellow","#b4b404"))+
  theme_bw(base_size=18)+  
  theme(plot.title = element_text(size=18))+
  theme(axis.text.x = element_text(angle=90, hjust=1,vjust=0.5, size=18))+
  theme(axis.text.y = element_text(hjust=1, size=18))+
  labs(y = "Count")+
  scale_x_discrete("Faunal Cluster")+  
  facet_wrap(~FaunalClusterPred,nrow=6,scales="free_y")+#nrow=6
  theme(strip.text.x = element_text(size = 18))

bgRSMPModFC


#### 54. FIGURE 5: COMBINE PLOTS 5a-5d ####

require(cowplot)

## Fist create a combined plot for Eunis levels 3 & 4
EUNISPlots=plot_grid(bgE32,EUNIS4Plots2,labels = c("a)","b)"),nrow = 2,label_size = 28,hjust=0.04,rel_heights = c(0.4, 1.8),align='v')#removealign='v' 
EUNISPlots

## Now join the combined Eunis plot with phy and faunal cluster plot
png("OUTPUTS/Figure_5.png",width=90, height=37, units="cm", res=400)#width=70.4,height=126
#tiff("OUTPUTS/Figure_5.tif",width=90, height=37, units="cm", res=250,compression = "lzw")#width=70.4,height=126
plot_grid(EUNISPlots,NULL,bgPhyClus,NULL,bgRSMPModFC,labels = c("","c)","","d)"),nrow = 1,label_size = 28,hjust=-0.6,rel_widths = c(1,0.04,0.55,0.04,0.55))
dev.off()


#### 55. TABLE 4: SIMPER ANALYSIS BY HABITAT CLASS ####

## Test data with sample labels
names(testdat8)

## Create a df for test data: sample and faunal variables
names(testdat)
data5SIM=testdat[,1:636]
names(data5SIM)

## Select samples common to df 'tdata' and df 'data5SIM' 
data4SIM=data5SIM[data5SIM$sample %in% testdat8$sample, ]
#dim(data4SIM)
#View(data4SIM)

# Now get dfs 'data4SIMOrder' and 'tdata' ordered by sample column
data4SIMOrder <- data4SIM[order(data4SIM$sample),]
#View(data4SIMOrder)
#dim(data4SIMOrder)

tdataOrder <- testdat8[order(testdat8$sample),]
#View(tdataOrder)
dim(tdataOrder)# 3989   11

## Now stitch faunal data to factors
SIMdat2=cbind(data4SIMOrder,tdataOrder)
names(SIMdat2)

## Take only required cols (faunal data (2:636), UKSM2016_Level3 (644),FaunalClusterPred4 (647),UKSM2016_Level4 (645),PhyClusterPred (643),FaunalClusterPred12 (642)                                   
                      
SIMdat3=SIMdat2[,c(1:636,644,647,645,643,642)]
names(SIMdat3)


## Savedata for use in Primer
write.csv(SIMdat3, file="OUTPUTS/DataForSimper2.csv",row.names=FALSE)

## Introduce a col between faunal and treatment groups in csv file. Open in Primer and do SIMPER by group following 4th root transformation.


#### 56. FIGURE 3: VARIATION PLOTS (CV RICH, CV ABUN, SIM) ####

## Data from power table (Table 4)
E3BIO4CV=data.frame(h=c("EUNIS 3","EUNIS 3","EUNIS 3","EUNIS 3","BIO","BIO","BIO","BIO"),Richness=c(0.62,0.71,0.53,0.52,0.47,0.36,0.50,0.71),Abundance=c(2.15,3.07,1.12,1.40,1.29,0.67,1.97,2.93),Similarity=c(23.82,19.23,27.06,25.60,33.46,42.80,28.57,18.64))

## Change order of hab classification systems
E3BIO4CV$h = factor(E3BIO4CV$h,c("EUNIS 3","BIO"))

## Call library 'reshape2' - required for melting data
library(reshape2)
library(ggplot2)

## Melt data into suitable form for facetting 
E3BIO4CV2 <- melt(E3BIO4CV,"h") # Produces a df of 3 cols (PhyCluster, variable, value)

## Change Abund label
levels(E3BIO4CV2$variable)
levels(E3BIO4CV2$variable) <- c("Richness", "Abundance","Similarity")

#Create a df for each variable
E3BIORICH=subset(E3BIO4CV2,E3BIO4CV2$variable=="Richness")
E3BIOABUN=subset(E3BIO4CV2,E3BIO4CV2$variable=="Abundance")
E3BIOSIM=subset(E3BIO4CV2,E3BIO4CV2$variable=="Similarity")

## Plotting
E3BIORICHPLOT=ggplot(E3BIORICH, aes(h,value)) + #fill=PhyCluster
  geom_point() + # don't display outliers#outlier.shape=NA
  theme_bw(base_size = 21)+
  facet_wrap(~variable,scales="free_y",ncol = 1)+
  labs(y="CV")

Fig3aRich=E3BIORICHPLOT+theme(legend.position="none",axis.title.x=element_blank())#+labs(x ="Classification")#,axis.title.x=element_blank()

## Plotting
E3BIOABUNPLOT=ggplot(E3BIOABUN, aes(h,value)) + #fill=PhyCluster
  geom_point() + # don't display outliers#outlier.shape=NA
  theme_bw(base_size = 21)+
  facet_wrap(~variable,scales="free_y",ncol = 1)+
  labs(y="CV")

Fig3aAbun=E3BIOABUNPLOT+theme(legend.position="none",axis.title.x=element_blank())

## Plotting
E3BIOSIMPLOT=ggplot(E3BIOSIM, aes(h,value)) + #fill=PhyCluster
  geom_point() + # don't display outliers#outlier.shape=NA
  theme_bw(base_size = 21)+
  facet_wrap(~variable,scales="free_y",ncol = 1)+
  labs(y="%")

Fig3aSim=E3BIOSIMPLOT+theme(legend.position="none")+labs(x ="Classification")#,axis.title.x=element_blank()

require(cowplot)
library(ggplot2)

## Fist create a combined plot for Eunis levels 3 & 4
Fig3Left=plot_grid(Fig3aRich,Fig3aAbun,Fig3aSim,labels = c("a)","b)","c)",""),nrow = 3,label_size = 21,rel_widths=c(1,1,1))#rel_widths 

Fig3Left

## Data from power table
datanreq=data.frame(h=c("EUNIS 4","EUNIS 4","EUNIS 4","EUNIS 4","EUNIS 4","EUNIS 4","EUNIS 4","EUNIS 4","EUNIS 4","EUNIS 4","EUNIS 4","EUNIS 4","PHY","PHY","PHY","PHY","PHY","PHY","PHY","PHY","PHY","PHY","PHY","PHY","BIO","BIO","BIO","BIO","BIO","BIO","BIO","BIO","BIO","BIO","BIO","BIO"),Richness=c(0.86,0.63,0.60,0.74,0.78,0.64,0.55,0.46,1.40,0.49,0.68,0.41,0.70,0.66,0.64,0.42,0.68,0.56,0.62,0.53,0.49,0.49,0.93,1.10,0.42,0.54,0.41,0.35,0.37,0.58,0.42,0.50,0.52,0.35,0.95,0.52),Abundance=c(1.64,2.07,1.97,3.30,3.43,1.63,0.84,0.86,2.41,1.19,1.43,0.95,2.50,1.38,2.60,0.90,2.07,0.86,1.01,0.85,1.26,0.87,1.77,3.06,0.90,1.17,1.48,0.73,0.57,1.60,2.22,2.19,3.28,1.42,3.11,2.71),Similarity=c(21.76,24.21,25.34,21.93,20.18,20.07,24.97,33.94,12.59,30.96,22.37,29.64,22.81,23.81,23.70,31.34,26.12,29.49,25.71,27.17,28.06,33.42,19.11,19.09,41.75,32.91,37.67,42.06,44.36,26.11,33.68,28.38,26.25,32.14,16.30,29.48))

## Change order of hab classification systems
datanreq$h = factor(datanreq$h,c("EUNIS 4","PHY","BIO"))

## Call library 'reshape2' - required for melting data
library(reshape2)
library(ggplot2)

## Melt data into suitable form for facetting 
dat <- melt(datanreq,"h") # Produces a df of 3 cols (PhyCluster, variable, value)

## Change Abund label
levels(dat$variable)
levels(dat$variable) <- c("Richness", "Abundance","Similarity")

#Create a df for each variable
E4PHYBIORICH=subset(dat,dat$variable=="Richness")
E4PHYBIOABUN=subset(dat,dat$variable=="Abundance")
E4PHYBIOSIM=subset(dat,dat$variable=="Similarity")

## Plotting
E4PHYBIORICHPLOT=ggplot(E4PHYBIORICH, aes(x=h,y=value)) + #fill=PhyCluster
  geom_point() + # don't display outliers#outlier.shape=NA
  theme_bw(base_size = 21)+
  facet_wrap(~variable,scales="free_y")+
  labs(y="")

Fig3bRich=E4PHYBIORICHPLOT+theme(legend.position="none",axis.title.x=element_blank())

Fig3bRich


E4PHYBIOABUNPLOT=ggplot(E4PHYBIOABUN, aes(x=h,y=value)) + #fill=PhyCluster
  geom_point() + # don't display outliers#outlier.shape=NA
  theme_bw(base_size = 21)+
  facet_wrap(~variable,scales="free_y")+
  labs(y="")

Fig3bAbun=E4PHYBIOABUNPLOT+theme(legend.position="none",axis.title.x=element_blank())

Fig3bAbun

E4PHYBIOSIMPLOT=ggplot(E4PHYBIOSIM, aes(x=h,y=value)) + #fill=PhyCluster
  geom_point() + # don't display outliers#outlier.shape=NA
  theme_bw(base_size = 21)+
  facet_wrap(~variable,scales="free_y")+
  labs(y="")

Fig3bSim=E4PHYBIOSIMPLOT+theme(legend.position="none")+labs(x ="Classification")

Fig3bSim

require(cowplot)

## Fist create a combined plot for Eunis levels 3 & 4
Fig3Right=plot_grid(Fig3bRich,Fig3bAbun,Fig3bSim,nrow = 3,label_size = 21,rel_widths=c(1,1,1),align = 'v',hjust=-0.9)#rel_widths #labels = c("d)","e)","f)","")

Fig3Right

## Now stitch left and right parts of Figure 3
require(cowplot)

## Fist create a combined plot for Eunis levels 3 & 4
Fig3=plot_grid(Fig3Left,Fig3Right,nrow = 1,label_size = 21,rel_widths=c(1,1))#rel_widths 

png("OUTPUTS/Figure_3.png",width=30, height=40, units="cm", res=600)
#tiff("OUTPUTS/Figure_3.tif",width=30, height=40, units="cm", res=250,compression = "lzw")
Fig3
dev.off()


#### 57. FIGURE 4: POWER PLOTS TOTAL NREQ ####

## Data from power table
E3BIO4=data.frame(h=c("EUNIS 3","BIO"),Richness=c(759,587),Abundance=c(4329,3240))

## Change order of hab classification systems
E3BIO4$h = factor(E3BIO4$h,c("EUNIS 3","BIO"))

## Call library 'reshape2' - required for melting data
library(reshape2)
library(ggplot2)

## Melt data into suitable form for facetting 
E3BIO4LF <- melt(E3BIO4,"h") # Produces a df of 3 cols (PhyCluster, variable, value)

## Change Abund label
levels(E3BIO4LF$variable)
levels(E3BIO4LF$variable) <- c("Richness", "Abundance")

## Plotting
E3BIO4boxplots=ggplot(E3BIO4LF, aes(x=h,y=value)) + #fill=PhyCluster
  geom_bar(stat="identity")+ 
  theme_bw(base_size = 21)+
  facet_wrap(~variable,scales="free_y",nrow = 2)+
  labs(y=expression(Total~italic(n[req])))

Fig3b=E3BIO4boxplots+theme(legend.position="none") +labs(x ="Classification")

Fig3b

## TOTAL N REQ E4 (12) vs PHY (12) vs BIO (12)

## Data from power table
datanreq=data.frame(h=c("EUNIS 4","PHY","BIO"),Richness=c(3375,2907,1717),Abundance=c(13539,11982,9368))

## Change order of hab classification systems
datanreq$h = factor(datanreq$h,c("EUNIS 4","PHY","BIO"))

## Call library 'reshape2' - required for melting data
library(reshape2)
library(ggplot2)

## Melt data into suitable form for facetting 
dat <- melt(datanreq,"h") # Produces a df of 3 cols (PhyCluster, variable, value)

## Change Abund label
levels(dat$variable)
levels(dat$variable) <- c("Richness", "Abundance")

## Plotting
phyboxplots=ggplot(dat, aes(x=h,y=value)) + #fill=PhyCluster
  geom_bar(stat="identity")+
  theme_bw(base_size = 21)+
  facet_wrap(~variable,scales="free_y",nrow = 2)+
  labs(y=expression(italic(n[req])))
 
Fig3d=phyboxplots+theme(legend.position="none",axis.title.y=element_blank()) +labs(x ="Classification")

Fig3d

## Stitch plots together
require(cowplot)

## Fist create a combined plot for Eunis levels 3 & 4
partcd=plot_grid(Fig3b,Fig3d,labels = c("a)","b)"),nrow = 1,label_size = 21,hjust=-0.9,align='v')#rel_widths 

## Now cowplot all parts together
png("OUTPUTS/Figure_4.png",width=34, height=30, units="cm", res=600)
#tiff("OUTPUTS/Figure_4.tif",width=34, height=30, units="cm", res=250,compression = "lzw")
partcd
dev.off()
