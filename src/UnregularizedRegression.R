rm(list=ls())
require("fields")
require("graphics")
require("lattice")
require("MASS")
require("stringr")
require("glmnet")
require("base")
require('caret')
setwd('WorkingDirectory');
#-------------------------------------------------------------------------
# QUESTION 1.1
#-------------------------------------------------------------------------
LocData = read.csv("Locations.txt",header=TRUE, sep = " ");
TempData = read.csv("Oregon_Met_Data.txt",header=TRUE, sep = " ")
LocData = LocData[,c(1,7:8)] 	#only concerned with the UTM units for lattitude and longitude
TempData = TempData[,c(1,5,6)]	#only concerned with the min temperature and the day of the year

#merge the data, so that we have the cols for location as well as min temp... 
#note: there are only 112 locations, so they will be repeated
data = merge(TempData,LocData); 
#remove any missing or impossible temperatures
data = data[data$Tmin_deg_C < 200 & !is.na(data$Tmin_deg_C),] 
full_data = data; 

#Smoothing will be performed on a processed data set which is obtained by averaging the 
#reported T_min at each location over all times (i.e one interpolate). 
#agg_data will hold this processed data
agg_data = aggregate(data$Tmin_deg_C,by=list(data$SID,data$East_UTM,data$North_UTM),mean)
names(agg_data) = c("SID","East_UTM","North_UTM","Tmin_deg_C");
rm(data);

#calculate the distance between every pair of stations
xmat <- as.matrix( LocData [, c(2,3)] );	#the latitude and longitude of each station
#Calculate the distance between each pair of stations
spaces <- dist ( xmat , method = "euclidean" , diag = FALSE, upper = FALSE);  
msp <- as.matrix( spaces ); #square matrix describes dist between every pair
mean_dist_btn_pts = mean(msp);
  
setwd('WorkingDirectory/img')

#Only one scale... set it to be the mean distance between points. 
scale = mean_dist_btn_pts;

#Each station is used as a base point to generate a bump function. 

#For each location calculate the blend weight (for that location) of every other bump (one centered at each station) 
wmat <- exp(-msp^2/(2*scale^2 ) ); #The blending weights of the bumps... 
wmat.df = data.frame(wmat);
wmat.df$SID = 1:ncol(msp); #add SID column to enable merging with x matrix below

#because we aren't cross validating, use all of the processed data to train the model.
trainData = agg_data;
trainData <- merge(trainData,wmat.df); #add features from wmat.df (merge based on SID) - repeat features for same station
ndat <- dim( trainData ) [ 1 ]

#regress to find a model... the features used are the blend weights (stored in cols 5:ncol(trainData) )
wmod<-glmnet (as.matrix(trainData[,5:ncol(trainData)]) , as.vector ( trainData [ , 4 ] ) , lambda=0)

####################################################
#Now generate a heatmap for the annual mean min temp 
#over an evenly spaced grid of 100x100 points.
####################################################
#Create a 2D grid of evenly spaced points held in pmat (points matrix) in the region housing all of the stations
#1st: determine the min and max boundaries defining the region containing all the stations
xmin<-min( xmat [ , 1 ] )
xmax<-max( xmat [ , 1 ] )
ymin<-min( xmat [ , 2 ] )
ymax<-max( xmat [ , 2 ] )
#2nd: split that region up into 100x100 evenly spaced grid
xvec<-seq ( xmin , xmax , length=100)
yvec<-seq ( ymin , ymax , length=100)
#3rd: populate pmat as a 10000x2 matrix where each row is one of the evenly spaced coordinates 
#(ie: the cols are latitude and longitude) and the rows are Major Ordered representing the 100x100 matrix.
pmat<-matrix (0 , nrow=100*100 , ncol=2) #pmat: a grid of values evenly spaced in longitude and latitude
ptr<-1
for ( i in 1: 100)
{
	for ( j in 1: 100)
	{
		pmat [ ptr , 1 ]<-xvec [ i ]
		pmat [ ptr , 2 ]<-yvec [ j ]
		ptr<-ptr+1
	}
}

#Create gaussian kernels for the even grid of base points... 
diffij <- function ( i , j ) sqrt ( rowSums ( ( pmat [ i , ]-xmat [ j , ])^2 ) ) #define a helper function
distsampletopts <- outer( seq_len (10000) , seq_len (dim( xmat ) [ 1 ] ) , diffij ) #squared differences
#The actual kernels... these will be the features for the new evenly spaced grid points. 
wmat<-exp(-distsampletopts^2 /(2*scale^2 ) ) #The actual kernels... distsampletopts has squared distances, scale controls the size.

#predict a min temperature for each of the evenly spaced grid points... 
preds <-predict.glmnet(wmod, wmat, s='lambda.min' )
  
pred_means = rowMeans(preds,na.rm=T);

#Generate a map of min temperature.
zmat<-matrix (0 , nrow=100, ncol=100)  
ptr<-1
for ( i in 1: 100)
{ 
	for ( j in 1: 100)
	{
		zmat [i,j ]<-pred_means [ ptr ]
		ptr<-ptr+1
	}
}

#Save the map to an image. Think of the final image as a heightmap or heatmap 
#and the z coordinate is the predicted value (min temp)
png(str_c("Unregularized_scale=",round(scale,0),"_final.png"))

image.plot(xvec,yvec, t(zmat),
     xlab='East' , ylab='North', 
	 col=heat.colors(12),
     useRaster=TRUE,
	 main=str_c("part a scale=",round(scale,0)))
    
dev.off()



