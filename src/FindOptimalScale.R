rm(list=ls())
require("fields")
require("graphics")
require("lattice")
require("MASS")
require("stringr")
require("glmnet")
require("base")
require('caret')
setwd('/WorkingDirectory')

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


#Use cross validation to find the optimal single scale from 6 possible scales
k = 8; #use 8 folds for each cross validation
folds = sample(rep(1:k,each=nrow(agg_data)/k),nrow(agg_data)); #randomly split the data into validation sets that do not overlap
scales = c(.25,.5,1,1.5,2,3)*mean_dist_btn_pts; #define scales corresponding to multiples of avg. dist between stations
  
#Contianer to store the mean squared error for each fold at every scale... ie: mse(s,k) will be the mse of the kth fold at scale s
mse = matrix(rep(NA,k*length(scales)),nrow=length(scales),ncol=k); 

setwd('WorkingDirectory/img') #define a directory to save images

#For each scale, perform an 8 fold cross validation. 
for (scale_idx in 1:length(scales))
{
    scale = scales[scale_idx];
	
	#For each fold, predict the annual mean min temperature at each point on a 100x100 evenly spaced grid spanning 
	#the subset of weather stations comprising the training set for the fold. Generate predictions for the fold’s 
	#test set using the built model and compare the predicted values with the actual min temperature of the validation 
	#stations. Record the difference as the MSE.
    for (fold in folds)
	{
      
		#Each station is used as a base point to generate a bump function. 

		#For each location calculate the blend weight (for that location) of every other bump (one centered at each station) 
		wmat <- exp(-msp^2/(2*scale^2 ) ); #The blending weights of the bumps... 
		wmat.df = data.frame(wmat);
		wmat.df$SID = 1:ncol(msp); #add SID column to enable merging with x matrix below

		#Populate the training and validation sets for this fold 
		trainData = agg_data[folds!=fold,];
		testData = agg_data[folds==fold,];

		trainData <- merge(trainData,wmat.df); #add features from wmat.df (merge based on SID) - repeat features for same station
		testData <- merge(testData,wmat.df);
		ndat <- dim( trainData ) [ 1 ]

		#regress to find a model... the features used are the blend weights (stored in cols 5:ncol(trainData) )
		wModel<-glmnet (as.matrix(trainData[,5:ncol(trainData)]) , as.vector ( trainData [ , 4 ] ) , lambda=0)

		#Predict the min temp for the validation set (a subset of the  original base station locations, for which min temp data is known)
		validation_preds = predict.glmnet(wModel, as.matrix(testData[,5:ncol(testData)]), s='lambda.min' ) #the features used are the blend weights (stored in cols 5:ncol(testData) )
		#Record the mse for this scale and fold
		mse[scale_idx,fold] = mean((testData[,4] - validation_preds)^2);

		####################################################
		#Now generate a heatmap for the annual mean min temp 
		#over an evenly spaced grid of 100x100 points.
		####################################################
		#Create a 2D grid of evenly spaced points held in pmat (points matrix) in the region housing all of the training stations
		#1st: determine the min and max boundaries defining the region containing all the training stations
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
		diffij <- function ( i , j ) sqrt ( rowSums ( ( pmat [ i , ]-xmat [ j , ])^2 ) ); #define a helper function
		distsampletopts <- outer( seq_len (10000) , seq_len (dim( xmat ) [ 1 ] ) , diffij ); #squared differences
		#The actual kernels... these will be the features for the new evenly spaced grid points. 
		wmat<-exp(-distsampletopts^2 /(2*scale[1]^2 ) ) #The actual kernels... distsampletopts has squared distances, scale[i] controls the size.

		#predict a min temperature for each of the evenly spaced grid points... 
		preds <-predict.glmnet(wModel, wmat, s='lambda.min' )

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
		png(str_c("Unregularized_scale=",round(scale,0),"_fold=",fold, ".png"))

		image.plot(xvec,yvec, t(zmat),
		   xlab='East' , ylab='North', 
		   col=heat.colors(12),
		   useRaster=TRUE,
		   main=str_c("Unregularized Reg, scale=",round(scale,0)," fold=",fold))

		dev.off()

    }
}

save(list=c("mse"),file="crossfold_mse");

ave_mse = rowMeans(mse);
ave_mse.df = data.frame(scales,ave_mse);
write.csv(ave_mse.df,"average_mse_across_folds_for_each_scale.csv",row.names=F)
View(ave_mse.df)



