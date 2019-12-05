###########################################################################################################################
# Necessary packages: 
library(ggpubr)
library(reshape2)

###########################################################################################################################
# Function for summary statistics, including standard error: 

## Gives count, mean, standard deviation, standard error of the mean, and confidence interval (default 95%).
##   data: a data frame.
##   measurevar: the name of a column that contains the variable to be summariezed
##   groupvars: a vector containing names of columns that contain grouping variables
##   na.rm: a boolean that indicates whether to ignore NA's
##   conf.interval: the percent range of the confidence interval (default is 95%)
summarySE <- function(data=NULL, measurevar, groupvars=NULL, na.rm=FALSE,
                      conf.interval=.95, .drop=TRUE) {
  library(plyr)
  
  # New version of length which can handle NA's: if na.rm==T, don't count them
  length2 <- function (x, na.rm=FALSE) {
    if (na.rm) sum(!is.na(x))
    else       length(x)
  }
  
  # This does the summary. For each group's data frame, return a vector with
  # N, mean, and sd
  datac <- ddply(data, groupvars, .drop=.drop,
                 .fun = function(xx, col) {
                   c(N    = length2(xx[[col]], na.rm=na.rm),
                     mean = mean   (xx[[col]], na.rm=na.rm),
                     sd   = sd     (xx[[col]], na.rm=na.rm)
                   )
                 },
                 measurevar
  )
  
  # Rename the "mean" column    
  datac <- rename(datac, c("mean" = measurevar))
  
  datac$se <- datac$sd #/ sqrt(datac$N)  # Calculate standard error of the mean
  
  # Confidence interval multiplier for standard error
  # Calculate t-statistic for confidence interval: 
  # e.g., if conf.interval is .95, use .975 (above/below), and use df=N-1
  ciMult <- qt(conf.interval/2 + .5, datac$N-1)
  datac$ci <- datac$se * ciMult
  
  return(datac)
}

###########################################################################################################################
# Import evidence file, you will have to adjust file path!
ev<-read.delim('D:/hs-mq/fp102/combined/txt/evidence.txt')
ev<-ev[ev$PEP<0.02,]
ev<-ev[ev$Potential.contaminant!="+",]

###########################################################################################################################
# Per experiment:  Normalize all reporter ion signal in the 1 to 6 cell equivalent channels to the mean signal 
# (across peptides) per cell signal. 
# So essentially, every RI value is divided by mean( mean(RI_1-cell / 1),  mean(RI_2-cell / 2), ..., mean(RI_6-cell / 6) )

for(X in levels(ev$Raw.file)[1:9]){
  
  mean1<-c()
  for(i in 1:6){
    
    mean.t<-mean( ev[ ev$Raw.file%in%X, c(57:62)[i] ]/i , na.rm = T )
    
    mean1<-c(mean1, mean.t )    
  }
  
  ev[ev$Raw.file==X, 57:62] <- ev[ev$Raw.file==X, 57:62] / mean(mean1, na.rm=T)
  
}

for(X in levels(ev$Raw.file)[10:18]){
  
  mean1<-c()
  for(i in 1:6){
    
    mean.t<-mean( ev[ ev$Raw.file%in%X, c(57:62)[i] ]/ c(4,2,3,1,5,6)[i] , na.rm = T )
    
    mean1<-c(mean1, mean.t )    
  }
  
  ev[ev$Raw.file==X, 57:62] <- ev[ev$Raw.file==X, 57:62] / mean(mean1, na.rm=T)
  
}

###########################################################################################################################
# Reorganize the quantitative data into a more convenient data structure: 
ev1<-ev[ev$Raw.file%in%levels(ev$Raw.file)[1:9],]
colnames(ev1)[57:62]<-c(1:6)
ev1.m<-melt(ev1[,c(16,57:62)])

ev2<-ev[ev$Raw.file%in%levels(ev$Raw.file)[10:18],]
colnames(ev2)[57:62]<-c(4,2,3,1,5,6)
ev2.m<-melt(ev2[,c(16,57:62)])

ev.m<-rbind(ev1.m,ev2.m)

ev.m$variable<-as.numeric(as.character(ev.m$variable))

###########################################################################################################################
# Perform total least squares regression on the normalized reporter ions: 
# Below "value" is the normalized reporter ions, "variable" is the cell-equivalent level (1-6)

ev.x<-ev.m[,c("value","variable")]

# Calculate the mean rRI per cell equivalent level per experiment
ev.x3<-aggregate(data = ev.m, value ~ variable + Raw.file, mean)

# Perform TLS regression on the mean rRI per cell equivalent level per experiment across all experiments
xm<-data.matrix(ev.x3)
ev.svd<-svd((xm[,c(1,3)]))

# The slope of TLS regression is given by: 
slopes<-ev.svd$v[1,2]/(-ev.svd$v[2,2])

# Perform TLS regression on the mean rRI per cell equivalent level per experiment so we can calculate the mean slope estimate: 
raw.reps<-data.frame( c( rep(1:3,3), rep(4:6,3), rep(7:9,3), rep(10:12,3), rep(13:15,3), rep(16:18,3) ), 
                      c( rep(1:18,3) )[order( rep(1:18,3))] )

colnames(raw.reps)<-c("raw","ref")

for(i in unique(xm[,2])){
  
  xm.t<-xm[!xm[,2]%in%c(raw.reps$raw[raw.reps$ref%in%i]), ]
  
  ev.svd<-svd((xm.t[,c(1,3)]))
  
  slopes<-c(slopes, ev.svd$v[1,2]/(-ev.svd$v[2,2]) )
  
}

###########################################################################################################################
# Calculate mean and standard error across all replicates: 

tgc <- summarySE(ev.x3, measurevar="value", groupvars=c("variable"))
pd <- position_dodge(0.1) # move them .05 to the left and right

###########################################################################################################################
# Plot! 

ggplot(tgc, aes(x=variable, y=value)) + 
  geom_errorbar(aes(ymin=value-ci, ymax=value+ci), colour="black", width=.1, position=pd) +
  geom_point(position=pd, size=3) +
  ylim(1,6.5) + 
  xlim(1,6.5) +
  ylab("Normalized Reporter Ions\n") + 
  xlab("Cell count") + 
  theme_pubr() + 
  font("xylab", size=20) + 
  scale_x_continuous(breaks=1:6) +
  scale_y_continuous(breaks=1:6) +
  font("xy.text", size=15)+ 
  annotate(geom="text", x=2.5, y=5, label=paste0("Slope =\n", formatC( round(slopes[1],2), format='f', digits=2 ) ),
           color="red", size=8)
