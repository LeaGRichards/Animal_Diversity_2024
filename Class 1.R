### Animal Diversity and Distribution ###

  ## Class 1 - camera trapping workshop
  ## Reference book : Camera Trapping for Wildlife Research - Francesco Rovero and FridolinZimmermann
  # Tools for organising data : Wildlife Insights, Agouti
    # 1st step = filter out empty pictures (with a high confidence)
      # Done by MegaDetector, an algorythm that filters pictures into human, animal, vehicles and empty
    # 2nd step = train an algorythm
    # 3rd step = AI trained animal identification unsupervised
    # 4th step = double check 

## Setting things up ##

# Install the necessary packages (erase the # before doing so or else it will stay a comment):
  # install.packages(c("chron","reshape","vegan","plotrix","ggplot2","unmarked","AICcmodavg","MuMIn")

# Call the installed packages with the library function : 
  library(chron)
  library(reshape)
  library(vegan)
  library(plotrix)
  library(ggplot2)
  library(unmarked)
  library(AICcmodavg)
  library(MuMIn)
  source("TEAM library 1.8.R")

# Set working dorectory
  setwd("C:/Users/user/Desktop/R)
  
# Load and visualise dataset
  team <- read.csv(file="team_udz_2009.csv", sep=",",h=T,stringsAsFactors=F)
    # Identifies as team the csv table
    # Seperate data by ',' 
    # h = T --> first row is a headline
    # stringsAsFactors=F --> character vectors should NOT be converted to factors
  head(team)
  
# Fixing data formats for analysis
  data <- fix.dta(team)
  head(data)

  unique(data$Genus)
  unique(data$Species)
  unique(data$bin)

  data<- droplevels(data[data$bin!="Homo sapiens", ]) # remove Homo sapiens from data set
  unique(data$bin)
  unique(data$Class)

# Select mammals and drop Homo, and explore
  names(data)
  data<-data[data$Class=="MAMMALIA",]

  unique(data$Sampling.Unit.Name) #58 camera traps
  unique(data$bin) # species of mammals
  unique(data$Camera.Start.Date.and.Time) # sampling period start
  unique(data$Camera.End.Date.and.Time) # sampling period end

## Descriptive analyses ##

# Camera trap days
  camera_days <- cam.days(data,2009.01)
  head(camera_days)
  summary(camera_days[,2:4]) # 31 average days of trapping --> useful to calculate RAI
    # RAI - relative abundance index, number of independent events/ number days * 100
  write.table(camera_days, file="camera_days2009.txt",quote=F, sep="\t",row.names = F)
    # How to save a document on you working directory, can be in .txt or .csv

# Independent events by chosen time interval
  events_hh<-event.sp(dtaframe=data, year=2009.01, thresh=60)  
    # Thresh works in minutes, is the range of time an event becomes independant
    # Ex: an animal stays 30 min in front of camera - within 60 min, will be 1 event
  events_dd<-event.sp(dtaframe=data, year=2009.01, thresh=1440) 
    # 1440 min = 24h

# Saving tables with events by species and camera site as files in the working drectory
  write.table(events_hh, file="events_hh.txt",quote=F, sep="\t")
  write.table(events_dd, file="events_dd.txt",quote=F, sep="\t")

# Cumulative events per species 
  events_hh_species <- colSums(events_hh)
  events_hh_species
  write.table(events_hh_species, file="events_hh_species.txt", quote=F, sep="\t")
  write.table(events_dd_species, file="events_dd_species.txt",quote=F, sep="\t")

# Sampling effort
  effort <- colSums(camera_days["ndays"])
  effort

# RAI, site-specific and for each species
  rai<-events_hh[,2:27]/effort*100
  rai

# Cumulative RAI 
  # One RAI value for each species (rather than for eahc camera and each spp)
  rai_total <- data.frame(events_hh_species, events_hh_species/sum(effort)*100)
  names(rai_total)[1] <- "events"
  names(rai_total)[2] <- "RAI"
  rai_total
  write.table(rai_total, file="rai_total.txt",quote=F, sep="\t")

# Accumulation curve  
  accumulation <- acc.curve(data,2009.01)
  write.table(accumulation, file="accsp_2009.txt",quote=F, sep="\t")

  ggplot(accumulation, aes(x=Camera.trap.days, y=species)) +
   geom_line(aes(y=species-sd), colour="grey50", linetype="dotted") +
   geom_line(aes(y=species+sd), colour="grey50", linetype="dotted") +
   theme_bw() +
   geom_line() +
   xlab("Camera days") +
   ylab("Num. of species") 
      # The black line is the mean
      # The dotted line are the standard errors
      
# Activity pattern of species 
  activity_24h <- events.hours(data)
  activity_24h
  write.table(activity_24h, file="events_24hour_2009.txt",quote=F, sep="\t",row.names = F)

# Example of plotting activity pattern of selected species (3 forest antelope)
  clock<-c(0:23) 
  clock24.plot(activity_24h$Cephalophus.harveyi,clock,show.grid=T,lwd=2,line.col="blue", main="Cephalophus.harveyi",cex.lab=0.5)
  par(mfrow=c(1,3),cex.lab=0.5, cex.axis=0.5)
  clock24.plot(activity_24h$Cephalophus.spadix,clock,show.grid=T,lwd=2,line.col="green", main="Cephalophus.spadix")
  clock24.plot(activity_24h$Cephalophus.harveyi,clock,show.grid=T,lwd=2,line.col="blue", main="Cephalophus.harveyi")
  clock24.plot(activity_24h$Nesotragus.moschatus,clock,show.grid=T,lwd=2,line.col="red", main="Nesotragus.moschatus")
  par(mfrow=c(1,1))
  clock24.plot(activity_24h$Panthera.pardus, clock, show.grid=T, lwd=2, line.col="yellow", main="leopard")

# Naive occupancy
  mat <- f.matrix.creator(data) # list of matrices camera x days for each species
  mat
  naive_occu_2009<-naive(mat) # get naive occupancy for each species
  naive_occu_2009
  write.table(naive_occu_2009, file="naive_occu_2009.txt",quote=F, sep="\t",row.names = F)

  mat_es <- mat[["Cricetomys gambianus"]]
  mat2 <- mat[["Rhynchocyon udzungwensis"]]

# Save mat
  mat <- save(mat, file="detection_matrix.RDATA")
  
  # to load it later on
  # load("detection_matrix.RDATA")
