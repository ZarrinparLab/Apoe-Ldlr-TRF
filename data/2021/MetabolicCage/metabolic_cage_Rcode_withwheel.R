# LOAD ALL FUNCTIONS FIRST!

load_metadata <- function(file_name) {
  # load metadata based on filename
  # file MUST conform to example discussed in R session
  
  mdat <- read.csv(file_name,
                   colClasses = c(Animal="factor"),
                   strip.white = TRUE)
  mdat <- na.omit(mdat)
  return(mdat)
}



load_metaboliccage_data <- function(file_name) {
  # loads and process input data
  # MUST be csv from 1st tab of metabolic cage data
  
  # Read input data
  df <- read.csv(dat_file, na.strings = ".")
  
  # Remove rows where Animal is NA
  df <- df[!is.na(df$Animal),]
  
  # Convert time
  df$DateTime <- as.POSIXct(df$DateTime, tz = "GMT")
  
  return(df)
}

add_new_col <- function(dat, col1, col2, op, newcol) {
  
  # given dat, performs operation between col1 and col2
  
  if(op == "+") {
    
    dat[,newcol] <- dat[,col1] + dat[,col2]
    
  }
  
  if(op == "-") {
    
    dat[,newcol] <- dat[,col1] - dat[,col2]
    
  }
  
  return(dat)
  
}

#create per hour columns for cumulative metrics
subtract_food_water <- function(dat) {
  
  dat1 <- group_by(dat, Animal)
  dat2 <- mutate(dat1, 
                 food_consumption= FoodInA_M - lag(FoodInA_M, default=first(FoodInA_M), order_by = DateTime),
                 water_consumption= WaterInA_M - lag(WaterInA_M, default=first(WaterInA_M), order_by = DateTime),
                 general_activity= AllMeters_M - lag(AllMeters_M, default=first(AllMeters_M), order_by = DateTime), 
                 wheel_activity= WheelMeters_M - lag(WheelMeters_M, default=first(WheelMeters_M), order_by = DateTime))
  
  res <- as.data.frame(dat2)
  res <- res[with(res, order(Animal, DateTime)),]
  return(res)
}

filter_last_X_hours <- function( dat, x_hours = myhours) {
  
  # x_hours is the last x number of hours to keep
  
  # get all unique timepoints from DateTime
  all_times <- sort(unique(df$DateTime)) 
  
  # Define time to keep
  first_time_id <- length(all_times) - (x_hours - 1)
  last_time_id <- length(all_times)
  time_tokeep <- all_times[first_time_id:last_time_id]
  
  #print(length(time_tokeep))  <---------- this prints the numbers in case one wants to check it
  
  res <- dat[dat$DateTime %in% time_tokeep,]
  return(res)
  
}



get_means_bytime <- function(x) {
  # this function takes data from metabolic cage 
  # and returns mean for each hour of day across all days
  
  df1 <- mutate(x, DateTime=hour(DateTime))
  df2 <- group_by(df1, DateTime, Animal)
  df3 <- as.data.frame(summarise_all(df2, mean, na.rm = TRUE))

  return(df3)
}


merge_data2metadata <- function(dat = means_df , metadata = mdat) {
  # combine metadata amd processed data 
  
  df4 <- merge(x = dat, y = metadata, by="Animal")
  df4$Animal <- as.factor(df4$Animal)
  
  return(df4)
}


plot_mcage_data <- function(dat, i) {
  
  p <- ggplot(data = dat, aes(x = DateTime, y = get(i), colour = Group, fill = Group)) +
    stat_summary(fun = mean, geom="line") +
    stat_summary(fun.data = mean_se, geom="ribbon", alpha=0.25, aes(colour=NULL)) +
    ylab(label = i)
  
  return(p)
  
}



####################################
library(dplyr)
library(lubridate)
library(ggplot2)
library(ggpubr)

setwd("~/Desktop/UCSD/ZarrinparLab/ShortPaper/data/MetabolicCage/")
# Setup
mainDir <- "~/Desktop"
mdat_file <- "./athero_metaboliccagedata_metadata.csv"
dat_file <- "./AthroMetabolicCageRun_Macro1.csv"

myhours <- 72
prefix <- "apoe_ldlr_cagedata"
mycolors <- ""

# main calculations

#use function
#mdat <- load_metadata(file_name = mdat_file)

mdat <- read.csv(mdat_file,
                 colClasses = c(Animal="factor"),
                 strip.white = TRUE)
mdat <- na.omit(mdat)
str(mdat)


df <- load_metaboliccage_data(file_name = dat_file)
colnames(df)
df_add1 <- add_new_col(dat=df, col1="AllMeters_R",
                       col2="PedMeters_R", op="-", newcol="fine_motor_activity")
df_add2 <- subtract_food_water(dat=df_add1)
df_add2

#df_x_hours <- filter_last_X_hours(dat=df_add2, x_hours=myhours)
#dont need for this file

means_df <- get_means_bytime(x=df_add2)
str(means_df)

mdat$Animal
str(mdat)

means_df$Animal
str(means_df)
means_df$Animal <- as.factor(means_df$Animal)
str(means_df)

df4 <- merge(x = means_df, y = mdat, by="Animal")
df4$Animal
str(df4)
df4$VO2_M
write.csv(df4, "./merged.csv", row.names = FALSE)

dat_toplot <- df4
str(dat_toplot)

#########

# Plot data

plot_mcage_data(dat=dat_toplot, i="VO2_M")


mycols <-c("general_activity", "fine_motor_activity", "food_consumption", "water_consumption", "VO2_M", "VCO2_M", "RER_M", "WheelMeters_R", "wheel_activity") 

lapply(mycols, function(x) {

  p <- plot_mcage_data(dat=dat_toplot, i = x)
  return(p)
}) -> myplots

#common.legend = TRUE for one universal key, arrange column numbers for height/width of images
ggarrange(plotlist = myplots, ncol=2, nrow=4, labels="AUTO") %>%
  ggexport(filename = paste0("./", prefix, "_mcage_plots.pdf"))

