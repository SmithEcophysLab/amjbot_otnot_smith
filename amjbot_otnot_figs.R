# amjbot_otnot_figs
## script to create figures for American Journal of Botany "On the Nature of Things" manuscript on photosynthetic acclimation

## load up packages
library(dplyr)
library(zoo)
library(lubridate)
library(plantecophys)
library(R.utils)
library(ggplot2)

## load up functions
source('model_code/optimal_vcmax_R/calc_optimal_vcmax.R')
sourceDirectory('model_code/optimal_vcmax_R/functions')

## load NEON data
HF_temp = read.csv('/Users/nicksmith/Documents/Research/Timescale_acclimation/Sensitivity/NEON_HarvardForest/NEON_temp-air-single/stackedFiles/SAAT_30min.csv')
HF_par = read.csv('/Users/nicksmith/Documents/Research/Timescale_acclimation/Sensitivity/NEON_HarvardForest/NEON_par/stackedFiles/PARPAR_30min.csv')
HF_rh = read.csv('/Users/nicksmith/Documents/Research/Timescale_acclimation/Sensitivity/NEON_HarvardForest/NEON_rel-humidity/stackedFiles/RH_30min.csv')
nrow(HF_temp)
nrow(HF_par)
nrow(HF_rh)

###############################################################################################
## 1. figure showing photosynthesis and costs of acclimating to different time scales
###############################################################################################

## join average data together
HF_temp_mean = HF_temp[, 6:7]
HF_par_mean = HF_par[, 6:7]
HF_rh_mean = HF_rh[, 6:7]
HF_mean = left_join(HF_temp_mean, HF_par_mean)
HF_mean = left_join(HF_mean, HF_rh_mean)

## calculate VPD from T and rh
vpd_from_t <- function(t, rh){ 
	
	svp = 610.7 * 10^(7.5*t/(237.3+t))
	vpd = ((100-rh) / 100)*svp
	vpd #Pa
	
}

HF_mean$VPD = vpd_from_t(HF_mean$tempSingleMean, HF_mean$RHMean) / 1000

## keep only daytime values
HF_mean_day = subset(HF_mean, PARMean > 0)

## aggregate by day
### separate out dates and times
### see: http://www.neonscience.org/dc-convert-date-time-POSIX-r
HF_mean_day$date = as.factor(as.Date(HF_mean_day$endDateTime)) 
### aggregate
HF_group_by_date = group_by(HF_mean_day[, 2:6], date)
HF_mean_date = summarize(HF_group_by_date, 
                         par = mean(PARMean, na.rm = T), 
                         temp = mean(tempSingleMean, na.rm = T), 
                         vpd = mean(VPD, na.rm = T))

## calculate running means
### previous 30 minutes: i = 1
### previous 3 months: i = 4320 (48 * 90 days)
MA <- list()
for (i in 1:90){
	
	temp_MA = rollapply(HF_mean_date$temp, i, mean, fill = NA, align ='right')
	par_MA = rollapply(HF_mean_date$par, i, mean, fill = NA, align ='right')
	vpd_MA = rollapply(HF_mean_date$vpd, i, mean, fill = NA, align ='right')
	
	temporary = cbind(temp_MA, par_MA, vpd_MA)
	colnames(temporary) = c(paste('temp_MA', i, sep = '_'), paste('par_MA', i, sep = '_'), paste('vpd_MA', i, sep = '_'))
	
	MA[[i]] = temporary
	
}

## calculate optimal traits throughout the summer (June to September)
days_interest = 515:636 # June through September
acclimation_values = list()
for (i in 1:90){
	
	vars = calc_optimal_vcmax(cao = 400, 
	                              tg_c = MA[[i]][days_interest, 1], 
	                              paro = MA[[i]][days_interest, 2], 
	                              z = 300, 
	                              vpdo = MA[[i]][days_interest, 3])
	
	vars_list = cbind(vars, paste(i - 1))
	
	acclimation_values[[i]] = vars_list
	
}

## calculate actual photosynthesis
HF_mean_date_2017 = HF_mean_date[days_interest, ]
HF_patm = (acclimation_values[[1]]$patm / 1000)[1]

photosynthesis_list = list()
for (i in 1:90){
	
	photosynthesis = Photosyn(VPD = HF_mean_date_2017$vpd, # vpd on day of interest
	                          Ca = 400, 
	                          PPFD = HF_mean_date_2017$par, # par on day of interest
	                          Tleaf = HF_mean_date_2017$temp, # temp on day of interest
	                          Patm = HF_patm, # patm on day of interest
	                          Jmax = acclimation_values[[i]]$jmax, # run through vcmax and jmax values
	                          Vcmax = acclimation_values[[i]]$vcmax)
	
	photosynthesis_list[[i]] = photosynthesis
	
}

## get seasonal photosynthesis total
photosynthesis_season = c()
for (i in 1:90){
	
	temp = mean(photosynthesis_list[[i]]$ALEAF)
	photosynthesis_season = c(photosynthesis_season, temp)
	
}

plot(photosynthesis_season ~ seq(1, 90, 1))

## calculate the cost of acclimation as the variability in Vcmax
vcmax_var = c()
for (i in 1:90){
	
	var = sd(acclimation_values[[i]]$vcmax)
	
	vcmax_var = c(vcmax_var, var)
	
}

plot(vcmax_var ~ seq(1, 90, 1))

## make a pretty plot with them together
timescale_dataframe <- data.frame(cbind(photosynthesis_season, vcmax_var))
timescale_plot <- ggplot(data = timescale_dataframe, aes(y = photosynthesis_season/photosynthesis_season[1], x = seq(1,90,1))) +
  theme(legend.position = "right", 
        plot.title = element_text(size = rel(2.2)),
        legend.title = element_text(size = rel(1.5)),
        legend.text = element_text(size = rel(1)),
        plot.tag = element_text(size = rel(2)),
        axis.title.y=element_text(size=rel(2.2), colour = 'black'),
        axis.title.x=element_text(size=rel(2.2), colour = 'black'),
        axis.text.x=element_text(size=rel(2), colour = 'black'),
        axis.text.y=element_text(size=rel(2), colour = 'black'),
        panel.background = element_rect(fill = 'white', colour = 'black'),
        panel.grid.major = element_line(colour = "grey")) +
  geom_line(aes(color='black'), linewidth = 2) +
  geom_line(aes(y=vcmax_var/vcmax_var[1], color='red'), linewidth = 2) +
  ylab('Relative seasonal value') +
  xlab('Acclimation timescale (days)') +
  scale_colour_manual(name = NULL, values =c('black'='black','red'='red'), labels = c('Photosynthesis','Cost')) +
  ylim(c(0.2, 1)) +
  xlim(c(0, 90))

###############################################################################################
## 2. figure showing 1 week of diurnal photosynthesis of average and midday acclimation
###############################################################################################

## calculate max values for acclimation
HF_mean_day$date = as.factor(as.Date(HF_mean_day$endDateTime)) 
HF_group_by_date = group_by(HF_mean_day[, 2:6], date)
HF_max_date = summarize(HF_group_by_date, 
                         par = max(PARMean, na.rm = T), # change needed to get acclimation to max PAR
                         temp = mean(tempSingleMean, na.rm = T), 
                         vpd = mean(VPD, na.rm = T))

## calculate running means
### previous 30 minutes: i = 1
### previous 3 months: i = 4320 (48 * 90 days)
MA_max <- list()
for (i in 1:90){
  
  temp_MA = rollapply(HF_mean_date$temp, i, mean, fill = NA, align ='right')
  par_MA = rollapply(HF_mean_date$par, i, mean, fill = NA, align ='right')
  vpd_MA = rollapply(HF_mean_date$vpd, i, mean, fill = NA, align ='right')
  
  temporary = cbind(temp_MA, par_MA, vpd_MA)
  colnames(temporary) = c(paste('temp_MA', i, sep = '_'), paste('par_MA', i, sep = '_'), paste('vpd_MA', i, sep = '_'))
  
  MA_max[[i]] = temporary
  
}

## calculate optimal traits throughout the summer (June to September)
days_interest = 515:636 # June through September
acclimation_values_max = list()
for (i in 1:90){
  
  vars = calc_optimal_vcmax(cao = 400, 
                            tg_c = MA[[i]][days_interest, 1], 
                            paro = MA[[i]][days_interest, 2], 
                            z = 300, 
                            vpdo = MA[[i]][days_interest, 3])
  
  vars_list = cbind(vars, paste(i - 1))
  
  acclimation_values_max[[i]] = vars_list
  
}

## make a dataframe with 30-min HF data that includes acclimated vcmax and jmax values
HF_30min_data <- subset(HF_mean_day, as.Date(date) > as.Date("2017-05-31") & as.Date(date) < as.Date("2017-10-31"))
head(HF_30min_data)
tail(HF_30min_data)

## create data frame to merge in acclimated values
merge_doys <- data.frame(as.factor(as.Date(seq(as.Date("2017-06-01"), as.Date("2017-09-30"), 1))), days_interest)
colnames(merge_doys) <- c('Date', 'doy')

# STELL TO DO: MERGE THE DATA FRAMES AND MAKE DIURNAL PLOTS, ALSO NEED TO MAKE BETA PLOT


# calculate seasonal photosynthesis assuming LAI = 6, seconds to season = 60 * 60 * 24 * (636 - 515) = 10454400, and 1 µmol CO2 = 1e-6 mol, and 12 g C / mole CO2
photosynthesis_season_gc = photosynthesis_season * 10454400 * 1e-6 * 12

par(mfrow = c(1, 1), mar = c(8, 8, 1, 8))
plot(photosynthesis_season_gc ~ seq(1, 90, 1), 
     type = 'l', lwd = 6, lty = 1, ylim = c(600, 800), xlim = c(1, 90), yaxt = 'n', xaxt = 'n', ylab = '', xlab = '')
axis(2, seq(200, 325, 25), cex.axis = 1.5, las = 2)
par(new = T)
plot(vcmax_var ~ seq(1, 90, 1), 
     type = 'l', lwd = 6, lty = 2, col = 'red', ylim = c(4, 20), xlim = c(1, 90), yaxt = 'n', xaxt = 'n', ylab = '', xlab = '')
axis(1, seq(0, 90, 10), cex.axis = 1.5)
axis(4, seq(2, 10, 2), cex.axis = 1.5, las = 2)

mtext(side = 1, 'Acclimation timescale (days)', cex = 2, line = 4)
mtext(side = 2, expression('Leaf C assimilation (gC m'^'-2'*' yr'^'-1'*')'), cex = 2, line = 4)
mtext(side = 4, expression('Photosynthetic cost (µmol m'^'-2'*' s'^'-1'*')'), cex = 2, line = 4)

legend('topright', c('Assimilation', 'Cost'), col = c('black', 'red'), lwd = 6, lty = c(1, 2), cex = 2)











