############# CALCULATE FROST INDEX FOR DIFFERENT PHENOPHASES ################
rm(list=ls())

### packages
library(ncdf4)
library(abind)
library(doParallel)

### data path
tempfiles <- list.files("../swiss_recon/temp/", pattern = "CH_temp_", full.names = T)

### phenology file
dir <- "data/calibrations/bt_run_2023-09-05_class123_2000burnin_18000iter/pheno_recon/"
phenofiles <- list(mprua65d = paste0(dir,"CH_PTT_mprua65d_1763-01-01-2020-12-31_2023-09-05.nc")) ## add more files if needed

th <- 0
daysahead <- 3
pv <- 1

## read pheno
print(names(phenofiles)[pv])
nc <- nc_open(phenofiles[[pv]])
pheno <- ncvar_get(nc,"DOY")
pheno_N <- ncvar_get(nc,"N")
pheno_E <- ncvar_get(nc,"E")
pheno[pheno == 1] <- NA
pheno[,,1] <- NA

newpheno <- pheno - daysahead ### choose days before phenology

## create nc files with number of days after pheno blossom
tmean_after_pheno <- array(NA, dim(pheno))

cores <- 25
registerDoParallel(cores=cores)
i <- 20

for(i in 1:length(tempfiles)){
  print(yr <- i + 1762)
  nc <- nc_open(tempfiles[i])
  tval <- ncvar_get(nc,ifelse(yr < 1961,"temp","TabsD"))[,,1:212]
  temp_0 <- aperm(apply((tval <= th) + 0 ,1:2,"*", c(1:212) ), c(2,3,1)) ## this turns array
  
  # sum of temperature with tmean below 0 °C 30 days after flowering
  temp_af <- mclapply(1:nrow(tval), function(x) {
    temp2d = tval[x,,]
    temp2d[(temp_0[x,,] < newpheno[x,,i])] = NA 
    ll = sapply(1:nrow(temp2d),function(j) if(!all(is.na(temp2d[j,]))){
      if(!is.na(newpheno[x,j,i])){
        cherrind <- newpheno[x,j,i]:212 # select all doys from onset to end of July
        td <- temp2d[j,cherrind]
        sum(td[td <= th], na.rm = T)
      } else {NA}
    } else {NA})
    return(ll)
  }, mc.cores = cores)
  temp_af <- do.call(rbind, temp_af)
  temp_af[!is.na(tval[,,1]) & is.na(temp_af)] <- 0 # set values inside Switzerland to zero again
  tmean_after_pheno[,,i] <- temp_af
  nc_close(nc)
}

###  create new files
filename <- paste(dir,"/CH_frostindex_acc",th,"_",names(phenofiles)[pv],"_1763-01-01-2020-12-31_",Sys.Date(),".nc",sep="")
nx <- length(pheno_E);ny <- length(pheno_N)
lon1 <- ncdim_def("E", "meter", pheno_E);lat2 <- ncdim_def("N", "meter", pheno_N)

dates <- seq(as.Date("1763-01-01"),as.Date("2020-12-31"), by = "day")
time_def <- ncdim_def("time","days since 1763-01-01", which(lubridate::yday(dates) == 1), unlim=TRUE)
var_temp <- ncvar_def("degrees", "accumulated degrees", list(lon1, lat2, time_def), longname="accumulated degrees", -999)
ncnew <- nc_create(filename, list(var_temp),force_v4=T)
ncvar_put(ncnew, var_temp, tmean_after_pheno)

print("files saved")


