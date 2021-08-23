# check if the devtools package is already installed
if (!require("devtools")) install.packages("devtools")

# install the soilspec package from GitHub
devtools::install_github("AlexandreWadoux/soilspec")

library(soilspec)
data("datEPO")
data("datsoilspc")
soilC <- datEPO$soilC
spectra0 <- datEPO$spectra0
spectra1 <- datEPO$spectra1
spectra2 <- datEPO$spectra2

wavelengths<-seq(350, 2500,by = 1)
#spectra0: 100 absorbance spectra from air-dried soil, data frame 
#with 100 rows and 2151 columns
matplot(wavelengths, t(spectra0), type = "l", 
        xlab = "wavelength", ylab = "Reflectance")
spectra0_snv<-scale(t(spectra0),center=TRUE,scale=TRUE)
matplot(wavelengths, spectra0_snv, type = "l", 
        xlab = "wavelength", ylab = "Reflectance")
#spectra1: 100 absorbance spectra from wet-soil, data frame 
#with 100 rows and 2151 columns.
matplot(wavelengths, t(spectra1), type = "l", 
        xlab = "wavelength", ylab = "Reflectance")
#spectra2: 100 absorbance spectra from wet-soil after being air-dried 
#for 1 day, data frame 
#with 100 rows and 2151 columns.
matplot(wavelengths, t(spectra2), type = "l", 
        xlab = "wavelength", ylab = "Reflectance")