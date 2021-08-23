data('datsoilspc')
library(dplyr)
library(caret)
##RAW SPECTRA datsoilspc_abs<- tibble(clay, silt, sand, TotalCarbon)
parameters<-as.matrix(datsoilspc[,1:4])
spectra<-as.matrix(datsoilspc$spc)    # The spectra is in reflectance
spc_abs<- apply(spectra, 2, log1R) # spectra converted to Absorbance
# Plotting the raw spectra
matplot(wavelengths, t(datsoilspc_abs$spc_abs), type = "l", 
        xlab = "wavelength", ylab = "log 1/R")
### SNV SPECTRA
spectra_snv<-scale(t(spc_abs),center=TRUE,scale=TRUE)
spectra_snv<- t(spectra_snv)
colnames(spectra_snv)<- wavelengths
matplot(as.numeric(colnames(spectra_snv)), t(spectra_snv), type = "l", 
        xlab = "wavelength", ylab = "Log 1/R")
### SNV+DT SPECTRA
library(prospectr)
spectra_snvdt<- detrend(spc_abs, as.numeric(colnames(spectra)))
matplot(as.numeric(colnames(spectra_snvdt)), t(spectra_snvdt), type = "l", 
                   col = "blue", xlab = "wavelength", ylab = "Log 1/R")
### SNV+DT+1D SPECTRA
spectra_1d<-gapDer(spectra_snvdt,m=1,w=21,s=4)
matplot(as.numeric(colnames(spectra_1d)),t(spectra_1d) ,type="l",
        col = "blue", xlab = "wavelength", ylab = "Log 1/R")
### SNV + 2D SPECTRA
spectra_2d<-gapDer(spectra_snv,m=2,w=21,s=4)
matplot(as.numeric(colnames(spectra_2d)), t(spectra_2d) ,type="l",
        col = "blue", xlab = "wavelength", ylab = "Log 1/R")
### Math treatment spectra correlation with Clay
cor_specsnv_clay<-cor(datsoilspc$clay,spectra_snv[,1:ncol(spectra_snv)])
cor_specsnvdt_clay<-cor(datsoilspc$clay,spectra_snvdt[,1:ncol(spectra_snvdt)])
cor_spec1d_clay<-cor(datsoilspc$clay,spectra_1d[,1:ncol(spectra_1d)])
cor_spec2d_clay<-cor(datsoilspc$clay,spectra_2d[,1:ncol(spectra_2d)])
matplot(as.numeric(colnames(spectra_snv)),t(cor_specsnv_clay),
        type = "l",xlab="nm",col = "red", ylab="Correlation",
        ylim = c(-1, 1),xlim = c(350,2500),
        main = "Mathtreatments Correlations with Clay")
par(new = TRUE)
matplot(as.numeric(colnames(spectra_snvdt)),t(cor_specsnvdt_clay),
        type = "l",xlab=" ",col = "brown", ylab=" ",
        main = " ", ylim = c(-1, 1), xlim = c(350,2500))
par(new = TRUE)
matplot(as.numeric(colnames(spectra_1d)),t(cor_spec1d_clay),
        type = "l",xlab=" ",col = "green", ylab=" ",
        main = " ", ylim = c(-1, 1),xlim = c(350,2500))
par(new = TRUE)
matplot(as.numeric(colnames(spectra_2d)),t(cor_spec2d_clay),
        type = "l",xlab=" ",col = "orange", ylab=" ",
        main = " ", ylim = c(-1, 1),xlim = c(350,2500))

### Correlation of 2D spectra with Sand
cor_spec2d_sand<-cor(datsoilspc$sand,spectra_2d[,1:2098])
matplot(as.numeric(colnames(spectra_2d)),t(cor_spec2d_sand),
        type = "l",xlab="nm",col = "red", ylab="Correlation",
        main = "2D Sand Corr. plot", ylim = c(-1, 1),xlim = c(350,2500))
### Correlation of 2D spectra with Sand
cor_spec2d_silt<-cor(datsoilspc$silt,spectra_2d[,1:2098])
matplot(as.numeric(colnames(spectra_2d)),t(cor_spec2d_silt),
        type = "l",xlab="nm",col = "red", ylab="Correlation",
        main = "2D Silt Corr. plot")
### TEXTURE Correlation plots
matplot(as.numeric(colnames(spectra_2d)),t(cor_spec2d_clay),
        type = "l",xlab="nm",col = "red", ylab="Correlation",
        main = "Texture Correlations", ylim = c(-1, 1))
par(new = TRUE)
matplot(as.numeric(colnames(spectra_2d)),t(cor_spec2d_sand),
        type = "l",xlab=" ",col = "orange", ylab=" ",
        main = " ", ylim = c(-1, 1))
par(new = TRUE)
matplot(as.numeric(colnames(spectra_2d)),t(cor_spec2d_silt),
        type = "l",xlab=" ",col = "brown", ylab=" ",
        main = " ", ylim = c(-1, 1))

#### LOOKING TO THE "Y" MATRIX
## HISTOGRAMS
hist(datsoilspc$clay, main = "Clay")
hist(datsoilspc$sand, main = "Sand")
hist(datsoilspc$silt, main = "Silt")
hist(datsoilspc$TotalCarbon, main = "Total Carbon")
## CORRELATION MATRIX
library(corrplot)
cor<-cor(parameters)
corrplot(cor, method = "circle")

#### SORT THE SPECTRA BY CONSTITUENTS
## By Sand
sand_sort <- datsoilspc_abs %>%
    arrange(desc(sand))
matplot(wavelengths, t(sand_sort$spc_abs[1:5, ]), type = "l", 
        xlab = "wavelength", ylab = "Log 1/R")
## By Clay
clay_sort <- datsoilspc_abs %>%
  arrange(desc(clay))
matplot(wavelengths, t(clay_sort$spc_abs[1:5, ]), type = "l", 
        xlab = "wavelength", ylab = "Log 1/R")
## By Silt
silt_sort <- datsoilspc_abs %>%
  arrange(desc(silt))
matplot(wavelengths, t(silt_sort$spc_abs[1:5, ]), type = "l", 
        xlab = "wavelength", ylab = "Log 1/R")
## By Total Carbon
tcarb_sort <- datsoilspc %>%
  arrange(desc(TotalCarbon))
matplot(wavelengths, t(tcarb_sort$spc[1:5, ]), type = "l", 
        xlab = "wavelength", ylab = "Reflectance")

### Principal Components Analysis with CARET
spc_snvdt_pca<-preProcess(spectra_snvdt,
                       method = c("center", "scale","pca"),
                       thresh = 0.95)
PC_scores<-predict(spc_snvdt_pca,spectra_snvdt)
plot(PC_scores[,1],PC_scores[,2],col="blue",
     xlim=c(min(PC_scores[,1]),max(PC_scores[,1])),
     ylim = c(min(PC_scores[,2]),max(PC_scores[,2])),
     xlab = "PC1",ylab = "PC2")
#standard deviation & mean spectra
matplot(wavelengths,spc_snvdt_pca$std, type = "l", col = "blue")
par(new=TRUE)
matplot(wavelengths,spc_snvdt_pca$mean, type = "l", col = "red")
#Loadings
matplot(wavelengths,spc_snvdt_pca$rotation[,5], type = "l", col = "blue")
#Plot projections
matplot(seq(1,5,by=1),t(PC_scores),pch=1, col = "blue",
        xlab = "PCs", ylab = "Projections")


##########################################################################
#################  DEVELOPING THE REGRESSION #############################
##########################################################################
