log1R<- function(x) {
  return(log(1/x))
} 
spectra_abs<- apply(spectra, 2, log1R)