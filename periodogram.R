calcium <- read.csv2(file.choose()) # select file

calcium <- read.csv2('20x20.csv') # 24x24 automatic ROI grid
calcium <- read.csv2('fluo4d.csv') # 40 cells are manually selected

# periodogram

par(mfrow=c(1,1), mar=c(2,10,1,1))

 

?analyze.wavelet

# gray

wt.image(my.w, n.levels = 250,
         legend.params = list(lab = "wavelet power levels", label.digits = 2),
         color.palette = "gray((n.levels):1/n.levels)",
         col.ridge = "blue",
         main = "Wavelet power spectrum")