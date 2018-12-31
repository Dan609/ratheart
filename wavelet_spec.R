library(WaveletComp)
?WaveletComp
# insert ROI

x <- ts(calcium$ROI416, start = calcium$Time[1])
y <- ts(calcium$ROI415, start = calcium$Time[1])

# The wavelet transform of x is computed as follows:
par(mfrow=c(1,1))
my.data <- data.frame(x = x)
my.w <- analyze.wavelet(my.data, "x",
                        loess.span = 0,
                        dt = 1, 
                        dj = 1/250,
                        lowerPeriod = 2,
                        upperPeriod = 128,
                        make.pval = TRUE, n.sim = 10)

# Plot the wavelet power spectrum, a series with a constant period

wt.image(my.w, color.key = "quantile", n.levels = 250, 
         legend.params = list(lab = "wavelet power levels", mar = 4.7))

# A series with a variable period
par(mfrow=c(2,1))
wt.image(my.w, n.levels = 250,
         legend.params = list(lab = "wavelet power levels"))

# The reconstructed series

my.rec <- reconstruct(my.w)
x.rec <- my.rec$series$x.r   #x: name of original series

# A series with two periods

my.data <- data.frame(x = x)
my.w <- analyze.wavelet(my.data, "x",
                        loess.span = 0,
                        dt = 1, 
                        dj = 1/250,
                        lowerPeriod = 2,
                        upperPeriod = 128,
                        make.pval = TRUE, n.sim = 10)
wt.image(my.w, n.levels = 250,
         legend.params = list(lab = "wavelet power levels"))

# Reconstruction, using only period 8

reconstruct(my.w, sel.period = 8, plot.waves = TRUE, lwd = c(1,2), legend.coords = "bottomleft")

reconstruct(my.w, sel.period = 16, plot.waves = TRUE, lwd = c(1,2), legend.coords = "bottomleft")

# actually uses the period closest to 10

my.w$Period[(my.w$Period > 9) & (my.w$Period < 11)]

# Compare two series with average powers calling:

x <- ts(calcium$ROI415, start = calcium$Time[1])
y <- ts(calcium$ROI416, start = calcium$Time[1])

my.data <- data.frame(x = x, y = y)
my.wx <- analyze.wavelet(my.data, "x", loess.span = 0,
                         dt = 1, dj = 1/20,
                         lowerPeriod = 2, upperPeriod = 64,
                         make.pval = TRUE,n.sim = 10)
my.wy <- analyze.wavelet(my.data, "y", loess.span = 0,
                         dt = 1, dj = 1/20,
                         lowerPeriod = 2, upperPeriod = 64,
                         make.pval = TRUE,n.sim = 10)

maximum.level = 1.001*max(my.wx$Power.avg, my.wy$Power.avg)
par(mfrow=c(1,2))
wt.avg(my.wx, maximum.level = maximum.level)
wt.avg(my.wy, maximum.level = maximum.level)

# White noise method

par(mfrow=c(1,1))

my.data <- data.frame(x = x)
my.w <- analyze.wavelet(my.data, "x",
                        method = "white.noise",
                        loess.span = 0,
                        dt = 1, 
                        dj = 1/250,
                        lowerPeriod = 2,
                        upperPeriod = 256,
                        make.pval = TRUE, n.sim = 10)
wt.image(my.w, color.key = "interval", n.levels = 250,
         legend.params = list(lab = "wavelet power levels"))

# Fourier randomization method

my.data <- data.frame(x = x)
my.w <- analyze.wavelet(my.data, "x",
                        method = "Fourier.rand",
                        loess.span = 0,
                        dt = 1, 
                        dj = 1/250,
                        lowerPeriod = 2,
                        upperPeriod = 256,
                        make.pval = TRUE, n.sim = 10)
wt.image(my.w, color.key = "interval", n.levels = 250,
         legend.params = list(lab = "wavelet power levels"))

# Plotting the power spectrum

par(mfrow=c(1,1))

wt.image(my.w, color.key = "quantile", n.levels = 250,
         legend.params = list(lab = "wavelet power levels", label.digits = 2))


# Grayscale

wt.image(my.w, n.levels = 250,
         legend.params = list(lab = "wavelet power levels", label.digits = 2),
         color.palette = "gray((n.levels):1/n.levels)",
         col.ridge = "blue")

# Time axis

my.data <- data.frame(x = x)

my.w <- analyze.wavelet(my.data, "x",
                        method = "Fourier.rand",
                        loess.span = 0,
                        dt = 1, 
                        dj = 1/250,
                        lowerPeriod = 2,
                        upperPeriod = 256,
                        make.pval = TRUE, n.sim = 10)

wt.image(my.w, color.key = "quantile", n.levels = 250,
         legend.params = list(lab = "wavelet power levels", label.digits = 2),
         timelab = "time")

# Analysis of bivariate time series

x <- ts(calcium$ROI415, start = calcium$Time[1])
y <- ts(calcium$ROI416, start = calcium$Time[1])

par(mfrow=c(1,1))

my.data <- data.frame(x = x, y = y)

my.wc <- analyze.coherency(my.data, my.pair = c("x", "y"),
                           loess.span = 0,
                           dt = 1, 
                           dj = 1/100,
                           lowerPeriod = 4, upperPeriod = 32,
                           make.pval = TRUE, n.sim = 10)


wc.image(my.wc, n.levels = 250, color.key = "interval",
         siglvl.contour = 0.1, siglvl.arrow = 0.05, which.arrow.sig = "wt",
         legend.params = list(lab = "cross-wavelet power levels"),
         timelab = "")

######


wc.image(my.wc, n.levels = 250,
         legend.params = list(lab = "cross-wavelet power levels"),
         timelab = "", periodlab = "period")


# A plot of the time-averaged cross-wavelet power

wc.avg(my.wc, siglvl = 0.01, sigcol = "red", sigpch = 20,
       periodlab = "period")

# Coherence

wc.image(my.wc, which.image = "wc", color.key = "interval", n.levels = 250,
         siglvl.contour = 0.1, siglvl.arrow = 0.05,
         legend.params = list(lab = "wavelet coherence levels"),
         timelab = "")


my.wc2 <- analyze.coherency(my.data, my.pair = c("x", "y"),
                           loess.span = 0, dt = 1, dj = 1/100,
                           lowerPeriod = 2, upperPeriod = 32,
                           window.type.t = 1, window.type.s = 1,
                           window.size.t = 5, window.size.s = 1,
                           make.pval = TRUE, n.sim = 10)

wc.image(my.wc2, which.image = "wc", color.key = "interval", n.levels = 250,
         siglvl.contour = 0.1, siglvl.arrow = 0.05,
         legend.params = list(lab = "wavelet coherence levels"),
         timelab = "")

my.wc3 <- analyze.coherency(my.data, my.pair = c("x", "y"),
                            loess.span = 0, dt = 1, dj = 1/100,
                            lowerPeriod = 2, upperPeriod = 64,
                            window.type.t = 3, window.type.s = 3,
                            window.size.t = 5, window.size.s = 1,
                            make.pval = TRUE, n.sim = 10)

wc.image(my.wc3, which.image = "wc", color.key = "interval", n.levels = 250,
         siglvl.contour = 0.1, siglvl.arrow = 0.05,
         legend.params = list(lab = "wavelet coherence levels"),
         timelab = "")

# Phase difference

wc.phasediff.image(my.wc, which.contour = "wc", use.sAngle = TRUE,
                   n.levels = 250, siglvl = 0.1,
                   legend.params = list(lab = "phase difference levels",
                                        lab.line = 3),
                   timelab = "")