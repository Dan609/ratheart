# The wavelet transform of x is computed as follows:
xcal <- calcium$ROI24

x <- ts(xcal, start = calcium$Time[1])

my.data <- data.frame(x = x)
my.w <- analyze.wavelet(my.data, "x",
                        loess.span = 0,
                        dt = 1, 
                        dj = 1/250,
                        lowerPeriod = 2,
                        upperPeriod = 128,
                        make.pval = TRUE, n.sim = 10)

#Plot the wavelet power spectrum, a series with a variable period

par(mfrow=c(1,1))

wt.image(my.w, n.levels = 250,
         legend.params = list(lab = "wavelet power levels"))

# White noise method

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

# The reconstructed series

my.rec <- reconstruct(my.w)
x.rec <- my.rec$series$x.r   #x: name of original series

# Reconstruction, using only period 16

reconstruct(my.w, sel.period = 16, plot.waves = TRUE, 
            lwd = c(1,2), legend.coords = "bottomleft")

# Reconstruction, using only period 8

reconstruct(my.w, sel.period = 8, plot.waves = TRUE, 
            lwd = c(1,2), legend.coords = "bottomleft")

# Compare two series with average powers calling:

x <- ts(calcium$ROI1, start = calcium$Time[1])
y <- ts(calcium$ROI9, start = calcium$Time[1])

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

# Analysis of bivariate time series
# Coherence

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

wc.image(my.wc, n.levels = 250,
         legend.params = list(lab = "cross-wavelet power levels"),
         timelab = "", periodlab = "period")


#

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

#

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



