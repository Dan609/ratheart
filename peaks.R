library(seewave)
library(phonTools)
library(signal)
library(oce)
library(rpanel)
library(ggplot2)

par(mfrow=c(1,1))

s.rate <- 14.29745
lenght 109.6 s

409
412
414
415
416
422

z <- ts(calcium$ROI415, start = calcium$Time[1])

# peaks detection

fpeaks(spec(z, f=s.rate), nmax = 5,
       threshold = 10,
       plot = T,
       title = TRUE,
       xlab = "Frequency (Hz)", ylab = "Amplitude",
       labels = TRUE, legend = TRUE, digits = 2)*1000


###
z <- ts(calcium$ROI333, start = calcium$Time[1])

specgram(z, 250, 14.29745, 10,
         overlap = ceiling(length(window)/2))

specgram(z, 500, 14.29745, 11,
         overlap = ceiling(length(window)/2))

specgram(z, 500, 14.29745, 12, 10)

specgram(z, 500, 14.29745, 11)

specgram(z, 500, 14.29745, 11)


?specgram
col = heat.colors(256)
###

spectro3D(z, 14.29745, wl = 512, wn = "hanning", zp = 0,
          ovlp = 0, norm = TRUE, correction = "none", fftw = FALSE,
          plot = TRUE,
          magt = 100, magf = 100, maga = 20,
          palette = reverse.terrain.colors)

?spectro3D
###

spectrogram(z, 14.29745, windowlength = 5000, 0.1,
            timestep = 1000,
            maxfreq = 10,
            colors = TRUE)

spectrogram(z, 14.29745)

?spectrogram
###

spectro(z, 14.29745, wl = 512, wn = "hanning", zp = 0,
        ovlp = 0, complex = FALSE, norm = TRUE, correction="none",
        fftw = FALSE, plot = TRUE,
        flog = FALSE, grid = TRUE, osc = FALSE, scale = TRUE, cont = FALSE,
        collevels = NULL, palette = spectro.colors,
        contlevels = NULL, colcont = "black",
        colbg = "white", colgrid = "black",
        colaxis = "black", collab="black",
        cexlab = 1, cexaxis = 1, 
        tlab = "Time (s)",
        flab = "Frequency (kHz)",
        alab = "Amplitude",
        scalelab = "Amplitude",
        main = NULL,
        scalefontlab = 1, scalecexlab =0.75,
        axisX = TRUE, axisY = TRUE, tlim = NULL, trel = TRUE,
        flim = NULL, flimd = NULL,
        widths = c(6,1), heights = c(3,1),
        oma = rep(0,4),
        listen=FALSE)

?spectro
###

spectro(z, f=14.29745, osc=TRUE, flim=c(0,0.0001))

spectro(z, f=14.29745, osc=TRUE)

spectro(z,f=22050,ovlp=85,zp=16,osc=TRUE,
        cont=TRUE,contlevels=seq(-30,0,20),colcont="red",
        lwd=1.5,lty=2,palette=reverse.terrain.colors)

### не похоже на првильный результат:

which(Mod(fft(z)) == max(abs(Mod(fft(z)))))*s.rate/length(z)

##

ggspectro(z, 14.29745)

###

z <- ts(calcium$ROI416, start = calcium$Time[1])
plot(z)
plot(fft(z), ylim = c(-2000, 2000))







