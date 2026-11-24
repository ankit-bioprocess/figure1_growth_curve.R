setwd("/Users/ankitthakor/Downloads/untitled folder 5")

files <- list.files(pattern = "pio07.*\\.csv$", ignore.case = TRUE)
if (length(files) == 0) stop("No Pio07 CSV found in this folder. Copy the CSV here.")
d <- read.csv(files[1])

x <- as.numeric(d$hours_since_experiment_created)
y <- as.numeric(d$od_reading)
ok <- !is.na(x) & !is.na(y)
x <- x[ok]
y <- y[ok]

# 10-minute medians
bin <- floor(x * 6)
med <- tapply(y, bin, median)
t10 <- as.numeric(names(med)) / 6 + 1 / 12
od10 <- as.numeric(med)
ln10 <- log(od10)

# Steepest 2-hour window gives mu max
n <- 12
nwin <- length(t10) - n + 1
slopes <- rep(NA, nwin)
centres <- rep(NA, nwin)
icpts <- rep(NA, nwin)
for (i in 1:nwin) {
  idx <- i:(i + n - 1)
  cf <- coef(lm(ln10[idx] ~ t10[idx]))
  icpts[i] <- cf[1]
  slopes[i] <- cf[2]
  centres[i] <- mean(t10[idx])
}
k <- which.max(slopes)
mu_max <- slopes[k]
lag <- (min(ln10) - icpts[k]) / mu_max
cat("mu max =", round(mu_max, 4), " lag =", round(lag, 2), "h\n")

# Time growth slowed to 10% of mu max
after <- which(slopes < 0.1 * mu_max & seq_along(slopes) > k)
t_slow <- centres[after[1]]

# Gompertz fit from the OD minimum onwards
t_min <- t10[which.min(od10)]
sel <- t10 >= t_min
gx <- t10[sel]
gy <- ln10[sel] - min(ln10[sel])
gfit <- nls(gy ~ A * exp(-exp(mu * exp(1) / A * (lam - gx) + 1)),
            start = list(A = max(gy), mu = 0.25, lam = 6.5))

# Stage bands
s0 <- c(0.6, 3.6, 14.2, 16.5, 31.7)
s1 <- c(3.6, 14.2, 16.5, 31.7, 47.8)
cols <- c("#e6e6e6", "#dbe9f6", "#cfe8cf", "#f9e6c9", "#f0dcdc")

png("figure1.png", width = 7.2, height = 2.9, units = "in", res = 200)
par(mfrow = c(1, 2), mar = c(4, 4, 2, 1), cex = 0.75, las = 1)

# Panel (a): raw OD
plot(NA, xlim = c(0, 48), ylim = c(0, 0.40), xaxs = "i",
     xlab = "Time since experiment created (h)", ylab = "OD (90-degree scatter, a.u.)",
     main = "(a) Raw OD")
rect(s0, -1, s1, 1, col = cols, border = NA)
points(x, y, pch = ".", col = "#9bbbbb")
lines(t10, od10, col = "#1f4e79", lwd = 1.4)
abline(v = 26.7, lty = 3)
text(26.9, 0.03, "26.7 h", cex = 0.8, adj = 0)
legend("topleft", legend = c("Raw (5 s)", "10-min median"), col = c("#9bbbbb", "#1f4e79"),
       lty = c(NA, 1), pch = c(16, NA), bty = "n", cex = 0.8)

# Panel (b): log OD with fits
plot(NA, xlim = c(0, 48), ylim = range(od10), log = "y", xaxs = "i",
     xlab = "Time since experiment created (h)", ylab = "OD (log, a.u.)",
     main = "(b) Log OD with fits")
rect(s0, 1e-6, s1, 10, col = cols, border = NA)
lines(t10, od10, col = "#1f4e79", lwd = 1.4)
tg <- seq(lag - 1.5, 20, length.out = 50)
lines(tg, exp(icpts[k] + mu_max * tg), lty = 2, col = "#c0392b", lwd = 1.2)
tg2 <- seq(t_min, 47.8, length.out = 200)
lines(tg2, exp(predict(gfit, newdata = data.frame(gx = tg2)) + min(ln10[sel])), lty = 4, col = "#2e7d32", lwd = 1.2)
points(lag, exp(min(ln10)), pch = 16, col = "#c0392b")
abline(v = t_slow, lty = 3, col = "grey40")
legend("bottomright", legend = c("10-min median", expression(paste("Tangent at ", mu[max])),
                                 "Gompertz fit", paste("Lag =", round(lag, 1), "h")),
       col = c("#1f4e79", "#c0392b", "#2e7d32", "#c0392b"), lty = c(1, 2, 4, NA),
       pch = c(NA, NA, NA, 16), bty = "n", cex = 0.75)

dev.off()
cat("Done: figure1.png saved\n")
getwd()
file.exists("figure1.png")
normalizePath("figure1.png")
