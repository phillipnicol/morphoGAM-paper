setwd(here::here("analysis/ratio_diagnostic/code"))

set.seed(1)
n <- 1000

t <- seq(0, 2*pi, length.out = n)

figure8 <- data.frame(
  x = sin(t),
  y = sin(2 * t)
)

#Add gaussian noise
set.seed(123)
figure8$x <- figure8$x + rnorm(n, sd = 0.1)
figure8$y <- figure8$y + rnorm(n, sd = 0.1)

plot(
  figure8$x, figure8$y,
  pch = 16,
  cex = 0.4,
  asp = 1
)

library(MorphoGAM)
fit <- CurveFinderInteractive(figure8,loop=TRUE)

#Side by side curve and coordinate
library(ggpubr)

p <- ggarrange(fit$curve.plot, fit$coordinate.plot, nrow=1,ncol=2)

ggsave(p, filename="../plots/figure8.png", width=10, height=5)
