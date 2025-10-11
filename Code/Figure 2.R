# remove everything
rm(list=ls(all.names=TRUE))

# libraries
set.seed(1)
# source("Z:/Project Tutorial dRMST vs RD/Code/RD.dRMST.R")
library(dplyr)
library(tidyverse)
library(survival)
library(survminer)
library(patchwork)
library(data.table)
# remotes::install_github("CHMMaas/PredictionTools")
library(PredictionTools)

# sample size
n <- 1000000

###
### create risk
###
# Coefficient of risk that controls the C-index
beta <- 0.5

# draw risk for each individual
x <- rnorm(n, mean=0, sd=1)

# calculate the linear predictor
lp <- beta*x

###
### create treatment
###
# assign one half to control and other to treatment
z <- c(rep(0, n/2), rep(1, n/2))

###
### create hazard by combining risk and treatment
###
# Baseline hazard that controls the event rate
bh <- 0.15

# treatment effect
TE <- log(0.8)

# combine risk and treatment into one hazard
h <- bh * exp(lp + z*TE)

###
### sample event times and make survival curve
###
# sample event times
# time <- rexp(n=n, rate=h)
s <- 2
time <- rweibull(n=n, shape=s,
                 scale=1/(gamma(1+1/s)*exp(TE)*h)) # increasing failure rate

# survival curve
S <- Surv(time=time, event=rep(1, n))
cat("C-index:", round(concordance(S ~ I(-lp), subset=z==0)$concordance, 2), "\n")

# limit S to horizon
max.horizon <- 25
S.max <- S
S.max[S[, 1] > max.horizon, 2] <- 0
S.max[S[, 1] > max.horizon, 1] <- max.horizon

###
### make risk groups
###
# sort on lp and make groups on lp
groups <- as.numeric(cut(lp, breaks=quantile(lp,
                                             probs=seq(0, 1, by=0.25),
                                             include.lowest=TRUE)))

###
### make Kaplan-Meier plot, stratify by treatment and risk group
###
horizons <- c(3, 5, 10, 15)
KM <- survival::survfit(S.max ~ z + groups, data=data.frame(S.max, z, groups))
space <- 0.3
plot.KM <- survminer::ggsurvplot(KM,
                                 times=seq(0, max.horizon, by=1),
                                 data=data.frame(time, z, x),
                                 palette=c("#04bca2", "#04b4f4", "#9494fc", "#e46cf4",
                                           "#04bca2", "#04b4f4", "#9494fc", "#e46cf4"),
                                 xlim=c(0, max.horizon),
                                 break.x.by = 5,
                                 legend="none")$plot

# shade area in between treatment and risk stratified Kaplan-Meier curves
# extract data from layer_data
dat <- plot.KM |>
  layer_data() |>
  arrange(colour, group) |>
  mutate(trt = consecutive_id(group), .by = colour) |>
  select(stage = colour, trt, time = x, surv = y)

# duplicate the observations where the step occurs
dat_step <- dat |>
  arrange(stage, trt, time) |>
  mutate(surv_lag = lag(surv, default = 1)) |>
  filter(surv_lag > surv, .by = c(stage, trt)) |>
  mutate(surv = surv_lag, .keep = "unused")

# arranging the data in the correct order, i.e.,
# arrange one treatment group by time and the second in reverse order
dat_merge <- bind_rows(dat, dat_step) |>
  split(~stage) |>
  lapply(\(x) {
    bind_rows(
      x |> filter(trt == 1) |> arrange(time, desc(surv)),
      x |> filter(trt == 2) |> arrange(desc(time), surv)
    )
  }) |>
  bind_rows()

# add polygon to KM plot
plot.KM <- plot.KM +
  geom_polygon(data = dat_merge, aes(time, surv, fill = I(stage)), alpha = .2)

# add horizons
for (i in 1:length(horizons)){
  plot.KM <- plot.KM +
    ggplot2::geom_vline(xintercept=horizons[i], linetype="dashed") +
    ggplot2::annotate("text", x=horizons[i]-space, y=ifelse(i==1, space, 1), label=LETTERS[i])
}

###
### calculate RD and dRMST in each risk group at each time horizon
###
RD <- c()
dRMST <- c()
event.rates <- c()
C.indexes <- c()
for (horizon.sel in horizons){
  for (risk.group.i in 1:4){
    # calculate RD and dRMST in risk strata
    out <- PredictionTools::calculate.RD.dRMST(S=S.max[groups==risk.group.i, ],
                               W=z[groups==risk.group.i],
                               horizon=horizon.sel)
    RD <- c(RD, out$RD)
    dRMST <- c(dRMST, out$dRMST)
  }

  event.rates <- c(event.rates,
                   sum(time[z==0]<=horizon.sel)/length(time[z==0])*100)
  S.h.sel <- S
  S.h.sel[S[, 1] > horizon.sel, 2] <- 0
  S.h.sel[S[, 1] > horizon.sel, 1] <- horizon.sel
  C.indexes <- c(C.indexes,
                 concordance(S.h.sel ~ I(-lp), subset=z==0)$concordance)
}
df.RD.dRMST <- data.frame(horizon=rep(horizons, each=4),
                           risk.group=rep(1:4, length(horizons)),
                           RD, dRMST)

###
### plot RD and dRMST
###
s <- max(df.RD.dRMST$dRMST)/max(df.RD.dRMST$RD * 100)
titles.lab <- paste0(LETTERS[1:length(horizons)], ". ", horizons,
                     "-year horizon, event rate: ",
                     sprintf("%.0f", event.rates), "%")
plot.RD.dRMST <- ggplot2::ggplot(data=df.RD.dRMST,
                                  ggplot2::aes(x=as.numeric(risk.group)))+
  ggplot2::geom_line(aes(y=RD * 100), col="#AD002AFF", alpha=0.5) +
  ggplot2::geom_point(aes(y=RD * 100), size=2, shape=18, col="#AD002AFF") +
  ggplot2::geom_line(aes(y=dRMST / s), col="#42B540FF", alpha=0.5) +
  ggplot2::geom_point(aes(y=dRMST / s), size=2, shape=18, col="#42B540FF") +
  ggplot2::facet_wrap(~ horizon, ncol=length(horizons), drop=FALSE,
                      labeller=as_labeller(c("3"=titles.lab[1],
                                             "5"=titles.lab[2],
                                             "10"=titles.lab[3],
                                             "15"=titles.lab[4]))) +
  ggplot2::scale_y_continuous(
    name="RD in %",
    sec.axis=ggplot2::sec_axis(~ . * s,
                               name="\u0394RMST in years")) +
  ggplot2::scale_x_discrete(limits=as.factor(1:4)) +
  ggplot2::xlab("Risk group") +
  ggplot2::theme_bw() +
  ggplot2::theme(axis.title.y=ggplot2::element_text(color="#AD002AFF"),
                 axis.title.y.right=ggplot2::element_text(color="#42B540FF"),
                 panel.grid.minor=ggplot2::element_blank(),
                 panel.grid.major=ggplot2::element_blank(),
                 plot.margin=unit(c(0, 0, 0, 0.8), "cm"))

###
### combine KM and RD, dRMST plot
###
ggsave(file="./Figure 2.png",
       ggpubr::ggarrange(plot.KM, plot.RD.dRMST,
                         nrow=2, ncol=1, heights=c(2, 1)),
       width=12, height=10, dpi=300)

