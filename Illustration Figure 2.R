# remove everything
rm(list=ls(all.names=TRUE))

# libraries
set.seed(1)
source("Z:/Project Tutorial dRMST vs ARD/Code/ARD.dRMST.R")
library(dplyr)
library(survival)
library(survminer)
library(patchwork)

# sample size
n <- 100000 # TODO: 1,000,000

# draw risk for each individual
x <- rnorm(n, mean=0, sd=1)

# vary the hazard ratio for treatment, the lower the HR the larger the treatment effect
# vary the coefficient of risk to vary C-index, the higher the beta the higher the C-index
# vary the baseline hazard to vary the outcome rate, the higher bh the higher the outcome rate
sens.df <- list(name.bh=rep(c("Low.OR", "Medium.OR", "High.OR"), each=3),
                bh=c(0.042, 0.037, 0.027,
                     0.135, 0.132, 0.122,
                     0.245, 0.265, 0.335),  # vary event rate, 20%, 50%, 70%
                name.C=rep(c("Low.C", "Medium.C", "High.C"), 3),
                beta=c(0.352, 0.750, 1.300,
                       0.355, 0.780, 1.480,
                       0.348, 0.798, 1.550), # vary C-index, 0.6, 0.7, 0.8
                name.TE=c("Small.TE", "Large.TE"),
                TE=c(log(0.8), log(0.5)))    # vary HR of treatment, 0.8, 0.5

# create new results
new.results <- TRUE
# set time horizon
horizon <- 5
# set up empty data frame for ARD and dRMST
df.ARD.dRMST <- c()
if (new.results){
  # calculate ARD and dRMST using different settings
  settings.df <- data.frame(c("Event rate for control",
                              "KM for control",
                              "C-index",
                              "HR for treatment",
                              paste("Risk control group", 1:4)))
  for (TE in sens.df$TE){
    for (nr.plot in 1:9){
      bh <- sens.df$bh[nr.plot]
      beta <- sens.df$beta[nr.plot]
      cat(sens.df$name.bh[sens.df$bh==bh], " and ",
          sens.df$name.C[sens.df$beta==beta], " and ",
          sens.df$name.TE[sens.df$TE==TE], "\n")

      # calculate the linear predictor
      lp <- beta*x

      # assign one half to control and other to treatment
      z <- c(rep(0, n/2), rep(1, n/2))

      # combine risk and treatment into one hazard
      h <- bh * exp(lp + z*TE)

      # sample event times
      time <- rexp(n=n, rate=h)

      # non-proportional
      # time.0 <- rweibull(n=n, shape=1, scale=1/h)
      # s <- 2
      # time.1 <- rweibull(n=n, shape=s, scale=1/(gamma(1+1/s)*exp(TE)*h))
      # time <- (1-z)*time.0 + z*time.1

      # survival curve
      S <- Surv(time=time, event=rep(1, n))

      # limit S to horizon
      S.h <- S
      S.h[S[, 1] > horizon, 2] <- 0
      S.h[S[, 1] > horizon, 1] <- horizon

      # Cox model
      cph <- coxph(S.h ~ x + z)
      overall.HR <- exp(cph$coefficients["z"])

      # overall ARD and dRMST
      overall <- calculate.ARD.dRMST(S=S.h, W=z, horizon=horizon)

      # sort on lp and make groups on lp
      groups <- as.numeric(cut(lp,
                               breaks=quantile(lp, probs=seq(0, 1, by=0.25),
                                               include.lowest=TRUE)))

      # save overall characteristics
      event.rate <- sum(time[z==0]<=horizon)/length(time[z==0])*100
      C.index <- concordance(S.h ~ I(-lp), subset=z==0)$concordance
      cat("Plot", nr.plot, "\n",
          "Event rate:", round(event.rate, 0), "\n",
          "C-index:", round(C.index, 2), "\n")

      # make KM plot for each risk group
      plot.KM <- list()
      plot.hazard <- list()
      HR <- c()
      ARD <- c()
      dRMST <- c()
      risk.control <- c()
      for (risk.group.i in 1:4){
        # calculate HR in risk strata
        HR <- c(HR, exp(coxph(S.h ~ z,
                              subset=groups==risk.group.i)$coefficients))

        # calculate ARD and dRMST in risk strata
        out <- calculate.ARD.dRMST(S=S.h[groups==risk.group.i, ],
                                   W=z[groups==risk.group.i],
                                   horizon=horizon)
        ARD <- c(ARD, out$ARD)
        dRMST <- c(dRMST, out$dRMST)

        # make KM plot
        KM.i <- survival::survfit(S.h ~ z, subset = groups==risk.group.i)
        plot.KM.i <- survminer::ggsurvplot(KM.i,
                                           data=data.frame(time, z, x),
                                           palette=c("#AD002AFF", "#42B540FF"),
                                           xlim=c(0, horizon),
                                           break.x.by = 2,
                                           legend="none")
        risk.control <- c(risk.control,
                          round((1-summary(KM.i, times=horizon, extend=TRUE)$surv[1])*100, 1))

        # edit theme of plot
        plot.KM.i$plot <- plot.KM.i$plot +
          theme(axis.text.x=element_text(size=12),
                axis.title.x=element_text(size=12),
                axis.line.x=element_line(colour="black"),
                axis.text.y=element_text(size=12),
                axis.title.y=element_text(size=12),
                axis.line.y=element_line(colour="black"),
                panel.grid.major=element_blank(),
                panel.grid.minor=element_blank(),
                panel.background=element_rect(fill='transparent'),
                plot.background=element_rect(fill='transparent', color=NA))

        # remove y-axis if not risk group 1
        if (risk.group.i!=1){
          plot.KM.i$plot <- plot.KM.i$plot +
            theme(axis.text.y=element_blank(),
                  axis.ticks.y=element_blank(),
                  axis.title.y=element_blank(),
                  axis.line.y=element_blank())
        }
        plot.KM[[risk.group.i]] <- plot.KM.i$plot

        # plot hazard over time
        hazard.df <- data.frame(cumhaz=summary(KM.i, extend=TRUE)$cumhaz,
                                time=summary(KM.i, extend=TRUE)$time,
                                z=as.numeric(summary(KM.i, extend=TRUE)$strata)-1)
        if (nrow(hazard.df) != 0){
          hazard.df$hazard <- c(hazard.df$cumhaz[1], diff(hazard.df$cumhaz))
          plot.hazard[[risk.group.i]] <- ggplot(data=hazard.df,
                                                aes(y=cumhaz, x=time,
                                                    color=as.factor(z))) +
            ggplot2::geom_line(linewidth=1) +
            ggplot2::theme_light(base_size=15) +
            ggplot2::ylim(0, ifelse(bh==sens.df$bh[1], 0.75, 7)) +
            ggplot2::xlim(0, horizon) +
            ggplot2::ylab("Cumulative hazard") +
            ggplot2::xlab("Time") +
            ggplot2::theme(legend.position="none")
        } else{
          print("No cumulative hazard can be computed")
          plot.hazard[[risk.group.i]] <- ggplot() + theme_void()
        }
      }
      # save results
      results <- c(event.rate,
                   round(summary(survival::survfit(S.h ~ z,
                                                   data=data.frame(S.h, z)),
                                 times=horizon,
                                 extend=TRUE)$surv[1]*100, 1),
                   C.index,
                   round(as.numeric(overall.HR), 2),
                   risk.control)
      settings.df <- cbind(settings.df, results)

      # KM plots
      KM.plots <- ggpubr::ggarrange(plot.KM[[1]], plot.KM[[2]],
                                    plot.KM[[3]], plot.KM[[4]],
                                    ncol=4, nrow=1, align="h")

      # HR plot
      HR.plot <- ggplot(data=data.frame(risk=1:4, HR=HR), aes(x=risk, y=HR)) +
        ggplot2::geom_hline(yintercept=0)+
        ggplot2::geom_hline(yintercept=overall.HR, col="darkgrey",
                            linetype="dashed", linewidth=1.2)+
        ggplot2::geom_point(size=5, shape=18, col="#AD002AFF") +
        ggplot2::ylab(paste0(horizon, "-year HR"))+
        ggplot2::ylim(c(min(0, HR), max(HR, 1)))+
        ggplot2::theme(axis.text.x=ggplot2::element_blank(),
                       axis.title.x=ggplot2::element_blank(),
                       axis.ticks.x=ggplot2::element_blank(),
                       axis.text.y=ggplot2::element_text(size=12),
                       axis.title.y=ggplot2::element_text(size=12),
                       axis.line.y=ggplot2::element_line(colour="black"),
                       panel.grid.major=ggplot2::element_blank(),
                       panel.grid.minor=ggplot2::element_blank(),
                       panel.background=ggplot2::element_rect(fill='transparent'),
                       plot.background=ggplot2::element_rect(fill='transparent', color=NA))

      # ARD plot
      # cat("ARD:", sort(round(ARD*100, 1)), "\n")
      ARD.plot <- ggplot(data=data.frame(risk=1:4, ARD=ARD), aes(x=risk, y=ARD)) +
        ggplot2::geom_hline(yintercept=0)+
        ggplot2::geom_hline(yintercept=overall$ARD, col="darkgrey",
                            linetype="dashed", linewidth=1.2)+
        ggplot2::geom_point(size=5, shape=18, col="#AD002AFF") +
        ggplot2::ylab(paste0(horizon, "-year ARD"))+
        ggplot2::ylim(c(min(0, ARD), max(ARD)))+
        ggplot2::theme(axis.text.x=ggplot2::element_blank(),
                       axis.title.x=ggplot2::element_blank(),
                       axis.ticks.x=ggplot2::element_blank(),
                       axis.text.y=ggplot2::element_text(size=12),
                       axis.title.y=ggplot2::element_text(size=12),
                       axis.line.y=ggplot2::element_line(colour="black"),
                       panel.grid.major=ggplot2::element_blank(),
                       panel.grid.minor=ggplot2::element_blank(),
                       panel.background=ggplot2::element_rect(fill='transparent'),
                       plot.background=ggplot2::element_rect(fill='transparent', color=NA))

      # dRMST plot
      # cat("dRMST:", sort(round(dRMST, 2)), "\n")
      dRMST.plot <- ggplot(data=data.frame(risk=1:4, dRMST=dRMST), aes(x=risk, y=dRMST)) +
        ggplot2::geom_hline(yintercept=overall$dRMST, col="darkgrey",
                            linetype="dashed", linewidth=1.2)+
        ggplot2::geom_point(size=5, shape=18, col="#AD002AFF") +
        ggplot2::ylab(paste0(horizon, "-year \u0394RMST"))+
        ggplot2::scale_y_continuous(limits=c(min(0, dRMST), max(dRMST)),
                                    labels=scales::number_format(accuracy=0.01))+
        ggplot2::xlab("Risk group")+
        ggplot2::theme(axis.text.x=ggplot2::element_text(size=12),
                       axis.title.x=ggplot2::element_text(size=12),
                       axis.line.x=ggplot2::element_line(colour="black"),
                       axis.text.y=ggplot2::element_text(size=12),
                       axis.title.y=ggplot2::element_text(size=12),
                       axis.line.y=ggplot2::element_line(colour="black"),
                       panel.grid.major=ggplot2::element_blank(),
                       panel.grid.minor=ggplot2::element_blank(),
                       panel.background=ggplot2::element_rect(fill='transparent'),
                       plot.background=ggplot2::element_rect(fill='transparent', color=NA))

      # save ARD and dRMST results
      norm.ARD <- (ARD-mean(ARD))/sd(ARD) #TODO: scale
      norm.dRMST <- (dRMST-mean(dRMST))/sd(dRMST)
      df.ARD.dRMST <- rbind(df.ARD.dRMST,
                            cbind(rep(sens.df$name.TE[sens.df$TE==TE], 4),
                                  rep(event.rate, 4),
                                  rep(C.index, 4),
                                  1:4,
                                  ARD, norm.ARD,
                                  dRMST, norm.dRMST))

      # save plot
      ggsave(file=paste0("Z:/Project Tutorial dRMST vs ARD/Illustration risk-stratified/Panel Figures/Fig",
                         ifelse(sens.df$name.TE[sens.df$TE==TE]=="Small.TE", "1", "2"),
                         letters[nr.plot], ".",
                         sens.df$name.bh[sens.df$bh==bh], ".",
                         sens.df$name.C[sens.df$beta==beta], ".",
                         sens.df$name.TE[sens.df$TE==TE], ".png"),
             plot=ggpubr::ggarrange(KM.plots, ARD.plot, dRMST.plot,
                                    nrow=3, ncol=1, align="h"),
             width=8, height=8, dpi=300)

      # hazard plots
      ggsave(file=paste0("Z:/Project Tutorial dRMST vs ARD/Illustration risk-stratified/Hazard plots/hazard.plot.",
                         sens.df$name.bh[sens.df$bh==bh], ".",
                         sens.df$name.C[sens.df$beta==beta], ".",
                         sens.df$name.TE[sens.df$TE==TE], ".png"),
             plot=ggpubr::ggarrange(plot.hazard[[1]], plot.hazard[[2]],
                                    plot.hazard[[3]], plot.hazard[[4]],
                                    ncol=2, nrow=2, align="h"),
             width=10, height=10, dpi=300)
    }
  }
  settings.df <- data.frame(settings.df)
  colnames(settings.df) <- c("Outcome", paste("Setting", 1:(ncol(settings.df)-1)))
  openxlsx::write.xlsx(data.frame(settings.df),
                       colNames=FALSE,
                       file=paste0("Z:/Project Tutorial dRMST vs ARD/Illustration risk-stratified/settings.xlsx"))
  print(settings.df)
  df.ARD.dRMST <- as.data.frame(df.ARD.dRMST)
  colnames(df.ARD.dRMST) <- c("HR", "OR", "Cindex", "risk.group",
                              "ARD", "norm.ARD", "dRMST", "norm.dRMST")
  save(df.ARD.dRMST,
       file="Z:/Project Tutorial dRMST vs ARD/Illustration risk-stratified/df.ARD.dRMST.Rdata")
} else{
  load("Z:/Project Tutorial dRMST vs ARD/Illustration risk-stratified/df.ARD.dRMST.Rdata")
}

# create plots
for (TE in sens.df$name.TE){
  # initialize data frame
  df.ARD.dRMST.select <- df.ARD.dRMST |>
    filter(HR==TE) |>
    select(-HR) |>
    mutate(setting=rep(1:9, each=4))

  # plot titles
  event.rate <- unique(df.ARD.dRMST.select[, "OR"])
  Cindex <- unique(df.ARD.dRMST.select[, "Cindex"])
  titles <- paste0(c("A. ", "B. ", "C. ",
                     "D. ", "E. ", "F. ",
                     "G. ", "H. ", "I. "),
                   "Event rate: ", sprintf("%.0f", as.numeric(event.rate)),
                   "%, C-index: ", sprintf("%.2f", as.numeric(Cindex)))
  title.lab <- c("1"=titles[1],
                 "2"=titles[2],
                 "3"=titles[3],
                 "4"=titles[4],
                 "5"=titles[5],
                 "6"=titles[6],
                 "7"=titles[7],
                 "8"=titles[8],
                 "9"=titles[9])

  # horizontal lines
  ordered.dRMST <- sort(c(df.ARD.dRMST.select |>
                            filter(setting==5) |>
                            select(dRMST) |>
                            mutate_all(as.numeric))$dRMST)
  scarce.treat <- ordered.dRMST[3]+diff(ordered.dRMST[3:4])/2
  costly.treat <- ordered.dRMST[2]+diff(ordered.dRMST[2:3])/2
  inexpensive.treat <- ordered.dRMST[1]+diff(ordered.dRMST[1:2])/2

  # make plot
  df <- data.frame(sapply(df.ARD.dRMST.select, as.numeric))
  s <- max(df$dRMST)/max(df$ARD * 100)
  if (TE=="Small.TE"){
    red.panels <- c(5, 6, 7, 8, 9)
  } else{
    red.panels <- c(6, 8, 9)
  }
  plot <- ggplot2::ggplot(data=df,
                          ggplot2::aes(x=as.numeric(risk.group)))+
    ggplot2::geom_line(aes(y=ARD * 100), col="#AD002AFF", alpha=0.5) +
    ggplot2::geom_point(aes(y=ARD * 100), size=2, shape=18, col="#AD002AFF") +
    ggplot2::geom_line(aes(y=dRMST / s), col="#42B540FF", alpha=0.5) +
    ggplot2::geom_point(aes(y=dRMST / s), size=2, shape=18, col="#42B540FF") +
    ggplot2::scale_y_continuous(
      name=paste0(horizon, "-year ARD in %"),
      sec.axis=ggplot2::sec_axis(~ . * s,
                                 name=paste0(horizon, "-year \u0394RMST in years"))) +
    geom_rect(data=subset(df, setting %in% red.panels),
              fill=NA, colour="#AD002AFF", linewidth=1,
              xmin=-Inf, xmax=Inf, ymin=-Inf, ymax=Inf) +
    ggplot2::facet_wrap(~ setting, ncol=3, drop=FALSE,
                        labeller=as_labeller(title.lab)) +
    ggplot2::scale_x_discrete(limits=as.factor(1:4)) +
    ggplot2::xlab("Risk group") +
    ggplot2::theme_bw() +
    ggplot2::theme(axis.title.y=ggplot2::element_text(color="#AD002AFF"),
                   axis.title.y.right=ggplot2::element_text(color="#42B540FF"),
                   panel.grid.minor=ggplot2::element_blank(),
                   panel.grid.major=ggplot2::element_blank())
  show(plot)
  ggsave(file=paste0("Z:/Project Tutorial dRMST vs ARD/Illustration risk-stratified/",
                     TE, ".png"),
         plot=plot,
         width=8, height=6, dpi=300)
}
