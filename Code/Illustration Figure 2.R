# remove everything
rm(list=ls(all.names=TRUE))

# libraries
set.seed(1)
source("Z:/Project Tutorial dRMST vs ARD/Code/ARD.dRMST.R")
library(dplyr)
library(survival)
library(survminer)
library(patchwork)

# set file path
file.path <- "Z:/Project Tutorial dRMST vs ARD/Illustration Figure 2/"

# sample size
n <- 1000000 # TODO: 1,000,000

# draw risk for each individual
x <- rnorm(n, mean=0, sd=1)

# vary the coefficient of risk to vary C-index, the higher the beta the higher the C-index
# vary the baseline hazard to vary the outcome rate, the higher bh the higher the outcome rate
sens.df.all <- list(const.FR.Small.TE=list(name.bh=rep(c("Low.OR", "Medium.OR", "High.OR"), each=3),
                                           bh=c(0.010, 0.007, 0.002,
                                                0.034, 0.027, 0.012,
                                                0.068, 0.063, 0.059),
                                           name.C=rep(c("Low.C", "Medium.C", "High.C"), 3),
                                           beta=c(0.352, 1.00, 2.2,
                                                  0.355, 1.05, 2.7,
                                                  0.348, 1.10, 3.1)),
                    const.FR.Large.TE=list(name.bh=rep(c("Low.OR", "Medium.OR", "High.OR"), each=3),
                                           bh=c(0.010, 0.0070, 0.002,
                                                0.034, 0.0268, 0.012,
                                                0.069, 0.063, 0.060),
                                           name.C=rep(c("Low.C", "Medium.C", "High.C"), 3),
                                           beta=c(0.352, 1.00, 2.2,
                                                  0.360, 1.05, 2.7,
                                                  0.348, 1.10, 3.1),
                                           name.TE=c("Large.TE")),
                    incr.FR.Small.TE=list(name.bh=rep(c("Low.OR", "Medium.OR", "High.OR"), each=3),
                                          bh=c(0.0441, 0.038, 0.021,
                                               0.082, 0.074, 0.049,
                                               0.116, 0.112, 0.106),
                                          name.C=rep(c("Low.C", "Medium.C", "High.C"), 3),
                                          beta=c(0.170, 0.487, 1.08,
                                                 0.185, 0.510, 1.4,
                                                 0.173, 0.54, 1.62)),
                    incr.FR.Large.TE=list(name.bh=rep(c("Low.OR", "Medium.OR", "High.OR"), each=3),
                                          bh=c(0.072, 0.06, 0.033,
                                               0.134, 0.118, 0.075,
                                               0.187, 0.18, 0.172),
                                          name.C=rep(c("Low.C", "Medium.C", "High.C"), 3),
                                          beta=c(0.170, 0.487, 1.1,
                                                 0.178, 0.512, 1.4,
                                                 0.177, 0.532, 1.5)))

# create new results
new.results <- TRUE
# set time horizon
horizon <- 10
# set up empty data frame for ARD and dRMST
df.ARD.dRMST <- c()
if (new.results){
  # calculate ARD and dRMST using different settings
  settings.df <- data.frame(c("Event rate for control",
                              "KM for control",
                              "C-index",
                              "HR for treatment",
                              paste("Risk control group", 1:4)))
  for (name.FR in c("const.FR", "incr.FR")){
    for (name.TE in c("Small.TE", "Large.TE")){
      sens.df <- eval(parse(text=paste0("sens.df.all$", name.FR, ".", name.TE)))
      TE <- ifelse(name.TE=="Small.TE", log(0.8), log(0.5))

      for (nr.plot in 1:9){
        bh <- sens.df$bh[nr.plot]
        beta <- sens.df$beta[nr.plot]
        cat(name.FR, " and ",
            sens.df$name.bh[sens.df$bh==bh], " and ",
            sens.df$name.C[sens.df$beta==beta], " and ",
            name.TE, "\n")

        # calculate the linear predictor
        lp <- beta*x

        # assign one half to control and other to treatment
        z <- c(rep(0, n/2), rep(1, n/2))

        # combine risk and treatment into one hazard
        h <- bh * exp(lp + z*TE)

        # sample event times
        if (name.FR=="const.FR"){
          # constant failure rate
          time <- rexp(n=n, rate=h)
        } else if (name.FR=="incr.FR"){
          # increasing failure rate
          shape <- 2
          time <- rweibull(n=n, shape=shape, scale=1/(gamma(1+1/shape)*exp(TE)*h))
        }

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
                                 breaks=quantile(lp,
                                                 probs=seq(0, 1, by=0.25),
                                                 include.lowest=TRUE)))

        # save overall characteristics
        event.rate <- sum(time[z==0]<=horizon)/length(time[z==0])*100
        C.index <- concordance(S.h ~ I(-lp), subset=z==0)$concordance
        cat("Plot", nr.plot, "\n",
            "Event rate:", round(event.rate, 0), "\n",
            "C-index:", round(C.index, 2), "\n")

        # make KM plot for each risk group
        plot.KM <- list()
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

          if (name.FR=="incr.FR" & name.TE=="Small.TE"){
            # make KM plot
            KM.i <- survival::survfit(S.h ~ z, subset = groups==risk.group.i)
            plot.KM.i <- survminer::ggsurvplot(KM.i,
                                               data=data.frame(time, z, x),
                                               palette=c("#04bca2", "#e46cf4"),
                                               xlim=c(0, horizon),
                                               break.x.by = 2,
                                               legend="none")
            control.event.rate <- (1-summary(KM.i, times=horizon, extend=TRUE)$surv[1])*100
            risk.control <- c(risk.control, sprintf("%.0f", control.event.rate))

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
                    plot.background=element_rect(fill='transparent', color=NA)) +
              annotate("text", x=ifelse(control.event.rate<80, 0, horizon-2),
                       y=ifelse(control.event.rate<80, 0, 1), size=4, hjust=0,
                       label=paste("CER:", sprintf("%.0f", control.event.rate), "%"))

            # remove y-axis if not risk group 1
            if (risk.group.i!=1){
              plot.KM.i$plot <- plot.KM.i$plot +
                theme(axis.text.y=element_blank(),
                      axis.ticks.y=element_blank(),
                      axis.title.y=element_blank(),
                      axis.line.y=element_blank())
            }
            plot.KM[[risk.group.i]] <- plot.KM.i$plot
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

        # save ARD and dRMST results
        df.ARD.dRMST <- rbind(df.ARD.dRMST,
                              cbind(rep(name.FR, 4),
                                    rep(name.TE, 4),
                                    rep(event.rate, 4),
                                    rep(C.index, 4),
                                    1:4, ARD, dRMST))

        if (name.FR=="incr.FR" & name.TE=="Small.TE"){
          # KM plots
          KM.plots <- ggpubr::ggarrange(plot.KM[[1]], plot.KM[[2]],
                                        plot.KM[[3]], plot.KM[[4]],
                                        ncol=4, nrow=1, align="h")

          # ARD and dRMST plot
          df.ARD.dRMST.plot.i <- data.frame(risk.group=1:4, ARD, dRMST)
          s <- max(df.ARD.dRMST.plot.i$dRMST)/max(df.ARD.dRMST.plot.i$ARD * 100)
          ARD.dRMST.plot <- ggplot2::ggplot(data=df.ARD.dRMST.plot.i,
                                            ggplot2::aes(x=as.numeric(risk.group)))+
            ggplot2::geom_line(aes(y=ARD * 100), col="#AD002AFF", alpha=0.5) +
            ggplot2::geom_point(aes(y=ARD * 100), size=3, shape=18, col="#AD002AFF") +
            ggplot2::geom_line(aes(y=dRMST / s), col="#42B540FF", alpha=0.5) +
            ggplot2::geom_point(aes(y=dRMST / s), size=3, shape=18, col="#42B540FF") +
            ggplot2::scale_y_continuous(
              name=paste0(horizon, "-year ARD in %"),
              sec.axis=ggplot2::sec_axis(~ . * s,
                                         name=paste0(horizon, "-year \u0394RMST in years"))) +
            ggplot2::scale_x_discrete(limits=as.factor(1:4)) +
            ggplot2::xlab("Risk group") +
            ggplot2::theme_bw() +
            ggplot2::theme(axis.title.y=ggplot2::element_text(color="#AD002AFF"),
                           axis.title.y.right=ggplot2::element_text(color="#42B540FF"),
                           axis.line=element_line(colour="black"),
                           plot.background=element_rect(fill='transparent', color=NA),
                           panel.background=element_rect(fill='transparent'),
                           panel.border=element_blank(),
                           panel.grid.minor=ggplot2::element_blank(),
                           panel.grid.major=ggplot2::element_blank(),
                           plot.margin=unit(c(0, 0, 0, 0.8), "cm"))

          # save plot
          ggsave(file=paste0(file.path, "Panel Figures/Fig",
                             ifelse(name.TE=="Small.TE", "2", "3"),
                             LETTERS[nr.plot], ".",
                             name.FR, ".",
                             sens.df$name.bh[sens.df$bh==bh], ".",
                             sens.df$name.C[sens.df$beta==beta], ".",
                             name.TE, ".png"),
                 plot=ggpubr::ggarrange(KM.plots, ARD.dRMST.plot,
                                        nrow=2, ncol=1, heights=c(2, 1), align="h"),
                 width=10, height=6, dpi=300)
        }
      }
    }
  }
  settings.df <- data.frame(settings.df)
  colnames(settings.df) <- c("Outcome", paste("Setting", 1:(ncol(settings.df)-1)))
  openxlsx::write.xlsx(data.frame(settings.df),
                       colNames=FALSE,
                       file=paste0(file.path, "settings.xlsx"))
  df.ARD.dRMST <- as.data.frame(df.ARD.dRMST)
  colnames(df.ARD.dRMST) <- c("FR", "HR", "OR",
                              "Cindex", "risk.group",
                              "ARD", "dRMST")
  save(df.ARD.dRMST,
       file=paste0(file.path, "df.ARD.dRMST.Rdata"))
} else{
  load(paste0(file.path, "df.ARD.dRMST.Rdata"))
}

# create plots
for (name.FR in c("const.FR", "incr.FR")){
  for (name.TE in c("Small.TE", "Large.TE")){
    # numeric TE
    TE <- ifelse(name.TE=="Small.TE", log(0.8), log(0.5))

    # initialize data frame
    df.ARD.dRMST.select <- df.ARD.dRMST |>
      filter(FR==name.FR & HR==name.TE) |>
      select(-HR, -FR) |>
      mutate(setting=rep(1:9, each=4))

    # plot titles
    event.rate <- df.ARD.dRMST.select[seq(1, 36, by=4), "OR"]
    Cindex <- df.ARD.dRMST.select[seq(1, 36, by=4), "Cindex"]
    titles <- paste0(LETTERS[1:9], ". ",
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
    # ordered.dRMST <- sort(c(df.ARD.dRMST.select |>
    #                           filter(setting==5) |>
    #                           select(dRMST) |>
    #                           mutate_all(as.numeric))$dRMST)
    # scarce.treat <- ordered.dRMST[3]+diff(ordered.dRMST[3:4])/2
    # costly.treat <- ordered.dRMST[2]+diff(ordered.dRMST[2:3])/2
    # inexpensive.treat <- ordered.dRMST[1]+diff(ordered.dRMST[1:2])/2

    # make plot
    df <- data.frame(sapply(df.ARD.dRMST.select, as.numeric))
    s <- max(df$dRMST)/max(df$ARD * 100)
    if (name.FR=="const.FR" & name.TE=="Small.TE"){
      red.panels <- c(6, 8, 9)
    } else if (name.FR=="const.FR" & name.TE=="Large.TE"){
      red.panels <- c(8, 9)
    } else if (name.FR=="incr.FR" & name.TE=="Small.TE"){
      red.panels <- c(6, 8, 9)
    } else if (name.FR=="incr.FR" & name.TE=="Large.TE"){
      red.panels <- c(9)
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
    ggsave(file=paste0(file.path, name.FR, ".", name.TE, ".png"),
           plot=plot,
           width=8, height=6, dpi=300)
  }
}
