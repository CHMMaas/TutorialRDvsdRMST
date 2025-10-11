#######
####### Perform risk-stratified treatment effect estimation for RCTs
#######
rm(list = ls(all.names = TRUE))
if(!is.null(dev.list())) dev.off()
cat('\014')

# load packages
library(survival)
library(patchwork)
library(ggplot2)
library(dplyr)
# remotes::install_github("CHMMaas/PredictionTools")
library(PredictionTools)

# initialize
set.seed(1)
new.results <- TRUE
file.path <- "Z:/."
nr.RCTs <- 1:2
g <- 4 # number of risk groups to stratify
new.imputation <- FALSE
data.names <- list(RCT1=list(name="SOLVD_P2",
                              name.long="SOLVD Prevention",
                              name.y="EPX",
                              name.time="EPXTIME",
                              names.X=c("age",
                                        "SEF25", # BUN
                                        "SBF38", # NYHA
                                        "diuretic",
                                        "female",
                                        "white",
                                        "cabg",
                                        "copd",
                                        "afib",
                                        "sysbp",
                                        "diabp",
                                        "pulse",
                                        "ace",
                                        "bblocker",
                                        "SEF26", # creatinine
                                        "SEF24", # potassium
                                        "SEF22Z1"), # total white blood count
                              name.W="treat"),
                   RCT2=list(name="SOLVD_T",
                              name.long="SOLVD Intervention",
                              name.y="EP1",
                              name.time="FUTIME",
                              names.X=c("SEF_AGE",
                                        "female",
                                        "wtkg",
                                        "SBF34", # heart rate
                                        "smk",
                                        "db",
                                        "mi",
                                        "edema",
                                        "copd",
                                        "SBF38", # NYHA
                                        "SBF27Z1", # ejection fraction percentage
                                        "ace",
                                        "bblocker",
                                        "creat2",
                                        "SBF35Z1", # sbp
                                        "SBF35Z2", # dbp
                                        "SEF23", # sodium
                                        "SEF25"), # BUN
                              name.W="treat"))

if (new.results){
  # orange: #EC7F12
  # purple: #9161BD
  # blue: #0D88AF
  plots <- data.frame(c())
  low.high.risk.df <- data.frame(c())
  low.high.risk.df2 <- data.frame(c())
  descriptives <- data.frame(c())
  for (nr.RCT in nr.RCTs){
    # load original data
    data.name <- eval(parse(text=paste0("data.names$RCT", nr.RCT)))
    original.data <- read.csv(paste0(file.path, "Data/", data.name$name, ".csv"))

    if (data.name$name.long=="DPP lifestyle"){
      # select lifestyle or placebo
      original.data <- original.data[-which(original.data$met==1),]
    } else if (data.name$name.long=="DPP metformin"){
      # select metformin or placebo
      original.data <- original.data[-which(original.data$life==1),]
    }

    if (data.name$name.long=="MTOPS doxazosin"){
      # select doxazosin or placebo
      original.data <- original.data[-which(original.data$treat2==1 | original.data$treat3==1),]
    } else if (data.name$name.long=="MTOPS finasteride"){
      # select finasteride or placebo
      original.data <- original.data[-which(original.data$treat1==1 | original.data$treat3==1),]
    } else if (data.name$name.long=="MTOPS combination"){
      # select combination or placebo
      original.data <- original.data[-which(original.data$treat1==1 | original.data$treat2==1),]
    }

    # select covariates
    y <- original.data[, data.name$name.y]
    X <- original.data[, data.name$names.X]
    W <- original.data[, data.name$name.W]

    # define time
    time <- original.data[, data.name$name.time]
    if (data.name$name=="DPP_DPPOSv3"){
      time <- time/12
    } else if (data.name$name!="ACCORD2" & data.name$name!="DCCT_P" & data.name$name!="DCCT_S"){
      time <- time/365.25
    }

    # omit observations that have no observed outcomes
    if (sum(is.na(y))>0){
      cat("Number of missings in y:", sum(is.na(y)), '\n')
      original.data <- original.data[-which(is.na(y)),]
    }

    # single imputation
    if (new.imputation){
      imputed.data <- mice::complete(mice::mice(cbind(y, X, W, time), m=1, seed=1, print=FALSE))
      save(imputed.data, file=paste0(file.path, "Data/",
                                     data.name$name, ".",
                                     data.name$name.W, ".Rdata"))
    } else{
      load(file=paste0(file.path, "Data/",
                       data.name$name, ".",
                       data.name$name.W, ".Rdata"))
    }

    # load imputed data
    cat("Number of observations loaded from", data.name$name, ":", nrow(imputed.data), '\n')

    # 1. select covariates
    y <- imputed.data$y
    time <- imputed.data$time
    W <- imputed.data$W
    X <- as.data.frame(scale(imputed.data[, data.name$names.X]))

    # 2. single imputation
    work.data <- as.data.frame(cbind(y, time, W, X))
    dd <- rms::datadist(work.data)
    options(datadist = 'dd')
    work.data$S <- survival::Surv(time=time, event=y)

    # 3. limit survival to median follow-up time
    # median follow-up time
    time.points <- c(0.1, 0.25, 0.5, 0.75, 0.9)
    q.f.time <- stats::quantile(prodlim::prodlim(prodlim::Hist(time=time, event=y) ~ 1, reverse=TRUE), time.points)
    FU.times <- c(q.f.time$quantile[q.f.time$q==(1-time.points[3])],
                  ifelse(is.na(q.f.time$quantile[q.f.time$q==(1-time.points[5])]),
                         max(time),
                         q.f.time$quantile[q.f.time$q==(1-time.points[5])]))

    # reverse Kaplan-Meier curve
    rev.Surv <- survminer::ggsurvplot(survival::survfit(survival::Surv(time=work.data$time,
                                                                       event=1-work.data$y) ~ W,
                                                        data=work.data),
                                      risk.table=TRUE,
                                      title=data.name$name.long)
    rev.Surv$plot <- rev.Surv$plot+
      ggplot2::geom_hline(yintercept=0.1)
    show(rev.Surv)

    for (FU.time in FU.times){
      cat("Follow-up time:", FU.time, "years \n")
      work.data$S.FU <- work.data$S
      work.data$S.FU[work.data$S[, 1] > FU.time, 2] <- 0
      work.data$S.FU[work.data$S[, 1] > FU.time, 1] <- FU.time

      # Occurrence of the event (taking into account censoring)
      KM.all <- summary(survival::survfit(S.FU ~ 1, data=work.data), times=FU.time)

      # 4. fit survival model
      metrics <- calculate.RD.dRMST(S=work.data$S.FU,
                                     W=work.data$W,
                                     horizon=FU.time)

      # 5. Kaplan-Meier at median follow-up time
      RD <- metrics$RD
      RD.lower <- sprintf("%.2f", RD-1.96*metrics$RD_se)
      RD.upper <- sprintf("%.2f", RD+1.96*metrics$RD_se)
      LRT.p <- metrics$LRT.p

      # estimate hazard ratio
      cox <- rms::cph(S.FU ~ W, data=work.data)
      HR <- as.numeric(exp(coef(cox)[1]))
      HR.lower <- sprintf("%.2f", as.numeric(exp(confint(cox)[1,1])))
      HR.upper <- sprintf("%.2f", as.numeric(exp(confint(cox)[1,2])))
      test.PH <- cox.zph(cox)$table["W", "p"]

      # 5. estimate RMST
      dRMST <- metrics$dRMST
      dRMST.lower <- sprintf("%.2f", dRMST-1.96*metrics$dRMST_se)
      dRMST.upper <-  sprintf("%.2f", dRMST+1.96*metrics$dRMST_se)

      # 6. fit risk model for risk stratification
      risk.formula <- eval(parse(text=paste("S.FU ~ ",
                                            paste(data.name$names.X, collapse=" + "),
                                            sep="")))
      risk.model <- rms::cph(risk.formula, data=work.data, se.fit=TRUE, x=TRUE, y=TRUE)
      C.index <- concordance(risk.model)$concordance

      if (FU.time==FU.times[2]){
        # save results
        descriptives["N", data.name$name.long] <- length(y)
        descriptives["KM", data.name$name.long] <- paste0(sprintf("%.1f", KM.all$surv*100), " [",
                                                          sprintf("%.1f", KM.all$lower*100), "; ",
                                                          sprintf("%.1f", KM.all$upper*100), "]")
        descriptives["p-value test PH assumption", data.name$name.long] <- sprintf("%.3f", test.PH)
        descriptives["RD", data.name$name.long] <- paste0(sprintf("%.2f", RD), " [",
                                                           RD.lower, "; ",
                                                           RD.upper, "]")
        descriptives["p-value log rank test", data.name$name.long] <- sprintf("%.3f", LRT.p)
        descriptives["dRMST", data.name$name.long] <- paste0(sprintf("%.2f", dRMST), " [",
                                                             dRMST.lower, "; ",
                                                             dRMST.upper, "]")
        descriptives["Median [IQR] follow-up", data.name$name.long] <- paste0(sprintf("%.1f", q.f.time$quantile[q.f.time$q==0.5]), " [",
                                                                              sprintf("%.1f", q.f.time$quantile[q.f.time$q==0.75]), "; ",
                                                                              sprintf("%.1f", q.f.time$quantile[q.f.time$q==0.25]), "]")
        descriptives["90th percentile", data.name$name.long] <- sprintf("%.1f", FU.time)
        descriptives["Max time", data.name$name.long] <- sprintf("%.1f", max(time))
        descriptives["C-index", data.name$name.long] <- sprintf("%.2f", as.numeric(Hmisc::rcorr.cens(predict(risk.model), y)["C Index"]))
      }

      # make predictions
      f.basehaz <- survival::basehaz(rms::cph(risk.formula, data=work.data, se.fit=TRUE, x=TRUE, y=TRUE))
      h0 <- f.basehaz$hazard[f.basehaz$time==max(f.basehaz$time[f.basehaz$time<=FU.time])]
      lp <- predict(risk.model, type="lp")
      surv.prob.risk.model <- 1-exp(-h0*exp(lp))

      # 7. stratify to risk group
      work.data$lp <- lp # - mean(lp) # center the LP?
      work.data$quartile.nr <- as.numeric(cut(work.data$lp,
                                              breaks=quantile(work.data$lp,
                                                              probs=seq(0, 1, by=1/g)),
                                              include.lowest=TRUE))
      table(work.data$quartile.nr)

      # test constant HR
      p.value.test.HR <- summary(survival::coxph(S ~ quartile.nr*W,
                                                 data=work.data))$coefficients[3, "Pr(>|z|)"]

      # 8. fit KM, Cox, and dRMST in each risk group
      KM.quartiles <- c()
      plot.KM <- list()
      table.KM <- list()
      HR.Q <- c()
      HR.Q.se <- c()
      RD.Q <- c()
      RD.Q.se <- c()
      dRMST.Q <- c()
      dRMST.Q.se <- c()
      for (i in 1:g){
        # select data from quartile
        quartile.data <- work.data[work.data$quartile.nr==i, ]
        quartile.data$S.q <- survival::Surv(time=quartile.data$time, event=quartile.data$y)
        quartile.data$S.q[quartile.data$S.q[, 1] > FU.time, 2] <- 0
        quartile.data$S.q[quartile.data$S.q[, 1] > FU.time, 1] <- FU.time

        # fit survival model
        KM.q <- survival::survfit(S.q ~ W, data=quartile.data)

        # plot.KM
        if (FU.time==FU.times[2]){
          # plot KM
          plot.KM.i <- survminer::ggsurvplot(KM.q,
                                             data=quartile.data,
                                             palette=c("#04bca2", "#e46cf4"),
                                             xlim=c(0, FU.time),
                                             conf.int=TRUE,
                                             legend="none",
                                             title=paste0("N=", nrow(quartile.data)),
                                             risk.table=TRUE,
                                             risk.table.fontsize=3,
                                             risk.table.title=ggplot2::element_blank())
          plot.KM.i$plot <- plot.KM.i$plot +
            theme(plot.title=element_text(hjust=0.5),
                  axis.text.x=element_text(size=12),
                  axis.title.x=element_text(size=12),
                  axis.line.x=element_line(colour="black"),
                  axis.text.y=element_text(size=12),
                  axis.title.y=element_text(size=12),
                  axis.line.y=element_line(colour="black"),
                  panel.grid.major=element_blank(),
                  panel.grid.minor=element_blank())
                  # panel.background=element_rect(fill='transparent'),
                  # plot.background=element_rect(fill='transparent', color=NA))
          plot.KM.i$table <- plot.KM.i$table
            # theme(panel.background=element_rect(fill='transparent'),
            #       plot.background=element_rect(fill='transparent', color=NA))
          if (i!=1){
            plot.KM.i$plot <- plot.KM.i$plot +
              theme(axis.text.y=element_blank(),
                    axis.ticks.y=element_blank(),
                    axis.title.y=element_blank(),
                    axis.line.y=element_blank())
            plot.KM.i$table <- plot.KM.i$table +
              theme(axis.text.y=element_blank(),
                    axis.ticks.y=element_blank(),
                    axis.title.y=element_blank(),
                    axis.line.y=element_blank())
          }
          plot.KM[[i]] <- plot.KM.i$plot


          # table KM
          table.KM[[i]] <- plot.KM.i$table
        }

        # Kaplan-Meier at median follow-up time
        metrics.q <- calculate.RD.dRMST(S=quartile.data$S.q,
                                         W=quartile.data$W,
                                         horizon=FU.time)
        RD.Qi <- metrics.q$RD
        RD.Q <- c(RD.Q, RD.Qi)
        RD.Q.se <- c(RD.Q.se, metrics.q$RD_se)
        RD.Qi.lower <- RD.Qi-1.96*metrics.q$RD_se
        RD.Qi.upper <- RD.Qi+1.96*metrics.q$RD_se

        # estimate dRMST
        dRMST.Qi <- metrics.q$dRMST
        dRMST.Q <- c(dRMST.Q, dRMST.Qi)
        dRMST.Q.se <- c(dRMST.Q.se, metrics.q$dRMST_se)
        dRMST.Qi.lower <- dRMST.Qi-1.96*metrics.q$dRMST_se
        dRMST.Qi.upper <- dRMST.Qi+1.96*metrics.q$dRMST_se

        # estimate hazard ratio
        quartile.cox <- survival::coxph(S.q ~ W, data=quartile.data)
        HR.Qi <- as.numeric(exp(coef(quartile.cox)[1]))
        HR.Q <- c(HR.Q, HR.Qi)
        HR.Q.se <- c(HR.Q.se, summary(quartile.cox)$coefficients[3])
        HR.Qi.lower <- as.numeric(exp(confint(quartile.cox)[1,1]))
        HR.Qi.upper <- as.numeric(exp(confint(quartile.cox)[1,2]))

        # truncate CAST and DCCT prevention upper limit of HR
        HR.Qi.upper <- ifelse(data.name$name.long=="CAST" & i == 1 & FU.time==FU.times[2], 3, HR.Qi.upper)
        HR.Qi.upper <- ifelse(data.name$name.long=="DCCT Prevention" & i == 1, 3, HR.Qi.upper)
      }
      # save for plots
      plots <- rbind(plots, cbind(rep(nr.RCT, 4),
                                  rep(data.name$name.long, 4),
                                  rep(ifelse(FU.time==FU.times[1], "median", "long"), 4),
                                  rep(KM.all$surv*100, 4),
                                  rep(C.index, 4),
                                  1:4,
                                  rep(RD, 4),
                                  RD.Q,
                                  RD.Q-1.96*RD.Q.se,
                                  RD.Q+1.96*RD.Q.se,
                                  rep(HR, 4),
                                  HR.Q,
                                  HR.Q-1.96*HR.Q.se,
                                  HR.Q+1.96*HR.Q.se,
                                  rep(p.value.test.HR, 4),
                                  rep(dRMST, 4),
                                  dRMST.Q,
                                  dRMST.Q-1.96*dRMST.Q.se,
                                  dRMST.Q+1.96*dRMST.Q.se))

      # make KM plots only for long follow-up
      if (FU.time==FU.times[2]){
        # Kaplan-Meier plots
        KM.plots <- ggpubr::ggarrange(plot.KM[[1]], plot.KM[[2]],
                                      plot.KM[[3]], plot.KM[[4]],
                                      ncol=4, nrow=1, align="h")

        # Kaplan-Meier tables
        KM.tables <- ggpubr::ggarrange(table.KM[[1]], table.KM[[2]],
                                       table.KM[[3]], table.KM[[4]],
                                       ncol=4, nrow=1, align="h")

        # RD and dRMST plot
        df.RD.dRMST <- data.frame(risk.group=1:4, RD=RD.Q, dRMST=dRMST.Q)
        s <- max(df.RD.dRMST$dRMST)/max(df.RD.dRMST$RD * 100)
        RD.dRMST.plot <- ggplot2::ggplot(data=df.RD.dRMST,
                                          ggplot2::aes(x=as.numeric(risk.group)))+
          ggplot2::geom_hline(yintercept=0, color="grey", linetype="dashed") +
          ggplot2::geom_line(aes(y=RD * 100), col="#AD002AFF", alpha=0.5) +
          ggplot2::geom_point(aes(y=RD * 100), size=3, shape=18, col="#AD002AFF") +
          ggplot2::geom_line(aes(y=dRMST / s), col="#42B540FF", alpha=0.5) +
          ggplot2::geom_point(aes(y=dRMST / s), size=3, shape=18, col="#42B540FF") +
          ggplot2::scale_y_continuous(
            name=paste0(sprintf("%.1f", FU.time), "-year RD in %"),
            sec.axis=ggplot2::sec_axis(~ . * s,
                                       name=paste0(sprintf("%.1f", FU.time), "-year \u0394RMST in years"))) +
          ggplot2::scale_x_discrete(limits=as.factor(1:4)) +
          ggplot2::xlab("Risk group") +
          ggplot2::theme_bw() +
          ggplot2::theme(axis.title.y=ggplot2::element_text(color="#AD002AFF"),
                         axis.title.y.right=ggplot2::element_text(color="#42B540FF"),
                         axis.line=element_line(colour="black"),
                         # plot.background=element_rect(fill='transparent', color=NA),
                         # panel.background=element_rect(fill='transparent'),
                         panel.border=element_blank(),
                         panel.grid.minor=ggplot2::element_blank(),
                         panel.grid.major=ggplot2::element_blank(),
                         plot.margin=unit(c(0, 0, 0, 0.8), "cm"))

        # combined plot
        ggplot2::ggsave(filename=paste0(file.path, "Results RCTs Kent/Plots/", data.name$name.long, ".png"),
                        plot=ggpubr::ggarrange(KM.plots, KM.tables,
                                               RD.dRMST.plot,
                                               nrow=3, ncol=1, align="h",
                                               heights=c(1, 0.5, 0.5)) + ggpubr::bgcolor("white"),
                        width=12, height=8, dpi=300)
      }
    }
  }
  # save descriptives
  openxlsx::write.xlsx(descriptives,
                       rowNames=TRUE,
                       file=paste0(file.path, "Results RCTs Kent/Descriptives.xlsx"))

  # make plots
  colnames(plots) <- c("nrData", "Dataset", "Follow-up", "KM", "C.index", "Group",
                       "RD", "RD.Q", "RD.Q.lower", "RD.Q.upper",
                       "HR", "HR.Q", "HR.Q.lower", "HR.Q.upper", "p-value test HR",
                       "dRMST", "dRMST.Q", "dRMST.Q.lower", "dRMST.Q.upper")
  openxlsx::write.xlsx(plots,
                       rowNames=TRUE,
                       file=paste0(file.path, "Results RCTs Kent/Plots.xlsx"))
  save(plots,
       file=paste0(file.path, "Results RCTs Kent/Plots/plots.Rdata"))
} else{
  load(paste0(file.path, "Results RCTs Kent/Plots/plots.Rdata"))
}

# make combined plot
for (length in c("median", "long")){
  # select data frame
  data.plot <- plots[plots$`Follow-up`==length,]

  # plot titles
  event.rate <- as.numeric(unique(data.plot[, "KM"]))
  Cindex <- unique(data.plot[, "C.index"])
  titles <- paste0(unique(data.plot[, "Dataset"]), ", ER: ",
                   sprintf("%.0f", 100-event.rate),
                   "%, C: ", sprintf("%.2f", as.numeric(Cindex)))

  # sort on event rate
  sort.event.rate <- order(-event.rate)
  data.plot$sort.nrData <- rep(sort.event.rate, each=4)
  title.lab <- c("1"=titles[sort.event.rate[1]],
                 "2"=titles[sort.event.rate[2]],
                 "3"=titles[sort.event.rate[3]],
                 "4"=titles[sort.event.rate[4]],
                 "5"=titles[sort.event.rate[5]],
                 "6"=titles[sort.event.rate[6]],
                 "7"=titles[sort.event.rate[7]],
                 "8"=titles[sort.event.rate[8]],
                 "9"=titles[sort.event.rate[9]],
                 "10"=titles[sort.event.rate[10]],
                 "11"=titles[sort.event.rate[11]],
                 "12"=titles[sort.event.rate[12]],
                 "13"=titles[sort.event.rate[13]],
                 "14"=titles[sort.event.rate[14]],
                 "15"=titles[sort.event.rate[15]],
                 "16"=titles[sort.event.rate[16]],
                 "17"=titles[sort.event.rate[17]])
  sorted.title.lab <- title.lab[sort.event.rate]

  # make plot
  df <- data.frame(sapply(data.plot |> group_by(nrData) |>
                            arrange(KM, C.index) |>
                            select(-Dataset, -`Follow-up`), as.numeric))
  s <- max(df$dRMST.Q)/max(df$RD.Q * 100)
  plot <- ggplot2::ggplot(data=df,
                          ggplot2::aes(x=as.numeric(Group)))+
    ggplot2::geom_line(aes(y=RD.Q * 100), col="#04bca2", alpha=0.5) +
    ggplot2::geom_point(aes(y=RD.Q * 100), size=2, shape=18, col="#04bca2") +
    ggplot2::geom_line(aes(y=dRMST.Q / s), col="#e46cf4", alpha=0.5) +
    ggplot2::geom_point(aes(y=dRMST.Q / s), size=2, shape=18, col="#e46cf4") +
    ggplot2::geom_hline(aes(yintercept=RD * 100),
                        col="#04bca2", linetype="dashed") +
    ggplot2::geom_hline(aes(yintercept=dRMST),
                        col="#e46cf4", linetype="dashed") +
    ggplot2::scale_y_continuous(
      name="RD in %",
      sec.axis=ggplot2::sec_axis(~ . * s,
                                 name=paste0("\u0394RMST in years"))) +
    ggplot2::facet_wrap(~ sort.nrData, ncol=3, drop=FALSE,
                        labeller=as_labeller(sorted.title.lab)) +
    ggplot2::scale_x_discrete(limits=as.factor(1:4)) +
    ggplot2::xlab("Risk group") +
    ggplot2::theme_bw() +
    ggplot2::theme(axis.title.y=ggplot2::element_text(color="#04bca2"),
                   axis.title.y.right=ggplot2::element_text(color="#e46cf4"),
                   panel.grid.minor=ggplot2::element_blank(),
                   panel.grid.major=ggplot2::element_blank())
  show(plot)
  ggsave(file=paste0(file.path, "Results RCTs Kent/",
                     length, ".FU.png"),
         plot=plot,
         width=8, height=8, dpi=300)
}

