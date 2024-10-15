# Load necessary library
library(ggplot2)

# baseline risk
baseline_risk <- seq(0.1, 0.9, by = 0.01)
n <- length(baseline_risk)

# constant HR
HR.CHR <- rep(0.66, n)
RD.CHR <- baseline_risk - 1 + (1- baseline_risk)^HR.CHR 
RR.CHR <- (1-(1-baseline_risk)^HR.CHR)/baseline_risk

# add to data frame
data <- data.frame(baseline_risk, 
                   RD.CRD, RR.CRD, HR.CRD,
                   RD.CRR, RR.CRR, HR.CRR,
                   RD.CHR, RR.CHR, HR.CHR)

# Plot RD
plot.RD.CHR <- ggplot(data, aes(x = baseline_risk, y = RD.CHR)) +
  geom_line() +
  scale_x_continuous(limits = c(0, 1),        
                     breaks = seq(0, 1, by=0.25),        
                     labels = seq(0, 1, by=0.25)) +      
  ylim(c(0, 0.2)) +
  theme_classic() +
  ylab("RD") +
  xlab("")

# Plot RR
plot.RR.CHR <- ggplot(data, aes(x = baseline_risk, y = RR.CHR)) +
  geom_line() +
  scale_x_continuous(limits = c(0, 1),              
                     breaks = seq(0, 1, by=0.25),     
                     labels = seq(0, 1, by=0.25)) +      
  ylim(c(0, 1)) +
  theme_classic() +
  ylab("RR") +
  xlab("Baseline risk")

# Plot HR
plot.HR.CHR <- ggplot(data, aes(x = baseline_risk, y = HR.CHR)) +
  geom_line() +
  scale_x_continuous(limits = c(0, 1),              
                     breaks = seq(0, 1, by=0.25),       
                     labels = seq(0, 1, by=0.25)) +      
  ylim(c(0, 1)) +
  theme_classic() +
  ylab("HR") +
  xlab("")

# Show the plots
ggsave(file="C:/Users/carol/OneDrive - Erasmus MC/Projects EMC/Project Tutorial dRMST vs RD/Supplemental Information/Figure 1.png",
       plot=ggpubr::ggarrange(plot.RD.CHR, plot.RR.CHR, plot.HR.CHR,
                              nrow=1, ncol=3),
       width=6, height=2, dpi=300)
