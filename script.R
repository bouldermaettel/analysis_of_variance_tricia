#################################################################################
#
# Equivalence testing between different Risk analysis methods using Boostrapping
#
#
################################################################################## 
# 
# Matthias Mueller
# Date: 20.05.2023
#
#################################################################################

# readClipboard()
remove(list=ls())

library(readxl)
library(dplyr)
library(tidyverse)

numeric_cols <- c("S", "D", "P")
mapping <- list(
  "S" = c('hard_severity', 'hard_severity_new'),
  "D" = c('hard_detectable', 'hard_detectable_new'),
  "S" = c('hard_probability', 'hard_probability_new'))

column_names <- colnames(read_excel("Varianzanalyse_v3_RESULTS.xlsx", 1))

col_types <- rep("guess", length(column_names))
col_types[column_names %in% numeric_cols] <- "numeric"
data_tot <- read_excel("Varianzanalyse_v3_RESULTS.xlsx", 1, col_types = col_types)
col_types_tricia <- rep('numeric', 9)
tricia_data <- read_excel('tricia_data.xlsx', 1, col_types = col_types_tricia)
head(tricia_data)
str(tricia_data)

data_tot$Testperson <- factor(data_tot$Testperson)
head(data_tot)
str(data_tot)


# define replicates as factor



# define alpha's for prediction and for the confidence interval of the mean difference
alpha_pred <- 0.05 
alpha_conf  <- 0.1

# number of future n number of within series differencies
n_test <- 5

tasks <- c("S", "D", "P")

# Estimation of standard deviation using bootstrap

## To test the script
# i <- 63793
# k <- "S"

for (k in unique(tasks)) {
  
  diff_vector <- c()

  # Loop over VK_IDs to get cases
  for (i in na.omit(unique(data_tot$VK_ID))) {
  
    id <- na.omit(data_tot[data_tot$VK_ID == i, ])
    
    case <- id %>% pull(k) %>% as.numeric()

    # Get all combinations of differences
    diff <- combn(case, 2, function(x) diff(x)) 
    
    diff_vector <- c(diff_vector, diff)

  }
    
    #####################################################################
    # using bootstrap to define standard deviation of mean differences  #
    #####################################################################
    
    # Only the sign is random (for the bootstrap), so work with absolute values and randomly allocate a plus or a minus
    abs_diff <- abs(diff_vector)
  
    # qqnorm(diff)
    # qqnorm(abs_diff)
    # qqline(abs_diff)
    
    mu_abs <- mean(abs_diff)
    sd_abs <- sd(abs_diff)
    
    # calculate mean and CI for abs diff (should be within the AC!)
    abs_diff_CI <- mu_abs + qt(c(0.05,0.95), length(abs_diff)-1)*sd_abs
    
    # define empty vectors
    sd_diff <- NA
    mean_diff <- NA
    diff_numb <- NA
    set.seed(123)
    for (z in 1:1e2){
      sign <- rbinom(length(abs_diff),1,0.5)*2-1 # binomial distribution with 50% 0 and 50% 1: factor 2 and plus 1 to transform to -1 and 1
      sign_diff <-sign*abs_diff
      sd_diff[z] <- sd(sign_diff)
      mean_diff[z] <- mean(sign_diff)
       
      if (z == 1) {
        
        diff_numb <- sign_diff
      } else {
        diff_numb <- append(diff_numb, sign_diff)
      }
    }
    
    sd(diff_numb)

    # calculate sd based on mean_diff and multipy with sqrt(n)
    sd_diff_boot <- sd(mean_diff)*sqrt(length(abs_diff))
    
    sd_diff_boot2 <- mean(sd_diff)
    # qqnorm(sd_diff)
    sd_diff_boot
    mean_diff_boot <- mean(mean_diff)
    
    # qqnorm(mean_diff)
    sd(mean_diff)
    mean_diff_boot # must be close to zero!
    
    
    ## do the sampling in tricia_data
    sd_diff_tricia <- NA
    mean_diff_tricia <- NA
    diff_numb_tricia <- NA
    set.seed(123)
    for (z in 1:1e2){
      tricia_sample <- tricia_data[sample(nrow(tricia_data), length(abs_diff)), ]
      cols = mapping[k]
      sign_diff <- tricia_data[cols[1]] - tricia_data[cols[2]]
      sd_diff_tricia[z] <- sd(sign_diff)
      mean_diff_tricia[z] <- mean(sign_diff)
      
      if (z == 1) {
        
        diff_numb_tricia <- sign_diff
      } else {
        diff_numb_tricia <- append(diff_numb_tricia, sign_diff_tricia)
      }
    }
    
    sd(diff_numb_tricia)
    
    # calculate sd based on mean_diff and multipy with sqrt(n)
    sd_diff_boot_tricia <- sd(mean_diff_tricia)*sqrt(length(abs_diff))
    
    sd_diff_boot2_tricia <- mean(mean_diff_tricia)
    # qqnorm(sd_diff)
    sd_diff_boot_tricia
    mean_diff_boot_tricia <- mean(mean_diff_tricia)
    
    # qqnorm(mean_diff)
    sd(mean_diff_tricia)
    mean_diff_boot_tricia # must be close to zero!
    
    
    
    ############################################################################################
    # The within standard deviation approach! do not forget to define Testperson as factor!   #
    ###########################################################################################
    
    data_anova <- data_tot[!is.na(data_tot[[k]]), ]
    
    summary(fit <- lm (data_anvoa[[k]] ~ Testperson, data = data_anova))
    anova(fit)
   
    sd_within <- sqrt(anova(fit)["Residuals", "Mean Sq"])
    sd_within*sqrt(2) # is close to sd
    
    # simulate data sets to calculate differencies
    set.seed(123)
    diffs <- rnorm(1e5,0,sd_within)-rnorm(1e5,0,sd_within)
    mean_diff <- mean(diffs)
    sd <- sd(diffs)  # is equivalent to the sd_within*sqrt(2)
    
    qqnorm(diffs)
    qqline(diffs)
    
    #######################################################################################
    # The "exact non-parametric approach"                                                 #
    #######################################################################################
    
    # define combinations
    comb4 <- expand.grid(x1 = c(1,-1), x2 = c(1,-1), x3 = c(1,-1), x4 = c(1,-1))
    comb5 <- expand.grid(x1 = c(1,-1), x2 = c(1,-1), x3 = c(1,-1), x4 = c(1,-1), x5 = c(1,-1))
    comb6 <- expand.grid (x1 = c(1,-1), x2 = c(1,-1), x3 = c(1,-1), x4 = c(1,-1), x5 = c(1,-1), x6 = c(1,-1))
    
    combs <- list(comb4,comb5, comb6)
    df <- data.frame(NA,NA,NA, NA, NA, NA)
    df <- df[,1:length(abs_diff)]
    CIs <- data.frame(NA,NA) # used for consistency check!
    CIs2 <- data.frame(NA,NA) # used for consistency check!
    
    for (l in 1: nrow(combs[[length(abs_diff)-3]])) {
      df[l,] <-  abs_diff*combs[[length(abs_diff)-3]][l,]
      CIs[l,] <- mean(as.numeric(df[l,])) + qt(c(alpha_conf/2,1-alpha_conf/2),length(abs_diff)-1)*sd(df[l,])/sqrt(length(abs_diff))
      CIs2[l,] <- c(t.test(as.numeric(df[l,]), conf_level = 0.9)$conf.int[1],t.test(as.numeric(df[l,]))$conf.int[2])

    }
      
    CIlower <- quantile(CIs[,1],p = 0.025)*sqrt(6)/sqrt(3)
    CIupper <- quantile(CIs[,2],p = 0.975)
    
    # calculate means, sd of means and sd of values
      means <- apply(df,1,mean )
      sd_means <- sd(means)
      sd_permut <- sd_means*sqrt(length(abs_diff))
      
      qqnorm(means)
      qqline(means)
      
      #calculate min max CI boundaries (too see if the worst case CI is in the AC)
      CImin <- min(CIs[,2])*2
      CImax <- max(CIs[,2])*2
      
      CIlower_avg <- mean(CIs[,1])
      CIupper_avg <- mean(CIs[,2])
      
      CIlower2_avg <- mean(CIs2[,1])
      CIupper2_avg <- mean(CIs2[,2])
      
      # calculate (1-alpha_conf)-confidence interval for the mean difference of the future measures (n.test)
      CI <- qt(c(alpha_conf/2,1-alpha_conf/2),n_test-1)*sd/sqrt(n_test)
      CIupper <- CI[2]
      CIlower <- CI[1]
      
      # calculate (1-alpha_pred)-predicion intervals from original samples
      pred <-  sqrt(1/length(abs_diff) + 1/n_test)*sd*qt(c(alpha_pred/2,1-alpha_pred/2),length(abs_diff)-1)
      
      # # calculation acc_ to pharmacopoe and according to A149 (not PI for mean but for single value)
      # pred_pharm <-  sqrt(1/length(abs_diff) + 1)*sd_diff_boot*qt(c(0.05,0.95),length(abs_diff)-1) # is used as AC, so compare it to the AC
      
      # calculate AC (PI + each side 1/2 (1-alpha_conf)-%-CI length)
      AC   <- pred + CI
      ACupper <- AC[2]
      AClower <- AC[1]
      AC
      
    
    setwd("/plots")
    
    png(filename=paste("Equivalence Limits_", i, k, ".png" ,sep = ""), 
        type="cairo",
        units="in", 
        width=6, 
        height=10, 
        pointsize=12, 
        res=200)
    
    a <- (AC[2]-AC[1])/10
    #Equivalence limit calculation
    plot(c(AC[1]-a,AC[2]+a),  type = "n",axes = F, xlim = c(2,18), 
         main = paste("Equivalence Limits for", i, "and" , k, sep = " "), ylab = data.sample$Units[1])
    axis(2, labels = c(round(AC[1],3),round(mean(pred),0),round(AC[2],3)), 
         at = c(AC[1],0,AC[2]))
    box()
    
    abline(AC[1],0, col = "red", lwd = 2)
    abline(AC[2],0, col = "red", lwd = 2)
    abline(0,0, col = "red", lwd = 2, lty = "dashed")
    
    # segments(8,CImin,12,CImin, lwd = 2, col = "blue")
    # segments(8,CImax,12,CImax, lwd = 2, col = "blue")
    # segments(10,CImin,10,CImax, lwd = 2)
    # 
    # segments(8,CIlower.avg,12,CIlower.avg, lwd = 2)
    # segments(8,CIupper.avg,12,CIupper.avg, lwd = 2)
    # 
    text (x = 5, y = AC[2]-a, labels = c("Upper acceptance limit"), col = "red")
    text (x = 5, y = AC[1]+a, labels = c("Lower acceptance limit"), col = "red")
    
    legend("topright", legend =c(paste("sd bootstrap approach     =" , round(sd_diff.boot, 3)), 
                                 paste("sd within approach            =", round(sd, 3)),
                                 paste("sd permutation approach =", round(sd.permut,3))))
    
    dev.off()
  }
    
}
