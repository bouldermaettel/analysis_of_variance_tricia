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

data_tot <- read_excel("Varianzanalyse_v3_RESULTS.xlsx", 1)
head(data_tot)
str(data_tot)


# define replicates as factor
data_tot$Testperson <- factor(data_tot$Testperson)


# define alpha's for prediction and for the confidence interval of the mean difference
alpha.pred <- 0.05 
alpha.conf  <- 0.1

# number of future n number of within series differencies
n.test <- 5

tasks = c("S", "D", "P")

# Estimation of standard deviation using bootstrap

## To test the script
# i <- 63793
# k <- s

# Loop over VK_IDs to get cases
for (i in na.omit(unique(data_tot$VK_ID))) {

  case <- na.omit(data_tot[data_tot$VK_ID == i,])
  
  # Loop over tasks (severity, detectability)
  
  for (k in unique(tasks)) {
    
    case_task <- as.numeric(case[[k]])
    
    # Get all combinations of differences
    all_combs <- combn(case_task, 2, function(x) diff(x)) 
    
    # add all combs with different sign
    # all_combs = c(all_combs, -1*(all_combs))
    
    # calculate within differencies
    diff <- na.omit(data.sample$result[data$replicate==1] - data.sample$result[data$replicate==2])
    sd(diff)
    
    #####################################################################
    # using bootstrap to define standard deviation of mean differences  #
    #####################################################################
    
    # Only the sign is random (for the bootstrap), so work with absolute values and randomly allocate a plus or a minus
    abs.diff <- abs(diff)
  
    # qqnorm(diff)
    # qqnorm(abs.diff)
    # qqline(abs.diff)
    
    mu.abs <- mean(abs.diff)
    sd.abs <- sd(abs.diff)
    
    # calculate mean and CI for abs diff (should be within the AC!)
    abs.diff.CI <- mu.abs + qt(c(0.05,0.95), length(abs.diff)-1)*sd.abs
    
    # define empty vectors
    sd.diff <- NA
    mean.diff <- NA
    diff.numb <- NA
    set.seed(123)
    for (z in 1:1e4){
      sign <- rbinom(length(abs.diff),1,0.5)*2-1 # binomial distribution with 50% 0 and 50% 1: factor 2 and plus 1 to transform to -1 and 1
      sign.diff <-sign*abs.diff
      sd.diff[z] <- sd(sign.diff)
      mean.diff[z] <- mean(sign.diff)
       
      if (z == 1) {
        
        diff.numb <- sign.diff
      } else {
        diff.numb <- append(diff.numb, sign.diff)
      }
    }
    

    # qqnorm(diff.numb)
    # qqline(diff.numb)
    # 
    # qqnorm(mean.diff)
    # qqline(mean.diff)
    
    sd(diff.numb)

    # calculate sd based on mean.diff and multipy with sqrt(n)
    sd.diff.boot <- sd(mean.diff)*sqrt(length(abs.diff))
    
    sd.diff.boot2 <- mean(sd.diff)
    # qqnorm(sd.diff)
    sd.diff.boot
    mean.diff.boot <- mean(mean.diff)
    # qqnorm(mean.diff)
    sd(mean.diff)
    mean.diff.boot # must be close to zero!
    
    
    #######################################################################################
    # The within standard deviation approach! do not forget to define series as factor!   #
    #######################################################################################
    
    summary(fit <- lm (result ~  as.factor(serie), data = data.sample))
    anova(fit)
   
    sd.within <- sqrt(anova(fit)["Residuals", "Mean Sq"])
    sd.within*sqrt(2) # is close to sd
    
    # simulate data sets to calculate differencies
    set.seed(123)
    diffs <- rnorm(1e5,0,sd.within)-rnorm(1e5,0,sd.within)
    mean.diff <- mean(diffs)
    sd <- sd(diffs)  # is equivalent to the sd.within*sqrt(2)
    
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
    df <- df[,1:length(abs.diff)]
    CIs <- data.frame(NA,NA) # used for consistency check!
    CIs2 <- data.frame(NA,NA) # used for consistency check!
    
    for (l in 1: nrow(combs[[length(abs.diff)-3]])) {
      df[l,] <-  abs.diff*combs[[length(abs.diff)-3]][l,]
      CIs[l,] <- mean(as.numeric(df[l,])) + qt(c(alpha.conf/2,1-alpha.conf/2),length(abs.diff)-1)*sd(df[l,])/sqrt(length(abs.diff))
      CIs2[l,] <- c(t.test(as.numeric(df[l,]), conf.level = 0.9)$conf.int[1],t.test(as.numeric(df[l,]))$conf.int[2])

    }
      
    CIlower <- quantile(CIs[,1],p = 0.025)*sqrt(6)/sqrt(3)
    CIupper <- quantile(CIs[,2],p = 0.975)
    
    # calculate means, sd of means and sd of values
      means <- apply(df,1,mean )
      sd.means <- sd(means)
      sd.permut <- sd.means*sqrt(length(abs.diff))
      
      qqnorm(means)
      qqline(means)
      
      #calculate min max CI boundaries (too see if the worst case CI is in the AC)
      CImin <- min(CIs[,2])*2
      CImax <- max(CIs[,2])*2
      
      CIlower.avg <- mean(CIs[,1])
      CIupper.avg <- mean(CIs[,2])
      
      CIlower2.avg <- mean(CIs2[,1])
      CIupper2.avg <- mean(CIs2[,2])
      
      # calculate (1-alpha.conf)-confidence interval for the mean difference of the future measures (n.test)
      CI <- qt(c(alpha.conf/2,1-alpha.conf/2),n.test-1)*sd/sqrt(n.test)
      CIupper <- CI[2]
      CIlower <- CI[1]
      
      # calculate (1-alpha.pred)-predicion intervals from original samples
      pred <-  sqrt(1/length(abs.diff) + 1/n.test)*sd*qt(c(alpha.pred/2,1-alpha.pred/2),length(abs.diff)-1)
      
      # # calculation acc. to pharmacopoe and according to A149 (not PI for mean but for single value)
      # pred.pharm <-  sqrt(1/length(abs.diff) + 1)*sd.diff.boot*qt(c(0.05,0.95),length(abs.diff)-1) # is used as AC, so compare it to the AC
      
      # calculate AC (PI + each side 1/2 (1-alpha.conf)-%-CI length)
      AC   <- pred + CI
      ACupper <- AC[2]
      AClower <- AC[1]
      AC
      
    
    setwd("\\\\euchbrnfil01\\Matthias.Mueller$\\MyDocs\\R Scripts\\Brigitte Equivalenztest Differenzen/WithinSeries/plots")
    
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
    
    legend("topright", legend =c(paste("sd bootstrap approach     =" , round(sd.diff.boot, 3)), 
                                 paste("sd within approach            =", round(sd, 3)),
                                 paste("sd permutation approach =", round(sd.permut,3))))
    
    dev.off()
  }
    
}


##### Testing

# Define a vector of numbers
numbers <- c(1, 2, 4, 3)

# Get all combinations of differences
all_combinations <- combn(numbers, 2, function(x) diff(x))

# Print the result
print(all_combinations)

