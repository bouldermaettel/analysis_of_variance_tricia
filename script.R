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
library(ggplot2)
# library(tidyverse)

numeric_cols <- c("S", "D", "P", "RBC")
mapping <- list(
  "S" = c('hard_severity', 'hard_severity_new'),
  "D" = c('hard_detectable', 'hard_detectable_new'),
  "P" = c('hard_probability', 'hard_probability_new'),
  "RBC" = c('old_hard_rat', "new_hard_rat"))

column_names <- colnames(read_excel("Varianzanalyse_v4_RESULTS.xlsx", 1))

col_types <- rep("guess", length(column_names))
col_types[column_names %in% numeric_cols] <- "numeric"
data_tot <- read_excel("Varianzanalyse_v4_RESULTS.xlsx", 1, col_types = col_types)

# remove case with no rat
data_tot <- data_tot[complete.cases(data_tot[, "S"]), ]

# # ausschluss wenn VK_NR 65106 ist (kein RAT)
# data_tot <- data_tot[data_tot$VK_ID != 65106,]

col_types_tricia <- rep('numeric', 9)
tricia_data <- read_excel('tricia_data.xlsx', 1, col_types = col_types_tricia)
head(tricia_data)
str(tricia_data)

data_tot$Testperson <- factor(data_tot$Testperson)
head(data_tot)
str(data_tot)


# define alpha's for pblackiction and for the confidence interval of the mean difference
alpha_pblack <- 0.05 
alpha_conf  <- 0.1

# number of future n number of within series differencies
n_test <- 10

tasks <- c("S", "D", "P", "RBC")

# Estimation of standard deviation using bootstrap

## To test the script assign the i
# i <- 63793
# k <- "RBC"
# k <- "P"

for (k in unique(tasks)) {
  
  diff_vector <- c()

  # Loop over VK_IDs to get cases
  for (i in na.omit(unique(data_tot$VK_ID))) {
  
    id <- na.omit(data_tot[data_tot$VK_ID == i, c(k, "VK_ID")])
    
    case <- id %>% pull(k) %>% as.numeric()
    vk_nr <- (id %>% pull('VK_ID') %>% as.numeric())[1]

    # Get all combinations of differences
    diff <- combn(case, 2, function(x) diff(x)) 
    print(c(diff, paste0('id: ', vk_nr)))
    
    diff_vector <- c(diff_vector, diff)

  }
    
    #####################################################################
    # using bootstrap to define standard deviation of mean differences  #
    #####################################################################
    
    # Only the sign is random (for the bootstrap), so work with absolute values and randomly allocate a plus or a minus
    abs_diff <- abs(diff_vector)
  
    plot(abs_diff)
  
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
    CIs2 <- data.frame(NA,NA)
    set.seed(123)
    for (z in 1:1e3){
      sign <- rbinom(length(abs_diff),1,0.5)*2-1 # binomial distribution with 50% 0 and 50% 1: factor 2 and plus 1 to transform to -1 and 1
      sign_diff <-sign*abs_diff
      sd_diff[z] <- sd(sign_diff)
      mean_diff[z] <- mean(sign_diff)
      CIs2[z,] <- c(t.test(sign_diff, conf_level = 0.9, var.equal = TRUE)$conf.int[1], t.test(sign_diff, conf_level = 0.9, var.equal = TRUE)$conf.int[2])
      
       
      if (z == 1) {
        
        diff_numb <- sign_diff
      } else {
        diff_numb <- append(diff_numb, sign_diff)
      }
    }
    
    sd(diff_numb)
    
    CIlower2 <- mean(CIs2[,1])
    CIupper2 <- mean(CIs2[,2])

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
    CIs_tricia <- data.frame(NA, NA)
    CIs2_tricia <- data.frame(NA, NA)
    cols <- mapping[k][[1]]
    set.seed(123)
    for (z in 1:1e3){
      tricia_sample_selection <- tricia_data[tricia_data[cols[2]] != 0, ]
      tricia_sample <- tricia_sample_selection[sample(nrow(tricia_sample_selection), length(abs_diff)), ]
      sign_diff_df <- tricia_sample[cols[1]] - tricia_sample[cols[2]]
      sign_diff <- sign_diff_df %>% pull(cols[1]) %>% as.numeric()
      sd_diff_tricia[z] <- sd(sign_diff)
      mean_diff_tricia[z] <- mean(sign_diff)
      # CIs_tricia[z,] <- mean(sign_diff) + qt(c(alpha_conf/2,1-alpha_conf/2),length(sign_diff)-1)*sd_diff_tricia[z]/sqrt(length(sign_diff))
      CIs2_tricia[z,] <- c(t.test(sign_diff, conf_level = 0.9, var.equal = TRUE)$conf.int[1], t.test(sign_diff, conf_level = 0.9, var.equal = TRUE)$conf.int[2])
      
      
      
      if (z == 1) {
        
        diff_numb_tricia <- sign_diff
      } else {
        diff_numb_tricia <- append(diff_numb_tricia, sign_diff)
      }
    }
    
    # CIlower_tricia <- mean(CIs_tricia[,1])
    # CIupper_tricia <- mean(CIs_tricia[,2])
    
    CIlower2_tricia <- mean(CIs2_tricia[,1])
    CIupper2_tricia  <- mean(CIs2_tricia[,2])
    
    sd(diff_numb_tricia)
    
    # Estimate density
    df_sme <- data.frame(results = diff_numb)
    df_tricia <- data.frame(results = diff_numb_tricia)
    df_sme$group <- "SME"
    df_tricia$group <- "TRICIA"
    
    df <- rbind(df_sme, df_tricia)
    
    
    ggplot(df, aes(x = results, fill = group)) +
      geom_histogram(position = "dodge", alpha = 0.6, binwidth = 1) +
      labs(title = paste0("Histogram  ", cols[1]), x = "Values", y = "Density") +
      scale_fill_manual(values = c('SME' = "grey", 'TRICIA' = "black")) +
      theme_minimal() +
      scale_x_continuous(breaks = sort(unique(c(diff_numb_tricia, diff_numb, -1*diff_numb_tricia, -1*diff_numb))))

    ggsave(paste0("./plots/Histogram", cols[1], ".jpeg"), device = "jpeg", plot = last_plot(), dpi = 300)


    # # Create density plot
    # p <- ggplot(df, aes(x = results, fill = group)) +
    #   geom_density(alpha = 0.8, color = "black") +
    #   labs(title = paste0("Density Plot ", cols[1]), x = "Values", y = "Density") +
    #   scale_fill_manual(values = c("grey", "black")) + 
    #   scale_x_continuous(breaks = sort(unique(diff_numb_tricia))) +
    #   geom_vline(xintercept = CIlower2, color = "grey", linetype = "dashed", linewidth = 1.5) +
    #   geom_vline(xintercept = CIupper2, color = "grey", linetype = "dashed", linewidth = 1.5) +
    #   geom_vline(xintercept = CIlower2_tricia, color = "black", linetype = "dashed", linewidth = 1.5) +
    #   geom_vline(xintercept = CIupper2_tricia, color = "black", linetype = "dashed", linewidth = 1.5)
    #   
    # png(paste0("./plots/density_plot_", cols[1], ".png"))
    # print(p)
    # dev.off()
    # Save the density plot as a PNG file
    # ggsave(paste0("density_plot", cols[1], ".png"), plot = last_plot())
    
  # Create density plot
    max_CI <- ceiling(max(abs(CIlower2), abs(CIlower2_tricia), abs(CIupper2), abs(CIupper2_tricia))) 
    
    ggplot(df, aes(x = results, fill = group)) +
      # geom_histogram(position = "dodge", alpha = 0.2, binwidth = 2) +
      xlim(-max_CI, max_CI) +
      labs(title = paste0("CI Plot ", cols[1]), x = "Values", y = "Density") +
      # scale_fill_manual(values = c('SME' = "grey", 'TRICIA' = "black")) +
      theme_minimal() +
      geom_vline(xintercept = CIlower2, color = "grey", linetype = "dashed", linewidth = 1.5) +
      geom_vline(xintercept = CIupper2, color = "grey", linetype = "dashed", linewidth = 1.5) +
      geom_vline(xintercept = CIlower2_tricia, color = "black", linetype = "dashed", linewidth = 1.5) +
      geom_vline(xintercept = CIupper2_tricia, color = "black", linetype = "dashed", linewidth = 1.5)
      
    # ggplot(df, aes(x = results, fill = group)) +
    #   geom_density(alpha = 0.4, color = "black") + 
    #   labs(title = paste0("Density Plot ", cols[1]), x = "Values", y = "Density") +
    #   scale_fill_manual(values = c("grey", "black")) + 
    #   xlim(-max_CI, max_CI) +
    #   geom_vline(xintercept = CIlower2, color = "grey", linetype = "dashed", linewidth = 1.5) +
    #   geom_vline(xintercept = CIupper2, color = "grey", linetype = "dashed", linewidth = 1.5) +
    #   geom_vline(xintercept = CIlower2_tricia, color = "black", linetype = "dashed", linewidth = 1.5) +
    #   geom_vline(xintercept = CIupper2_tricia, color = "black", linetype = "dashed", linewidth = 1.5)
    # 
    ggsave(paste0("./plots/CI_", cols[1], ".jpeg"), device = "jpeg", plot = last_plot(), dpi = 300)
    
    # calculate sd based on mean_diff and multipy with sqrt(n)
    sd_diff_boot_tricia <- sd(mean_diff_tricia)*sqrt(length(abs_diff))
    
    sd_diff_boot2_tricia <- mean(mean_diff_tricia)
    # qqnorm(sd_diff)
    sd_diff_boot_tricia
    mean_diff_boot_tricia <- mean(mean_diff_tricia)
    
    # qqnorm(mean_diff)
    sd(mean_diff_tricia)
    mean_diff_boot_tricia # must be close to zero!
}