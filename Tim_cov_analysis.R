cutoff_len = 1000
map_type = "_pairedcov_"
infileext = paste(map_type, "minlen=1000_contig_cov.txt", sep = "")

### libs

library(ggplot2)
library(stringr)
library(modeest)
library(plyr)
library(cowplot)
library(grid)


### function
peakfinder <- function(d){
  dh <- hist(d,plot=FALSE, breaks=1000)
  ins <- dh[["counts"]]
  nbins <- length(ins)
  ss <- which(rank(ins)%in%seq(from=nbins,to=nbins)) ## pick the top 3 intensities
  dh[["mids"]][ss]
}



########################################################################################################################################################################
### DATA

# raw

setwd("data/coverage/")

dat1_Tdi <- read.table (paste("Tdi", infileext, sep = ""), header = T, sep = ',')
head(dat1_Tdi)
### filter by len

dat1_Tdi_filt_1000 <- subset(dat1_Tdi, dat1_Tdi$length >= cutoff_len)
head(dat1_Tdi_filt_1000)


##########################################################################################################################
##### cov of each sample

cov_plot <- function(df, samp, min_cov, max_cov) {
  
  median_cov <- median(eval(parse(text=paste(df,"$",samp, sep = ""))))
  #modal_cov  <- mlv(eval(parse(text=paste(df,"$",samp, sep = ""))), method = "shorth")
  modal_cov  <- mlv(eval(parse(text=paste(df,"$",samp, sep = ""))), method = "hsm")
  q99 <- quantile(eval(parse(text=paste(df,"$",samp, sep = ""))), .99)
  
  
  p1 <- ggplot(eval(parse(text=df)), aes(x=eval(parse(text=samp)))) +
    theme_bw() +
    geom_histogram(color="blue", fill="blue", binwidth=0.2,alpha = 1) +
    xlim(c(min_cov, max_cov)) + 
    geom_hline(yintercept = 0) +
    xlab("Coverage") + 
    ggtitle(paste(samp, " | len >= ", cutoff_len, sep = "")) +
    geom_vline(xintercept = median_cov, linetype="dotted", color = "red", size=1.5) + 
    geom_vline(xintercept = q99, linetype="dotted", color = "orange", size=1.5)+    	
    geom_vline(xintercept = modal_cov, linetype="dotted", color = "green", size=1.5)    	
  
  out_df <- as.data.frame(rbind(c(samp, round(median_cov, digits=2), round(modal_cov, digits=2), round(q99, digits=2))))
  colnames(out_df) <- c("sample", "med_cov", "mod_cov", "q99")
  
  output = list("p1" = p1, "out_df" = out_df)
  return(output)
  
  
}




#################################################################################################################
#################################################################################################################
### Plot cov 

### Tdi
min_cov_1 = 1
max_cov_1 = 100

Tdi_a1 <- plot_grid(
  cov_plot("dat1_Tdi_filt_1000", "Tdi_M_18.3997", min_cov_1, max_cov_1)$p1,
  cov_plot("dat1_Tdi_filt_1000", "Tdi_M_18.3998", min_cov_1, max_cov_1)$p1,
  cov_plot("dat1_Tdi_filt_1000", "Tdi_F_ReSeq_Di02", min_cov_1, max_cov_1)$p1,
  cov_plot("dat1_Tdi_filt_1000", "Tdi_F_ReSeq_Di04", min_cov_1, max_cov_1)$p1,
  cov_plot("dat1_Tdi_filt_1000", "Tdi_F_ReSeq_Di06", min_cov_1, max_cov_1)$p1,
  cov_plot("dat1_Tdi_filt_1000", "Tdi_F_ReSeq_Di08", min_cov_1, max_cov_1)$p1,
  cov_plot("dat1_Tdi_filt_1000", "Tdi_F_ReSeq_Di10", min_cov_1, max_cov_1)$p1,
  ncol = 2, nrow = 5)

pdf(paste("Tdi_v8_cov_min=", min_cov_1, "_max=", max_cov_1,map_type, ".pdf", sep = ""), width = 	10, height = 15)
plot_grid(Tdi_a1, ncol = 1)
dev.off()
getwd() ## where has my plot gone....?

### sum coverage in males and females 



MF_cov_sum <- function(df, sp){
  
  female_covs <- df[,grep(paste("^",sp,"_F", sep = ""),colnames(df))]
  male_covs   <- df[,grep(paste("^",sp,"_M", sep = ""),colnames(df))]
  df$Female_cov_sum <- rowSums(female_covs)
  
  print(head(male_covs))
  
  print(length(male_covs))
  if (length(male_covs) > 100){
    df$Male_cov_sum   <- male_covs
  } else {
    df$Male_cov_sum   <- rowSums(male_covs)
  }	
  
  ##### filter contigs with 0 cov in females or males 
  
  df_filt  <- subset(df , df$Male_cov_sum   > 0 & df$Female_cov_sum > 0)
  
  print(length(df[,1]))
  print(length(df_filt[,1]))
  
  ##### contigs with 0 cov in females or males 
  
  MF_0   <- subset(df , df$Male_cov_sum   == 0 | df$Female_cov_sum == 0)
  print(length(	MF_0[,1]))
  
  ## norm by average cov
  
  male_median_cov   = median(df_filt$Male_cov_sum)
  female_median_cov = median(df_filt$Female_cov_sum)
  
  # male_mode_cov   = mlv(df_filt$Male_cov_sum, method = "hsm")
  # female_mode_cov = mlv(df_filt$Female_cov_sum, method = "hsm")
  
  male_mode_cov   = mlv(df_filt$Male_cov_sum, method = "shorth")
  female_mode_cov = mlv(df_filt$Female_cov_sum, method = "shorth")
  
  df_filt$Female_cov_sum_norm_mode = df_filt$Female_cov_sum / female_mode_cov
  df_filt$Male_cov_sum_norm_mode   = df_filt$Male_cov_sum /   male_mode_cov
  
  df_filt$Female_cov_sum_norm_med = df_filt$Female_cov_sum / female_median_cov
  df_filt$Male_cov_sum_norm_med   = df_filt$Male_cov_sum   / male_median_cov		
  
  #### calc M to F ratio
  ### Scaffolds were considered to be X candidates if they had Log2(M/F coverage) within the range [A_coord -1.1, A_coord -0.9];
  
  df_filt$M_F_mode <- log2(df_filt$Male_cov_sum_norm_mode / df_filt$Female_cov_sum_norm_mode)
  df_filt$M_F_med <- log2(df_filt$Male_cov_sum_norm_med / df_filt$Female_cov_sum_norm_med)
  
  
  print(head(df_filt))
  
  p1_M = ggplot(df_filt, aes(x=Male_cov_sum)) +
    theme_bw() +
    geom_histogram(color="darkblue", fill="blue", binwidth=0.5,alpha = 0.2) +
    coord_cartesian(xlim=c(0,200))  +
    geom_hline(yintercept = 0) +
    xlab("cov") +
    geom_vline(xintercept = male_median_cov, linetype="dotted", color = "red", size=1.5) + 
    geom_vline(xintercept = male_mode_cov, linetype="dotted", color = "green", size=1.5) + 
    ggtitle(paste(sp, " v8 total male cov, len >= ", cutoff_len, sep = ""))
  
  p1_F = ggplot(df_filt, aes(x=Female_cov_sum)) +
    theme_bw() +
    geom_histogram(color="darkblue", fill="blue", binwidth=0.5,alpha = 0.2) +
    coord_cartesian(xlim=c(0,200))  +
    geom_hline(yintercept = 0) +
    xlab("cov") +
    geom_vline(xintercept = female_median_cov, linetype="dotted", color = "red", size=1.5) + 
    geom_vline(xintercept = female_mode_cov, linetype="dotted", color = "green", size=1.5) +      	
    ggtitle(paste(sp, " v8 total female cov, len >= ", cutoff_len, sep = ""))
  
  p1 = ggplot(df_filt, aes(x=M_F_mode)) +
    theme_bw() +
    geom_histogram(color="darkblue", fill="blue", binwidth=0.01,alpha = 0.2) +
    coord_cartesian(xlim=c(-3,3))  +
    geom_hline(yintercept = 0) +
    xlab("log2(M cov / F cov)") + 
    ggtitle(paste(sp, " v8 len >= ", cutoff_len, sep = ""))
  
  p2 <- p1 + geom_rect(aes(xmin = -0.9,  xmax = -1.1,  ymin = -Inf, ymax = Inf), fill = "red", alpha = 0.01)
  p3 <- p2 + geom_rect(aes(xmin = -0.85, xmax = -1.15, ymin = -Inf, ymax = Inf), fill = "red", alpha = 0.01)
  
  #### adjusted peak grab
  df_filt_cut <- subset(df_filt, df_filt$M_F_mode < -0.5) 
  Sex_chr_peak <- peakfinder(df_filt_cut$M_F_mode)
  
  
  
  p2_adj <- p1 + geom_rect(aes(xmin = Sex_chr_peak + 0.1,  xmax = Sex_chr_peak - 0.1,  ymin = -Inf, ymax = Inf), fill = "red", alpha = 0.01) + geom_vline(xintercept = Sex_chr_peak, linetype="dotted", color = "yellow", size=0.3) 
  p3_adj <- p1 + geom_rect(aes(xmin = Sex_chr_peak + 0.15, xmax = Sex_chr_peak - 0.15, ymin = -Inf, ymax = Inf), fill = "red", alpha = 0.01) + geom_vline(xintercept = Sex_chr_peak, linetype="dotted", color = "yellow", size=0.3) 
  
  df_exp <- as.data.frame(cbind(
    as.character(df_filt$contig), 
    df_filt$Female_cov_sum,
    df_filt$Male_cov_sum,
    df_filt$Female_cov_sum_norm_mode ,
    df_filt$Male_cov_sum_norm_mode,
    df_filt$M_F_mode,
    df_filt$length
  ))	
  colnames(df_exp) <- c("contig_name", "Female_cov_sum","Male_cov_sum","Female_cov_sum_norm_mode","Male_cov_sum_norm_mode","M_F_mode","scaf_len")
  
  out_list = list("p1_M" = p1_M, "p1_F" = p1_F, "p2" = p2, "p2_adj" = p2_adj, "p3" = p3, "p3_adj" = p3_adj, "Sex_chr_peak" = Sex_chr_peak, "df_all" = df_filt, "df_exp" = df_exp)
  return(out_list)	
  
}




Tdi_out <- MF_cov_sum(dat1_Tdi_filt_1000, "Tdi")



######################################################################################################################
##### export tables and plots

### sex chromosome adjusted peak
sex_chr_peaks_df <- as.data.frame(
  cbind(c("Tdi"),
        c(
          Tdi_out$Sex_chr_peak
        )))

colnames(sex_chr_peaks_df) <- c("Sp", "sex_chr_peaks")
write.csv(sex_chr_peaks_df, file=paste("sex_chr_peaks_", cutoff_len, map_type, ".csv", sep = ""), row.names=FALSE)



###### full tables

head(Tdi_out$df_exp)
write.csv(Tdi_out$df_exp, file=paste("Tdi_v8_MFcov_filt_", cutoff_len, map_type, ".csv", sep = ""), row.names=FALSE)

####### plots
### with adjusted peak and summed cov 

# Tdi
pdf(paste("Tdi_v8_MFcov_filt_", cutoff_len, map_type, "_wAdjpeak2.pdf", sep = ""), width = 5, height = 5)
Tdi_out$p2_adj
dev.off()
png(filename = paste("Tdi_v8_MFcov_filt_", cutoff_len, map_type, "_wAdjpeak2.png", sep = ""), width = 5, height = 5, units = "in", bg = "white", res = 300)
Tdi_out$p2_adj
dev.off()

pdf(paste("Tdi_v8_MFcov_filt_", cutoff_len,  map_type,"_wAdjpeak3.pdf", sep = ""), width = 5, height = 5)
Tdi_out$p3_adj
dev.off()
png(filename = paste("Tdi_v8_MFcov_filt_", cutoff_len,  map_type,"_wAdjpeak3.png", sep = ""), width = 5, height = 5, units = "in", bg = "white", res = 300)
Tdi_out$p3_adj
dev.off()

pdf(paste("Tdi_v8_MFsumcov_filt_", cutoff_len,  map_type,".pdf", sep = ""), width = 10, height = 10)
plot_grid(Tdi_out$p1_M, Tdi_out$p1_F, ncol = 1, nrow = 2)
dev.off()
png(filename = paste("Tdi_v8_MFsumcov_filt_", cutoff_len, map_type, ".png", sep = ""), width = 10, height = 10, units = "in", bg = "white", res = 300)
plot_grid(Tdi_out$p1_M, Tdi_out$p1_F, ncol = 1, nrow = 2)
dev.off()







######################### peak adjust all vals 

head(Tdi_out$df_all)

Tdi_out_df_all <- Tdi_out$df_all
Tdi_out_df_all$M_F_mode_padj <- Tdi_out_df_all$M_F_mode +  (-1-Tdi_out$Sex_chr_peak)

head(Tdi_out_df_all)



p1 <- ggplot(Tdi_out_df_all, aes(x=M_F_mode_padj)) +
  theme_bw() +
  geom_histogram(color="darkblue", fill="blue", binwidth=0.01,alpha = 0.2) +
  coord_cartesian(xlim=c(-3,3))  +
  geom_hline(yintercept = 0) +
  xlab("log2(M cov / F cov)") + 
  ggtitle(paste("Tdi peak adj all", " v8 len >= ", cutoff_len, sep = ""))



png(filename = paste("Tdi_v8_MFcov_filt_", cutoff_len, map_type, "_all_peak_adjusted.png", sep = ""), width = 5, height = 5, units = "in", bg = "white", res = 300)
p1
dev.off()









########################################################################################################################################################################
####### output session info
print (sessionInfo())
writeLines(capture.output(sessionInfo()), paste("male_tim_cov_v8.R_sessionInfo_len_", cutoff_len, map_type,".txt"), sep = "")