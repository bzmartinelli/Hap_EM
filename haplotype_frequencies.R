
library(ggplot2)
library("RColorBrewer")
library(cowplot)

# Normalized frequencies

Props = data.frame()
Haps_Freqs = data.frame()
Ndepth = c()
for(depth in seq(10,60,10)) {
  Ndepth = c(Ndepth, depth)
  filename = paste("Frequencies_depth",depth,"_sorted_normalized", sep="")
  data = read.table(filename,h=T,sep="\t", colClasses="character", as.is=T)
  df = data.frame(Haplotype=data$Haplotype, 
                  Simulated_Frequency=as.numeric(data$Simulated_Frequency),
                  Estimated_Frequency=as.numeric(data$Estimated_Frequency), 
                  ID = factor(data$Sim_id)) 
  df$Depth = factor(c(rep(depth,nrow(df))))
  Haps_Freqs = rbind(Haps_Freqs, df)
  
  
  p1 = as.numeric(length(which(df$Estimated_Frequency>0 & df$Simulated_Frequency == 0.1))/1000)
  p2 = as.numeric(length(which(df$Estimated_Frequency>0 & df$Simulated_Frequency == 0.2))/1000)
  p3 = as.numeric(length(which(df$Estimated_Frequency>0 & df$Simulated_Frequency == 0.3))/1000)
  p4 = as.numeric(length(which(df$Estimated_Frequency>0 & df$Simulated_Frequency == 0.4))/1000)
  
  SEp1 = sqrt(p1 *(1-p1)/1000)
  SEp2 = sqrt(p2 *(1-p2)/1000)
  SEp3 = sqrt(p3 *(1-p3)/1000)
  SEp4 = sqrt(p4 *(1-p4)/1000)
  
  propdf = data.frame(Simulated_Frequency=as.factor(0.1),
                      Proportion_Haps=as.numeric(p1), SE=SEp1, Depth=as.factor(depth))
  propdf = rbind(propdf, data.frame(Simulated_Frequency=as.factor(0.2),
                                    Proportion_Haps=as.numeric(p2), SE=SEp2, Depth=as.factor(depth)))
  propdf = rbind(propdf, data.frame(Simulated_Frequency=as.factor(0.3),
                                    Proportion_Haps=as.numeric(p3), SE=SEp3, Depth=as.factor(depth)))
  propdf = rbind(propdf, data.frame(Simulated_Frequency=as.factor(0.4),
                                    Proportion_Haps=as.numeric(p4), SE=SEp4, Depth=as.factor(depth)))
  Props = rbind(Props, propdf)
}

p1 = ggplot(Props,aes(x=Simulated_Frequency, y=Proportion_Haps, color=Depth)) +
  geom_point(position=position_dodge(width=0.6)) +
  geom_errorbar(aes(ymax = Proportion_Haps+SE, ymin=Proportion_Haps-SE),
                position=position_dodge(width=0.6) ) +
  scale_y_continuous( limits=c(0, 1)) +
  # theme_gray() +
  scale_color_manual(values = brewer.pal(9, "Blues")[4:9]) +
  labs(title="Recovered Haplotypes",
       x="Simulated Frequency", y="Proportion of recovered haplotypes",
       color="Depth") +
  theme(plot.title=element_text(hjust=0.5,face="bold"),
        legend.position="bottom",
        legend.title=element_text(size=8) , legend.text=element_text(size=7)) 


p2 = ggplot(Haps_Freqs, aes(x=as.factor(Simulated_Frequency), 
                            y=Estimated_Frequency, fill=Depth)) + 
  geom_boxplot(notch=TRUE, outlier.size = 0.3, lwd=0.35) + 
  scale_y_continuous(breaks=seq(0, 1, by=0.1)) +
  scale_fill_manual(values = brewer.pal(9, "Blues")[4:9]) +
  
  # theme_gray() +
  labs(title="Normalized Haplotype Frequencies", x="Simulated Frequency", 
       y="Estimated Frequency", fill = "Depth")+
  theme(plot.title=element_text(hjust=0.5, face="bold"), 
        legend.position="bottom",
        legend.title=element_text(size=8) , legend.text=element_text(size=7))



# multiplot 
plot_grid(p1, p2, ncol=2, labels=c('A', 'B'))
jpeg("hap_freq_prop_normalized.jpg", res=300,  height=5, width=10, units="in")
plot_grid(p1, p2, nrow=1, labels=c('A', 'B'))
dev.off()



# Haplotype frequencis (all)

Props = data.frame()
Haps_Freqs = data.frame()
Ndepth = c()
for(depth in seq(10,60,10)) {
  Ndepth = c(Ndepth, depth)
  filename = paste("Frequencies_depth",depth,"_sorted", sep="")
  data = read.table(filename,h=T,sep="\t", colClasses="character", as.is=T)
  df = data.frame(Haplotype=data$Haplotype, 
                  Simulated_Frequency=as.numeric(data$Simulated_Frequency),
                  Estimated_Frequency=as.numeric(data$Estimated_Frequency), 
                  ID = factor(data$Sim_id)) 
  df$Depth = factor(c(rep(depth,nrow(df))))
  Haps_Freqs = rbind(Haps_Freqs, df)
}


jpeg("HapsFreqs.jpg", res=350, height=6, width=9, units="in")
ggplot(Haps_Freqs, aes(x=as.factor(Simulated_Frequency), 
                       y=Estimated_Frequency, fill=Depth)) + 
  geom_boxplot(notch=TRUE, outlier.size = 0.3, lwd=0.35) + 
  scale_y_continuous(breaks=seq(0, 1, by=0.1)) +
  scale_fill_manual(values = brewer.pal(9, "Blues")[4:9]) +
  
  theme_gray() +
  labs(title="Haplotype Frequencies", x="Simulated Frequency", 
       y="Estimated Frequency", fill = "Depth")+
  theme(plot.title=element_text(hjust=0.5, face="bold"), 
        legend.position="bottom",
        legend.title=element_text(size=8) , legend.text=element_text(size=7))
dev.off()
