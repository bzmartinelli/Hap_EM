library(ggplot2)
library("RColorBrewer")


SqEall = data.frame()
all = data.frame()
mean_diffs = data.frame()
pvalue = c()
diffmean_and_pvalues = data.frame()
positions = c()
all_reads_haps = data.frame()
Ndepth = c()

for(depth in seq(10,60,10)) {
  for(position in 0:4) {
    Ndepth = c(Ndepth, depth)
    positions = c(positions, position)
    filenameS = paste("Simulated_Meth_Proportion_depth",depth,sep="")
    filenameE = paste("Estimated_Meth_Proportion_depth",depth,sep="")
    dataS = read.table(filenameS,h=T, sep="\t")
    dataE = read.table(filenameE,h=T, sep="\t")
    dataE = dataE[dataE$Position==position,]
    dataS = dataS[dataS$Position==position,]
    df_reads = data.frame(Simulated_Frequency=dataS$Proportion, 
                          Estimated_Frequency=dataE$From_reads, 
                          Method=factor(c(rep("Reads"))), Sim_ID=dataS$Sim_id)
    df_haps = data.frame(Simulated_Frequency=dataS$Proportion, 
                         Estimated_Frequency=dataE$From_haplotypes, 
                         Method=factor(c(rep("Haplotypes"))), Sim_ID=dataS$Sim_id)
    df = rbind(df_reads, df_haps)
    df_hr = data.frame(dataE)
    df_hr = cbind(df_hr, Simulated = dataS$Proportion )
    all_reads_haps = rbind(all_reads_haps, df_hr)    
    df$Position = factor(c(rep(position,nrow(df))))
    
    Squared_error = (df$Estimated_Frequency-df$Simulated_Frequency)**2
    df = cbind(df, Squared_error)
    df$Depth = factor(c(rep(depth,nrow(df))))
    all = rbind(all,df)    
    SE_reads = (sd(df[df$Method=="Reads",6]))/sqrt(as.numeric(length(df$Method=="Reads")))
    
    SqE = data.frame(MeanSE=mean(df[df$Method=="Reads",6]), 
                     SD=sd(df[df$Method=="Reads",6]), SE=SE_reads,
                     Method=factor(c(rep("Reads"))))    
    
    SE_haps = (sd(df[df$Method=="Haplotypes",6]))/sqrt(as.numeric(length(df$Method=="Haplotypes")))    
    SqE = rbind(SqE,data.frame(MeanSE=mean(df[df$Method=="Haplotypes",6]), 
                               SD=sd(df[df$Method=="Haplotypes",6]), SE=SE_haps, 
                               Method=factor(c(rep("Haplotypes")))))    
    SqE$Position = factor(c(rep(position,nrow(SqE))))
    SqE$Depth = factor(c(rep(depth,nrow(SqE))))
    SqEall = rbind(SqEall,SqE)
    
    mean_diffs = data.frame(Meandiff=mean(df[df$Method=="Reads",6] - mean(df[df$Method=="Haplotypes",6])),
                            Position=df$Position)
    mean_diffs$pvalues = c(pvalue,t.test(df[df$Method=="Reads",6], df[df$Method=="Haplotypes",6])$p.value)
    diffmean_and_pvalues = rbind(diffmean_and_pvalues, mean_diffs)
  }
}



new_labels = c('10'="depth 10", '20'="depth 20", '30'="depth 30",
               '40'="depth 40", '50'="depth 50", '60'="depth 60")

#plots
jpeg("meth_proport_plots.jpg", res=300, height=6, width=10, units="in")
ggplot(SqEall, aes(x=as.factor(Position), y=MeanSE, group=Method, colour=Method)) +
  geom_pointrange(aes(ymin=MeanSE-SE, ymax=MeanSE+SE), fatten=1, 
                  position=position_dodge(width=0.4)) + scale_x_discrete(labels=1:5) +
  # scale_color_brewer(palette="Paired",name="Method")+
  scale_color_manual(name="Method",values = c(brewer.pal(3, "Blues")[3:3], brewer.pal(4, "Set1")[4:4])) +
  
  labs( title="Methylation proportions",x="CpG position within the haplotypes", y = "Squared error") + theme_grey()+
  theme(legend.position = "bottom", plot.title=element_text(hjust=0.5, face="bold",size=15),axis.title=element_text(size=12), 
        axis.text.y = element_text(size=10),
        axis.text.x = element_text(size=10))  +
  facet_wrap(~Depth, scales = "free",  ncol=3, labeller = as_labeller(new_labels))
dev.off()


