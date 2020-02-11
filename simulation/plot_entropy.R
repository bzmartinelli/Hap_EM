
library(ggplot2)
library(gridExtra)


# Entropy

Entropy = data.frame()
Ndepth = c()
for(depth in seq(10,60,10)) {
    Ndepth = c(Ndepth, depth)
    filenameS = paste("Simulated_Entropy_5cpgs_depth", depth, sep="")
    filenameE = paste("Estimated_Entropy_5cpgs_depth", depth, sep="")
    dataS = read.table(filenameS,h=T, sep="\t")
    dataE = read.table(filenameE,h=T, sep="\t")
    df = data.frame(Estimated_Entropy=dataE$Estimated_Entropy, 
                    Simulated_Entropy=as.numeric(dataS$Simulated_entropy),
                    ID = factor(dataE$Sim_id))
    df$Depth = factor(c(rep(depth,nrow(df))))
    Entropy = rbind(Entropy,df)
}

new_labels = c('10'="depth 10", '20'="depth 20", '30'="depth 30",
               '40'="depth 40", '50'="depth 50", '60'="depth 60")

# make plots and save as image
jpeg("entropy_plots.jpg", res=350, height=6, width=9, units="in")
ggplot(Entropy, aes(x=Estimated_Entropy)) + geom_density() +
  geom_vline(data=Entropy, aes(xintercept=Entropy$Simulated_Entropy, 
                               linetype="simulated entropy"), color="red") +
  scale_linetype_manual(name = "", values = c(1, 1),
                        guide = guide_legend(override.aes = list(color = c("red")))) + theme_grey()+
  labs(title="Entropy", x="Estimated entropy", y="Density") + 
  theme(legend.position = "bottom", plot.title=element_text(hjust=0.5, face="bold",size=15), axis.title=element_text(size=12))  + 
  facet_wrap(~Depth, labeller = as_labeller(new_labels)) 
dev.off()





# Normalized Entropy

Entropy = data.frame()
Ndepth = c()
for(depth in seq(10,60,10)) {
  Ndepth = c(Ndepth, depth)
  filenameS = paste("Simulated_Entropy_5cpgs_depth",depth,sep="")
  filenameE = paste("Estimated_Entropy_normalized_depth",depth,sep="") 
  dataS = read.table(filenameS,h=T, sep="\t")
  dataE = read.table(filenameE,h=T, sep="\t") 
  df = data.frame(Estimated_Entropy=dataE$Estimated_entropy, 
                  Simulated_Entropy=as.numeric(dataS$Simulated_entropy),
                  ID = factor(dataE$Sim_id))
  df$Depth = factor(c(rep(depth,nrow(df))))
  Entropy = rbind(Entropy,df)
}

new_labels = c('10'="depth 10", '20'="depth 20", '30'="depth 30",
               '40'="depth 40", '50'="depth 50", '60'="depth 60")

# make plots and save as image
jpeg("entropy_normalized_plots.jpg", res=300, height=6, width=10, units="in")
ggplot(Entropy, aes(x=Estimated_Entropy)) + geom_density() +
  geom_vline(data=Entropy, aes(xintercept=Entropy$Simulated_Entropy, 
                               linetype="simulated entropy"), color="red") +
  scale_linetype_manual(name = "", values = c(1, 1),
                        guide = guide_legend(override.aes = list(color = c("red")))) + theme_grey()+
  labs(title="Normalized Entropy", x="Estimated entropy", y="Density") + 
  theme(legend.position = "bottom", plot.title=element_text(hjust=0.5, face="bold",size=15), axis.title=element_text(size=12))  + 
  facet_wrap(~Depth, labeller = as_labeller(new_labels), ncol = 3) 
dev.off()


