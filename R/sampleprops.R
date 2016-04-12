library(ggplot2)
library(scales)
library(grid)
library(RColorBrewer)
library(reshape2)

args=commandArgs(trailingOnly=TRUE)

x = read.table(args[1], header=TRUE)
x = x[,c("sample","missing","het","homalt")]
counts = melt(x, id.vars=c("sample"))

txtFontSize=16
axisFontSize=16
axisTtlFontSize=22
lgdTtlFontSize=22
lgdFontSize=16
scienceTheme=theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank(), legend.key=element_blank(), legend.background=element_blank(), panel.background = element_blank(), panel.border=element_blank(), strip.background = element_blank(), axis.line.x=element_line(size=0.7, color="black"), axis.line.y=element_line(size=0.7, color="black"), axis.text.x=element_text(size=axisFontSize), axis.text.y=element_text(size=axisFontSize), axis.title.x=element_text(size=axisTtlFontSize), axis.title.y=element_text(size=axisTtlFontSize), legend.title=element_text(size=lgdTtlFontSize, face="bold"), legend.text=element_text(size=lgdFontSize), text=element_text(size=txtFontSize))


png("sampleprops.png", height=800, width=1200)
p=ggplot(data=counts, aes(x=variable, y=value)) + geom_boxplot() 
p=p + scienceTheme + xlab("Genotype") + ylab("Count")
p
dev.off()
print(warnings())
