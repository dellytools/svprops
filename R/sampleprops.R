library(ggplot2)
library(scales)
library(grid)
library(RColorBrewer)
library(reshape2)

args=commandArgs(trailingOnly=TRUE)

x = read.table(args[1], header=TRUE)
x = x[order(x$superpopulation),]
x$population = factor(x$population, levels=unique(x$population))

txtFontSize=16
axisFontSize=16
axisTtlFontSize=22
lgdTtlFontSize=22
lgdFontSize=16
scienceTheme=theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank(), legend.key=element_blank(), legend.background=element_blank(), panel.background = element_blank(), panel.border=element_blank(), strip.background = element_blank(), axis.line=element_line(size=0.7, color="black"), axis.text.x=element_text(size=axisFontSize), axis.text.y=element_text(size=axisFontSize), axis.title.x=element_text(size=axisTtlFontSize), axis.title.y=element_text(size=axisTtlFontSize), legend.title=element_text(size=lgdTtlFontSize, face="bold"), legend.text=element_text(size=lgdFontSize), text=element_text(size=txtFontSize))


png("sampleprops.png", height=800, width=1200)
grid.newpage()
pushViewport(viewport(layout=grid.layout(2,2)))

# Missing
p1=ggplot(data=x, aes(x=population, y=missing)) + geom_violin(aes(fill=superpopulation), adjust=2) 
p1=p1 + geom_point(position=position_jitter(w=0.05, h=0.05), alpha=0.1, size=0.1)
p1=p1 + ylab("#Missing Genotypes") + xlab("Populations") + coord_flip() 
p1=p1 + scienceTheme

# het
p2=ggplot(data=x, aes(x=population, y=het)) + geom_violin(aes(fill=superpopulation), adjust=2, show_guide=FALSE) 
p2=p2 + geom_point(position=position_jitter(w=0.05, h=0.05), alpha=0.1, size=0.1)
p2=p2 + ylab("#Heterozygous Genotypes") + xlab("Populations") + coord_flip() 
p2=p2 + scienceTheme

# hom
p3=ggplot(data=x, aes(x=population, y=homalt)) + geom_violin(aes(fill=superpopulation), adjust=2, show_guide=FALSE) 
p3=p3 + geom_point(position=position_jitter(w=0.05, h=0.05), alpha=0.1, size=0.1)
p3=p3 + ylab("#Homozygous Genotypes") + xlab("Populations") + coord_flip() 
p3=p3 + scienceTheme

# hetX
p4=ggplot(data=x, aes(x=gender, y=hetX)) + geom_boxplot(aes(group=gender)) 
p4=p4 + geom_point(aes(color=superpopulation), position=position_jitter(w=0.05, h=0.05), alpha=0.1, size=0.1)
p4=p4 + ylab("#Het. chrX") + xlab("Gender") + coord_flip() 
p4=p4 + scienceTheme


print(p1, vp = viewport(layout.pos.row=1, layout.pos.col=1))
print(p2, vp = viewport(layout.pos.row=1, layout.pos.col=2))
print(p4, vp = viewport(layout.pos.row=2, layout.pos.col=1))
print(p3, vp = viewport(layout.pos.row=2, layout.pos.col=2))
dev.off()
print(warnings())
