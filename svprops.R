library(ggplot2)
library(scales)
library(grid)
library(reshape2)

args = commandArgs(trailingOnly=TRUE)
x = read.table(args[1], header=T)

# Plotting SV types
#plotSVTypes=factor(c("DEL", "DUP", "DRDEL", "DIDEL", "INVDUP", "PDUP", "INVDEL", "INV"))
plotSVTypes=factor(c("DEL", "DUP"))

# Keep only desired SV types
x=x[x$svType %in% plotSVTypes,]
x$svType=factor(x$svType, levels=plotSVTypes)

# VAF binning
x$vafbin=cut(x$af, breaks=c(-1, 0.001, 0.01, 0.1, 1), labels=c("<0.001", "0.001-0.01", "0.01-0.1", ">0.1"))

# Theme
lSize=1.2
txtFontSize=10
axisFontSize=16
axisTtlFontSize=22
lgdTtlFontSize=22
lgdFontSize=16
scienceTheme=theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank(), legend.key=element_blank(), legend.background=element_blank(), panel.background = element_blank(), panel.border=element_blank(), strip.background = element_blank(), axis.line=element_line(size=0.7, color="black"), axis.text.x=element_text(size=axisFontSize), axis.text.y=element_text(size=axisFontSize), axis.title.x=element_text(size=axisTtlFontSize), axis.title.y=element_text(size=axisTtlFontSize), legend.title=element_text(size=lgdTtlFontSize, face="bold"), legend.text=element_text(size=lgdFontSize), text=element_text(size=txtFontSize))

png("svprops.png", height=800, width=1200)
grid.newpage()
pushViewport(viewport(layout=grid.layout(2,5)))

# Binned variant allele frequency
freqX=melt(table(x$vafbin, x[,c("svType")]), varnames=c("vafbin", "svType"))
p1=ggplot(data=freqX, aes(x=vafbin, y=value)) + geom_line(aes(group=svType, colour=svType), size=lSize, show_guide=FALSE) + geom_point(aes(colour=svType), size=lSize, show_guide=FALSE)
p1=p1 + xlab("Variant allele frequency") + ylab("#SV sites")
p1=p1 + scale_y_log10(breaks=c(1,10,100,1000,10000), limits=c(1,max(freqX$value)))
p1=p1 + scienceTheme + labs(colour="SV Type")

# Variant allele count frequency spectrum
freqX=melt(table(x$vac, x[,c("svType")]), varnames=c("vac", "svType"))
p2=ggplot(data=freqX, aes(x=vac, y=value)) + geom_line(aes(group=svType, colour=svType), size=lSize)
p2=p2 + xlab("Variant allele count") + ylab("#SV sites") 
p2=p2 + scale_x_log10(breaks=c(1,10,100,1000)) + scale_y_log10(breaks=c(1,10,100,1000), limits=c(1,max(freqX$value)))
p2=p2 + scienceTheme + labs(colour="SV Type") + theme(legend.position=c(1,1), legend.justification=c(1,1))

# Size vs. VAF
sizeVaf=melt(tapply(x$size, x[,c("vafbin", "svType")], median), varnames=c("vafbin", "svType"))
p3=ggplot(data=sizeVaf, aes(x=vafbin, y=value))  + geom_line(aes(group=svType, colour=svType), size=lSize, show_guide=FALSE) + geom_point(aes(colour=svType), size=lSize, show_guide=FALSE)
p3=p3 + xlab("Variant allele frequency") + ylab("Median SV Size") + ylim(0, max(sizeVaf$value))
p3=p3 + scienceTheme + labs(colour="SV Type")

# Size histograms
p4=ggplot(x, aes(x=size)) + geom_freqpoly(aes(group=svType, colour=svType), binwidth=0.1, size=lSize, show_guide=FALSE) 
p4=p4 + xlab('SV size') + ylab('#SV sites') + scale_y_log10(breaks=c(1,10,100,1000))
p4=p4 + scale_x_log10(limits=c(50, max(x$size)), expand=c(0,0), breaks=c(1e2, 1e3, 1e4, 1e5, 1e6), labels=c('100bp', '1kb', '10kb', '100kb', '1Mb'))
p4=p4 + scienceTheme + labs(colour="SV Type")

print(p2, vp = viewport(layout.pos.row=1, layout.pos.col=1:3))
print(p4, vp = viewport(layout.pos.row=2, layout.pos.col=1:3))
print(p1, vp = viewport(layout.pos.row=1, layout.pos.col=4:5))
print(p3, vp = viewport(layout.pos.row=2, layout.pos.col=4:5))
dev.off()
print(warnings())
