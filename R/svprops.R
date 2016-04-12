library(ggplot2)
library(scales)
library(grid)
library(reshape2)

args = commandArgs(trailingOnly=TRUE)
x = read.table(args[1], header=T)

if (!("svtype" %in% colnames(x))) { x$svtype = "Variant"; }

# Plotting SV types
#plotsvt = c("DEL", "INS", "DUP")
plotsvt = unique(x$svtype)
plotsvt = plotsvt[plotsvt %in% unique(x$svtype)]
plotSVTypes=factor(plotsvt)

# Keep only desired SV types
x=x[x$svtype %in% plotSVTypes,]
x$svtype=factor(x$svtype, levels=plotSVTypes)

# VAF binning
if (max(x$vac)<100) {
  x$vafbin=cut(x$vaf, breaks=c(-1, 0.01, 0.1, 1), labels=c("<0.01", "0.01-0.1", ">0.1"))
} else {
  x$vafbin=cut(x$vaf, breaks=c(-1, 0.001, 0.01, 0.1, 1), labels=c("<0.001", "0.001-0.01", "0.01-0.1", ">0.1"))
}

# Break vector
brVector=c(1,10,100,1000,10000,100000,1000000,10000000)

# Theme
lSize=1.2
txtFontSize=10
axisFontSize=16
axisTtlFontSize=22
lgdTtlFontSize=22
lgdFontSize=16
scienceTheme=theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank(), legend.key=element_blank(), legend.background=element_blank(), panel.background = element_blank(), panel.border=element_blank(), strip.background = element_blank(), axis.line.x=element_line(size=0.7, color="black"), axis.line.y=element_line(size=0.7, color="black"), axis.text.x=element_text(size=axisFontSize), axis.text.y=element_text(size=axisFontSize), axis.title.x=element_text(size=axisTtlFontSize), axis.title.y=element_text(size=axisTtlFontSize), legend.title=element_text(size=lgdTtlFontSize, face="bold"), legend.text=element_text(size=lgdFontSize), text=element_text(size=txtFontSize))




# Variant allele frequency spectrum
png("vaf.png", height=800, width=1200)
freqX=melt(table(x$vafbin, x[,c("svtype")]), varnames=c("vafbin", "svtype"))
p1=ggplot(data=freqX, aes(x=vafbin, y=value)) + geom_line(aes(group=svtype, colour=svtype), size=lSize, show.legend=FALSE) + geom_point(aes(colour=svtype), size=lSize, show.legend=FALSE)
p1=p1 + xlab("Variant allele frequency") + ylab("#Variants")
p1=p1 + scale_y_log10(breaks=brVector[max(freqX$value)>brVector/10], limits=c(1,max(freqX$value)))
p1=p1 + scienceTheme + labs(colour="SV Type")
p1
z=dev.off()
print(warnings())


# Variant allele count spectrum
png("vac.png", height=800, width=1200)
freqX=melt(table(x$vac, x[,c("svtype")]), varnames=c("vac", "svtype"))
p2=ggplot(data=freqX, aes(x=vac, y=value)) + geom_line(aes(group=svtype, colour=svtype), size=lSize)
p2=p2 + xlab("Variant allele count") + ylab("#Variants") 
p2=p2 + scale_y_log10(breaks=brVector[max(freqX$value)>brVector/10], limits=c(1,max(freqX$value)))
if (max(x$vac)<100) {
   p2=p2 + scale_x_continuous()
} else {
   p2=p2 + scale_x_log10(breaks=brVector[max(freqX$vac)>brVector/10])
}
p2=p2 + scienceTheme + labs(colour="SV Type") + theme(legend.position=c(1,1), legend.justification=c(1,1))
p2
z=dev.off()
print(warnings())

# Size statistics
if (length(unique(x$size)) > 100) {
   # Size vs. VAF
   png("sizevaf.png", height=800, width=1200)
   sizeVaf=melt(tapply(x$size, x[,c("vafbin", "svtype")], median), varnames=c("vafbin", "svtype"))
   p3=ggplot(data=sizeVaf, aes(x=vafbin, y=value))  + geom_line(aes(group=svtype, colour=svtype), size=lSize, show.legend=FALSE) + geom_point(aes(colour=svtype), size=lSize, show.legend=FALSE)
   p3=p3 + xlab("Variant allele frequency") + ylab("Median SV Size") + scale_y_log10(breaks=c(1,10,100,1000,10000,100000))
   p3=p3 + scienceTheme + labs(colour="SV Type")
   p3
   z=dev.off()
   print(warnings())

   # Size histograms
   png("size.png", height=800, width=1200)
   p4=ggplot(x, aes(x=size)) + geom_freqpoly(aes(group=svtype, color=svtype), binwidth=0.1, size=lSize, show.legend=FALSE) 
   p4=p4 + xlab('SV size') + ylab('#SV sites') + scale_y_log10(breaks=c(1,10,100,1000))
   p4=p4 + scale_x_log10(limits=c(50, max(x$size)), expand=c(0,0), breaks=c(1e2, 1e3, 1e4, 1e5, 1e6), labels=c('100bp', '1kb', '10kb', '100kb', '1Mb'))
   p4=p4 + scienceTheme + labs(colour="SV Type")
   p4
   z=dev.off()
   print(warnings())
}

# Plot missing rate
png("missingrate.png", height=800, width=1200)
p5=ggplot(data=x, aes(svtype, missingrate)) + geom_boxplot()
p5=p5 + xlab("SV Type") + ylab("Missing Genotype Rate")
p5=p5 + scienceTheme
p5
z=dev.off()
print(warnings())

# Plot singleton counts
png("singleton.png", height=800, width=1200)
freqX=melt(table(x[!is.na(x$singleton),]$singleton))
colnames(freqX) = c("sample", "count")
freqX=freqX[with(freqX, order(-freqX$count)), ]
freqX=freqX[1:50,]
freqX$sample = factor(freqX$sample, levels=freqX$sample)
p6=ggplot(data=freqX, aes(x=sample, y=count)) + geom_point()
p6=p6 + xlab("Sample") + ylab("#Singletons")
p6=p6 + scienceTheme + coord_flip()
p6
z=dev.off()
print(warnings())


