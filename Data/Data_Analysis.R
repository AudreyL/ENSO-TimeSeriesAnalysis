# Figure
m.data <- read.table('/home/audrey/Documents/Projects/Velayoudoum/2015-Enso-gegenbauer-models-ClimaticDynamic/EnsoData2015/monthly-SST-1950-2014.txt',header=TRUE)
tt <- NULL
for (t in seq(1,dim(m.data)[1])){
    tt <- c(tt,paste('28-',m.data[t,2],'-',m.data[t,1],sep=''))
}
tt <- as.Date(as.character(tt),"%d-%m-%Y")
m.data <- cbind(tt,m.data[,3:dim(m.data)[2]])

#png('/home/audrey/Documents/Projects/Velayoudoum/2015-Enso-gegenbauer-models-ClimaticDynamic/Submission2/Lustig_et_al_manuscript_June2015/ts_properties2.png',width=12, height=8, units="in", res=300)
pdf('/home/audrey/Documents/Projects/Velayoudoum/2015-Enso-gegenbauer-models-ClimaticDynamic/Submission2/Lustig_et_al_manuscript_June2015/Lustig_et_al_ts_properties2.pdf',width=12,height=8,paper='a4r')
op <- par(mfrow = c(4,4), oma = c(0,0,1,0) + 0.1, mar = c(4,4,1,0) + 0.1)
plot(m.data[,1],m.data[,2], type='l', xlab="", ylab="Nino 1-2", xlim=c(m.data[1,1],m.data[length(m.data[,1]),1]), ylim=c(min(m.data[,2]),max(m.data[,2])), main="Series")
tt<- acf(m.data[,2],lag=100,plot=FALSE)
plot(tt, xlab="",ylab="", main=" ")
title(outer=outer,main="ACF",line=0.2)
d <- density(m.data[,2])
plot(d, xlab="", ylab="", main="Density")
polygon(d, col="gray", border="black")
spectrum(m.data[,2],method="ar", main="Spectral density", xlab="",ylab="")

plot(m.data[,1],m.data[,4], type='l', xlab="", ylab="Nino 3", xlim=c(m.data[1,1],m.data[length(m.data[,1]),1]), ylim=c(min(m.data[,4]),max(m.data[,4])))
acf(m.data[,4], xlab="", ylab="", lag=100,main="")
d <- density(m.data[,4])
plot(d, xlab="", ylab="", main="")
polygon(d, col="gray", border="black")
spectrum(m.data[,4],method="ar", main=" ", xlab="",ylab="")

plot(m.data[,1],m.data[,6], type='l', xlab="", ylab="Nino 4", xlim=c(m.data[1,1],m.data[length(m.data[,1]),1]), ylim=c(min(m.data[,6]),max(m.data[,6])))
acf(m.data[,6], xlab="", ylab="", lag=100,main="")
d <- density(m.data[,6])
plot(d, xlab="", ylab="", main="")
polygon(d, col="gray", border="black")
spectrum(m.data[,6],method="ar", main=" ", xlab="",ylab="")

plot(m.data[,1],m.data[,8], type='l', xlab="Date", ylab="Nino 3.4", xlim=c(m.data[1,1],m.data[length(m.data[,1]),1]), ylim=c(min(m.data[,8]),max(m.data[,8])))
acf(m.data[,8], xlab="Time-lag", ylab="", lag=100,main="")
d <- density(m.data[,8])
plot(d, xlab="Temperature range", ylab="", main="")
polygon(d, col="gray", border="black")
spectrum(m.data[,8],method="ar", main=" ", xlab="Frequency",ylab="")
par(op)
dev.off()


# Figure 2 - time series agin itself at lag -1
#png('/home/audrey/Documents/Projects/Velayoudoum/2015-Enso-gegenbauer-models-ClimaticDynamic/Submission2/Lustig_et_al_manuscript_June2015/Lustig_et_al_scatterplot-lag1.png',width=12, height=8, units="in", res=300)
pdf('/home/audrey/Documents/Projects/Velayoudoum/2015-Enso-gegenbauer-models-ClimaticDynamic/Submission2/Lustig_et_al_manuscript_June2015/Lustig_et_al_scatterplot-lag1.pdf',width=12,height=8,paper='a4r')
par(mfrow=c(2,2))
plot(m.data[seq(1,(length(m.data[,2])-1)),2],m.data[seq(2,(length(m.data[,2]))),2], xlab='Nino 1-2', ylab='Nino 1-2 (Lag-1)', pch=20, col=colors()[184],cex.axis=1.5,cex.lab=1.5)
abline(0,1, col='blue')
plot(m.data[seq(1,(length(m.data[,2])-1)),4],m.data[seq(2,(length(m.data[,2]))),4], xlab='Nino 3', ylab='Nino 3 (Lag-1)', pch=20, col=colors()[184],cex.axis=1.5,cex.lab=1.5)
abline(0,1, col='blue')
plot(m.data[seq(1,(length(m.data[,2])-1)),6],m.data[seq(2,(length(m.data[,2]))),6], xlab='Nino 4', ylab='Nino 4 (Lag-1)', pch=20, col=colors()[184],cex.axis=1.5,cex.lab=1.5)
abline(0,1, col='blue')
plot(m.data[seq(1,(length(m.data[,2])-1)),8],m.data[seq(2,(length(m.data[,2]))),8], xlab='Nino 3.4', ylab='Nino 3.4 (Lag-1)', pch=20, col=colors()[184],cex.axis=1.5,cex.lab=1.5)
abline(0,1, col='blue')
dev.off()







######################### llplot

require('TSA')
spec.val.enso12 <- spectrum(m.data[,2], main="Spectral density")
perio.val.enso12 <- periodogram(m.data[,2],log='no',plot=TRUE,ylab="Periodogram", xlab="Frequency",lwd=2)

spec.val.enso3 <- spectrum(m.data[,4], main="Spectral density")
perio.val.enso3 <- periodogram(m.data[,4],log='no',plot=TRUE,ylab="Periodogram", xlab="Frequency",lwd=2)

spec.val.enso4 <- spectrum(m.data[,6], main="Spectral density")
perio.val.enso4 <- periodogram(m.data[,6],log='no',plot=TRUE,ylab="Periodogram", xlab="Frequency",lwd=2)

spec.val.enso34 <- spectrum(m.data[,8], main="Spectral density")
perio.val.enso34 <- periodogram(m.data[,8],log='no',plot=TRUE,ylab="Periodogram", xlab="Frequency",lwd=2)


require('longmemo')


pdf('/home/audrey/Documents/Projects/Velayoudoum/2015-Enso-gegenbauer-models-ClimaticDynamic/Submission2/Lustig_et_al_manuscript_June2015/Lustig_et_al_llplot.pdf',width=12,height=8,paper='a4r')
par(mfrow=c(2,2))
llplot(perio.val.enso12$spec, spec.val.enso12$spec)
title('Nino 1-2')
llplot(perio.val.enso3$spec, spec.val.enso3$spec)
title('Nino 3')
llplot(perio.val.enso4$spec, spec.val.enso4$spec)
title('Nino 4')
llplot(perio.val.enso34$spec, spec.val.enso34$spec)
title('Nino 3.4')
dev.off()











