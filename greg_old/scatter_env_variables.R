library(solaR)
library(hexbin)
 
 
splom(X,
 panel=panel.hexbinplot, diag.panel = function(x, ...){
	 yrng <- current.panel.limits()$ylim
	 d <- density(x, na.rm=TRUE)
	 d$y <- with(d, yrng[1] + 0.95 * diff(yrng) * y / max(y) )
	 panel.lines(d)
	 diag.panel.splom(x, ...)
	},
	#lower.panel = NULL,
	lower.panel = function(x, y, ...){
		panel.hexbinplot(x, y, ...)
		panel.loess(x, y, ..., col = 'red')
	},
	 varname.cex=0.7,type='p',cex=0.5)	
 

 
xy_dens(d$sst,d$t_avg)

xy_dens(d$sst,log(d$npp))

par(mar=c(2,2,2,2))
hexbinplot(d$npp ~ d$sst)

plot(d$sst,log(d$npp))
plot(d$sst,d$t_avg)
plot(d$o2_avg,d$t_avg)


 