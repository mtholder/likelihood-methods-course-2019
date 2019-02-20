qgamma(0.025, shape=6, rate=2.2)
qgamma(0.975, shape=6, rate=2.2)

s=0.01
x=seq(s,10,by=s)
b = 0.2
plot(x, b*exp(-b*x), type="line", col="blue", ylim=c(0,1), xlab="nu", ylab="")
lines(x, x^5*exp(-2*x), col="red")
lines(x, dgamma(x, shape=6, rate=2.2), col="black")
legend(6, .8, legend=c("prior", "likelihood", "posterior"),
       col=c("blue", "red", "black"), lty=1:2, cex=0.8)

# p = plot_ly(type='contour',  x=~volcano)
nu = seq(0,12,by=.01)
plot(nu, 5*log(nu), type="l", ylim=c(-5,15), ylab="5 ln(nu)", xlab="nu")
plot(nu, -2*nu, type="l", ylim=c(-25,0), ylab="-2nu", xlab="nu")
plot(nu, (5/nu)-2, type="l", ylim=c(-2,10), ylab="dlnL/dnu", xlab="nu")
abline(h=0, lty=3)
lines(nu, 5/nu, lty=2)
abline(h=2, lty=2)
plot(nu, -2*nu, type="l", ylim=c(-25,0), ylab="-2nu", xlab="nu")
plot(nu, -2*nu, type="l", ylim=c(-25,0), ylab="-2nu", xlab="nu")
plot(nu,  5*log(nu) -2*nu, type="l", ylim=c(-25,0), ylab="log(Pr(x=5|nu))", xlab="nu")

#require(sfsmisc)
s = 0.02
nu = seq(s,7,by=s)
x = nu
y = nu
# xyg = xy.grid(x, y)
nx = length(x)
ny = length(y)
df  = 
z = matrix(nrow=nx, ncol=ny)
zf = matrix(nrow=nx, ncol=ny)
for (i in 1:nx) {
  for (j in 1:ny) {
    z[i,j] = 6*log(x[i] + 2*y[j]) - x[i] - 2*y[j] ;
    zf[i,j] = 5*log(x[i]) -2*x[i] + 6*log(x[i] + 2*y[j]) - x[i] - 2*y[j] ;
  }
}

mz = max(z) -.00001
contour(x, y, z, levels=seq(mz, mz-10, by=-.5) )
mzf = max(zf)-.00001
contour(x, y, zf, levels=seq(mzf, mzf-10, by=-.5) )
