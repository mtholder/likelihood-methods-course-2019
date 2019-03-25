ns = 95
ni = 4
nv = 1
kap = 1
x = seq(0, .2, 0.001)

lnlf = function(nu, kappa) {
  a = exp(-4*nu/(2+kappa))
  b = 2*exp(-(2 + 2*kappa)*nu/(2+kappa))
  return(ns*log(1 + a + b) + ni*log(1 + a - b) + nv*log(1 - a) - (ns+ni+nv)*log(4))
}
plot(x, lnlf(x, kap), type="l")

kap1lnlf = function(nu) {
    return(lnlf(nu, 1.0))
}

optimize(kap1lnlf, interval = c(0, 10), maximum = TRUE)


s = 0.01
nu = seq(s, 0.25,by=s)
k = seq(10*s, 100 , by=10*s)
x = nu
y = k
nx = length(x)
ny = length(y)
z = matrix(nrow=nx, ncol=ny)
for (i in 1:nx) {
  for (j in 1:ny) {
    z[i,j] = lnlf(x[i], y[j]) ;
  }
}

mz = max(z) -.00001
contour(x, y, z, levels=seq(mz, mz-10, by=-.5) )
neg.lnl = function(theta) {
  return(-1*lnlf(theta[1], theta[2]))
}
INITIAL.PARAMETER.GUESS = c(0.5, 1)
optim(par=INITIAL.PARAMETER.GUESS,
      fn=neg.lnl,
      method="Nelder-Mead");

