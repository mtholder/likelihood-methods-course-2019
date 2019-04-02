ns = 95
ni = 4
nv = 1
kap = 1
n = ni + nv + ns
pdist = (ni + nv)/n
x = seq(0, 10*pdist, pdist/500)

lnlf = function(nu, kappa) {
  a = exp(-4*nu/(2+kappa))
  b = 2*exp(-(2 + 2*kappa)*nu/(2+kappa))
  return(ns*log(1 + a + b) + ni*log(1 + a - b) + nv*log(1 - a) - (ns+ni+nv)*log(4))
}
kap1lnlf = function(nu) {
  return(lnlf(nu, 1.0))
}
plot(x, kap1lnlf(x), type="l")



num.opt.answer = optimize(kap1lnlf, interval = c(0, 10), maximum = TRUE)

phat = pdist
stderr.phat = sqrt(phat*(1-phat)/n)
moe.p = 1.96*stderr.phat
p.ci = c(phat - moe.p, phat + moe.p)
jc.p.to.nu = function(p) {
  return(-0.75*log(1 - 4*p/3))
}
nu.ci = jc.p.to.nu(p.ci)
nu.hat = jc.p.to.nu(phat)
nu.hat - nu.ci

crit.lnL = num.opt.answer$objective - 1.92

abs.diff.from.crit.lnL = function(nu) {
  llrd = lnlf(nu, 1.0) - crit.lnL;
  return(abs(llrd));
}

mle.nu = num.opt.answer$maximum
ci.width = nu.ci[2] - nu.ci[1]
x = seq(nu.ci[1]/2, 2*nu.ci[2], ci.width/500)
plot(x, abs.diff.from.crit.lnL(x), type="l")
lower.answer = optimize(abs.diff.from.crit.lnL, interval = c(0, mle.nu), maximum = FALSE)
upper.answer = optimize(abs.diff.from.crit.lnL, interval = c(mle.nu, 10*nu.ci[2]), maximum = FALSE)

lrt.nu.ci = c(lower.answer$minimum, upper.answer$minimum)

squared.diff.from.crit.lnL = function(nu) {
  llrd = lnlf(nu, 1.0) - crit.lnL;
  return(llrd^2);
}

plot(x, squared.diff.from.crit.lnL(x), type="l")
lower.answer.from.sqdiff = optimize(squared.diff.from.crit.lnL, interval = c(0, mle.nu), maximum = FALSE)
upper.answer.from.sqdiff = optimize(squared.diff.from.crit.lnL, interval = c(mle.nu, 10*nu.ci[2]), maximum = FALSE)
lrt.nu.ci.from.sqdiff = c(lower.answer$minimum, upper.answer$minimum)


lrt.nu.ci
lrt.nu.ci.from.sqdiff
nu.ci


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

