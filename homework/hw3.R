p = 0.00125
p = 0.001
k = 5
M = 2000
n = 2
lnLwobc = k*log(p) + (n*M - k)*log(1-p) 
bc = 2000*1999*2000*1999*1998/(12)
lnL = lnLwobc + log(bc)
log(bc)
plnl = lnL

2*(plnl - lnL)
dbinom(2, 2000, 0.00125, log=T)+dbinom(3, 2000, 0.00125, log=T)
