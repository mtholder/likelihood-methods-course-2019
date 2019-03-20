y1 = c(19, 17)
y2 = c(25, 28)
y3 = c(39, 35)
plot(c(1,1,2,2,3,3), c(y1, y2, y3), xlim=c(0,3), ylim=c(0,40))
abline(a=0, b=MLE.null)
abline(a=muhat-bhat, b=bhat, col="red")
y = sum(y1 + y2 + y3)
MLE.null = y/12

m1 = MLE.null
m2 = 2*MLE.null
m3 = 3*MLE.null

lnl = sum(y1)*log(m1) + sum(y2)*log(m2) + sum(y3)*log(m3) - 2*(m1 + m2 + m3)

lnL.null = lnl
lnL.alt1 = lnl
lnL.alt2 = lnl

s1 = sum(y1)
s2 = sum(y2)
s3 = sum(y3)
sd = s3 -s1
f1 = (sd)/6
f2 = (.5 + s2/(2*s1+sd) + s3/(2*s1+2*sd))
bhat = f1*f2
muhat = 2*s1*bhat/(s3-s1)
m1 = muhat
m2 = m1 + bhat
m3 = m1 +2*bhat

m1 = mean(y1)
m2 = mean(y2)
m3 = mean(y3)