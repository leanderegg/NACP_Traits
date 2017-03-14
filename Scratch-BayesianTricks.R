## some interesting tools from:http://www.unc.edu/courses/2010fall/ecol/563/001/docs/lectures/lecture28.htm



#lme centering function
center.data2 <- function(k) {
  inverts$k <- k
  lme.models <- lme(log(PLD)~I(log(temp)-log(k)), random=~I(log(temp)-log(k))|species, data=inverts)
  c(k,as.numeric(VarCorr(lme.models)[2,3]))
}


# calculate the correlations for a bunch of centering constants
my.results2 <- matrix(NA, ncol=2, nrow=40)
for(k in 1:40) {
  my.results2[k,] <- center.data2(k)
}


#display results graphically
plot(my.results2[,2]~ my.results2[,1], type='l', xlab='centering constant, k', ylab='correlation')
abline(h=0, lty=2, col=4)
#locate smallest correlation
my.results2[abs(my.results2[,2])== min(abs(my.results2[,2])),]