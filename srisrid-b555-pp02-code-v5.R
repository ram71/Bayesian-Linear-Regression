# install.packages("ggplot2")
library(ggplot2)

# set working directory
setwd("E:\\iu\\fa18\\b555\\pp02\\pp2data")


######################################### TASK 1 #########################################

conv.matr = function(x)
{
  rownames(x) = NULL
  colnames(x) = NULL
  return(as.matrix(x))
}

lm.regul = function(lambda, y, X)
{
  return(solve(lambda * diag(1, dim(X)[2]) + t(X) %*% X, tol = 1e-80) %*% t(X) %*% y)
}

lm.ml = function(lambda, y, X)
{
  return(solve(t(X) %*% X, tol = 1e-80) %*% t(X) %*% y)
}

lm.pred = function(what, y, X)
{
  return((t(y - (X %*% what)) %*% (y - (X %*% what)))/dim(y)[1])
}

train.files.y = c("trainR-100-10.csv", "trainR-100-100.csv", "trainR-1000-100.csv", "trainR-crime.csv", "trainR-wine.csv")
train.files.x = c("train-100-10.csv", "train-100-100.csv", "train-1000-100.csv", "train-crime.csv", "train-wine.csv")
test.files.x = c("test-100-10.csv", "test-100-100.csv", "test-1000-100.csv", "test-crime.csv", "test-wine.csv")
test.files.y = c("testR-100-10.csv", "testR-100-100.csv", "testR-1000-100.csv", "testR-crime.csv", "testR-wine.csv")


lambda = 0:150
label.train = c()
mse.train = c()
label.test = c()
mse.test = c()

for(i in 1:length(train.files.x))
{
  train.x = conv.matr(read.csv(train.files.x[i], header = FALSE))
  train.y = conv.matr(read.csv(train.files.y[i], header = FALSE))
  what = sapply(0:150, train.y, train.x, FUN = lm.regul)
  mse.train = c(mse.train, apply(what, train.y, train.x, MARGIN = 2, FUN = lm.pred))
  label.train = c(label.train, rep(train.files.x[i], length(lambda)))
  
  test.x = conv.matr(read.csv(test.files.x[i], header = FALSE))
  test.y = conv.matr(read.csv(test.files.y[i], header = FALSE))
  mse.test = c(mse.test, apply(what, test.y, test.x, MARGIN = 2, FUN = lm.pred))
  label.test = c(label.test, rep(test.files.x[i], length(lambda)))
}

(df.train = data.frame(lambda = rep(lambda , length(train.files.x)), label.train = as.factor(label.train), mse.train = mse.train))

# MSE plots for train datasets

ggplot(df.train, aes(x = lambda, y = mse.train, color = label.train)) + geom_point() + geom_line() + ggtitle("plot of mse vs lambda for train") + ylim(0, 10)

(df.test = data.frame(lambda = rep(lambda , length(train.files.x)), label.test = as.factor(label.test), mse.test = mse.test))

# MSE plots for test datasets
# NOTE: I am not displaying the MSE value for 100-100 for the MLE case as it is very extreme and ruins the plot

ggplot(df.test, aes(x = lambda, y = mse.test, color = label.test)) + geom_point() + geom_line() + ggtitle("plot of mse vs lambda for test") +  ylim(0, 10)

# pick out least MSE values

merge(df.test, aggregate(mse.test ~ label.test, data = df.test, FUN = min))

######################################### END OF TASK 1 #########################################



######################################### TASK 2 #########################################

train.x = conv.matr(read.csv(train.files.x[3], header = FALSE))
train.y = conv.matr(read.csv(train.files.y[3], header = FALSE))
test.x = conv.matr(read.csv(test.files.x[3], header = FALSE))
test.y = conv.matr(read.csv(test.files.y[3], header = FALSE))

rng = seq(10, 800, 3)

smp = function(i, vec)
{
  replicate(20, sample(x = vec, size = i))
}

s = sapply(rng, 1:dim(train.x)[1], FUN = smp)

calc.mse = function(rowind, train.y, train.x, lambda, test.y, test.x)
{
  what = lm.regul(lambda, train.y[rowind, ], train.x[rowind, ])
  return(lm.pred(what, test.y, test.x))
}

calc.mean.mse = function(a, train.y, train.x, lambda, test.y, test.x)
{
  return(mean(apply(a, train.y, train.x, lambda, test.y, test.x, MARGIN = 2, FUN = calc.mse)))
}

dat1 = unlist(lapply(s, train.y, train.x, lambda = 6, test.y, test.x, FUN = calc.mean.mse))
dat2 = unlist(lapply(s, train.y, train.x, lambda = 27, test.y, test.x, FUN = calc.mean.mse))
dat3 = unlist(lapply(s, train.y, train.x, lambda = 80, test.y, test.x, FUN = calc.mean.mse))

df.learn = data.frame(mean.mse = c(dat1, dat2, dat3), train.size = rep(rng, 3), lambda = as.factor(rep(c(6, 27, 80), each = length(rng))))

# make plot to compare mean MSE against sample size for different lambda

ggplot(df.learn, aes(x = train.size, y = mean.mse, color = lambda)) + geom_point() + geom_line() + ggtitle("MSE vs train size plot by sampling")

# zoomed in version of above plot

ggplot(df.learn, aes(x = train.size, y = mean.mse, color = lambda)) + geom_point() + geom_line() + ggtitle("MSE vs train size plot by sampling - zoomed") + ylim(4, 8)

######################################### END OF TASK 2 #########################################



######################################### TASK 3 #########################################

evidence.method = function(X, y)
{
  alpha0 = 5
  beta0 = 5
  p = dim(X)[2]
  n = dim(X)[1]
    
  sn = solve(alpha0 *  diag(rep(1, p)) + beta0 * t(X) %*% X, tol = 1e-80)
  mn = beta0 * sn %*% t(X) %*% y
  e.values = eigen(beta0 * t(X) %*% X)$values
  gamm = sum(sapply(e.values, FUN = function(x){(x/(alpha0 + x))}))
  alpha = as.numeric(gamm/(t(mn) %*% mn))
  beta = as.numeric(1/((1/(n - gamm)) * t(y - X %*% mn) %*% (y - X %*% mn)))
  
  while((abs(alpha0 - alpha) >= 1e-4) || (abs(beta0 - beta) >= 1e-4))
  {
    alpha0 = alpha
    beta0 = beta
    sn = solve(alpha0 *  diag(rep(1, p)) + beta0 * t(X) %*% X, tol = 1e-80)
    mn = beta0 * sn %*% t(X) %*% y
    e.values = eigen(beta0 * t(X) %*% X)$values
    gamm = sum(sapply(e.values, FUN = function(x){(x/(alpha0 + x))}))
    alpha = as.numeric(gamm/(t(mn) %*% mn))
    beta = as.numeric(1/((1/(n - gamm)) * t(y - X %*% mn) %*% (y - X %*% mn)))
  }
  
  A = alpha *  diag(rep(1, p)) + beta * t(X) %*% X  
  sn = solve(A, tol = 1e-80)
  mn = beta * sn %*% t(X) %*% y
  emn = (beta/2) * t(y - X %*% mn) %*% (y - X %*% mn) + (alpha/2) * t(mn) %*% mn
  L = (p/2) * log(alpha) + (n/2) * log(beta) - emn - (1/2) * log(det(A)) - (n/2) * log(2 * pi)
  w.ml = solve(t(X) %*% X, tol = 1e-80) %*% t(X) %*% y
  
  return(list("mn" = mn, "sn" = sn, "alpha" = alpha, "beta" = beta, "log.evidence" = L, "w.ml" = w.ml))
}

mse.test = rep(NA, length(test.files.x))
label.test = c()

train.files.y = c("trainR-100-10.csv", "trainR-100-100.csv", "trainR-1000-100.csv", "trainR-crime.csv", "trainR-wine.csv")
train.files.x = c("train-100-10.csv", "train-100-100.csv", "train-1000-100.csv", "train-crime.csv", "train-wine.csv")
test.files.x = c("test-100-10.csv", "test-100-100.csv", "test-1000-100.csv", "test-crime.csv", "test-wine.csv")
test.files.y = c("testR-100-10.csv", "testR-100-100.csv", "testR-1000-100.csv", "testR-crime.csv", "testR-wine.csv")

for(i in 1:length(test.files.x))
{
  train.x = conv.matr(read.csv(train.files.x[i], header = FALSE))
  train.y = conv.matr(read.csv(train.files.y[i], header = FALSE))
  test.x = conv.matr(read.csv(test.files.x[i], header = FALSE))
  test.y = conv.matr(read.csv(test.files.y[i], header = FALSE))
  mse.test[i] = lm.pred(evidence.method(train.x, train.y)$mn, test.y, test.x)
  label.test = c(label.test, test.files.x[i])
}

df.evid = data.frame(label.test, mse.test, method = "evidence")
df.temp = merge(df.test, aggregate(mse.test ~ label.test, FUN = min, data = df.test))
df.temp$method = "minimum-MSE" 

# print table to compare

(df.compare = rbind(df.evid, df.temp[, c("mse.test", "label.test", "method")]))

# % change with respect to best result

data.frame(label.test, "abs percent diff" = abs(((df.temp$mse.test - df.evid$mse.test)/df.temp$mse.test) * 100))

# plot bar graph to compare performance

ggplot(df.compare, aes(x = label.test, y = mse.test, fill = method)) + geom_bar(stat = "identity", position = "dodge") + ggtitle("compare best test result with evidence method")

######################################### END OF TASK 3 #########################################



######################################### TASK 4 #########################################

gen.poly = function(d, X)
{
  f.pow = function(pow, x)
  {
    return(x ^ pow)
  }
  
  return(cbind(1, sapply(1:d, X, FUN = f.pow)))
}

train.files.y = c("trainR-f3.csv", "trainR-f5.csv")
train.files.x = c("train-f3.csv", "train-f5.csv")
test.files.x = c("test-f3.csv", "test-f5.csv")
test.files.y = c("testR-f3.csv", "testR-f5.csv")

log.evid = c()
label.test = c()
mse.test.reg = c()
mse.test.ml = c()
d = c()

for(i in 1:length(test.files.x))
{
  train.x = conv.matr(read.csv(train.files.x[i], header = FALSE))
  train.y = conv.matr(read.csv(train.files.y[i], header = FALSE))
  
  train.list.x = sapply(1:10, train.x, FUN = gen.poly)
  ev.list = lapply(train.list.x, train.y, FUN = evidence.method)
  
  test.x = conv.matr(read.csv(test.files.x[i], header = FALSE))
  test.y = conv.matr(read.csv(test.files.y[i], header = FALSE))
  test.list.x = sapply(1:10, test.x, FUN = gen.poly)
  
  for(j in 1:10)
  {
    log.evid = c(log.evid, ev.list[[j]]$log.evidence)
    label.test = c(label.test, test.files.x[i])
    d = c(d, j)
    mse.test.reg = c(mse.test.reg, lm.pred(ev.list[[j]]$mn, test.y, test.list.x[[j]]))
    mse.test.ml = c(mse.test.ml, lm.pred(ev.list[[j]]$w.ml, test.y, test.list.x[[j]]))
  }
}

df.plot.regul = data.frame("dataset" = label.test, "mse" = mse.test.reg, "d" = d, "type" = "regularized")
df.plot.ml = data.frame("dataset" = label.test, "mse" = mse.test.ml, "d" = d, "type" = "ml")
df.plot.1 = rbind(df.plot.regul, df.plot.ml)

# plot MSE for MLE and Bayesian

ggplot(df.plot.1[df.plot.1$dataset == "test-f3.csv",], aes(x = d, y = mse, color = type)) + geom_point() + geom_line() + facet_wrap(~ type) + ggtitle("plot of mse vs d for f3")
ggplot(df.plot.1[df.plot.1$dataset == "test-f5.csv",], aes(x = d, y = mse, color = type)) + geom_point() + geom_line() + facet_wrap(~ type) + ggtitle("plot of mse vs d for f5")

# plot log evidence

df.plot.2 = data.frame("dataset" = label.test, "log.evidence" = log.evid, "d" = d)
ggplot(df.plot.2, aes(x = d, y = log.evidence)) + geom_point() + geom_line() + facet_wrap(~ dataset) + ggtitle("plot of log evidence vs d")


# look at differences b/w mle and bayesian methods

df.comp = merge(df.plot.regul, df.plot.ml, by = c("dataset", "d"))
df.comp$delta = df.comp$mse.x - df.comp$mse.y
df.comp$percent.diff = (df.comp$delta/df.comp$mse.y) * 100
df.comp = df.comp[order(df.comp$dataset, df.comp$d), ]
df.comp[, c("dataset", "d", "percent.diff")]

# see if evidence method works

merge(df.plot.1[df.plot.1$type == "regularized",], aggregate(mse ~ dataset + type, data = df.plot.1[df.plot.1$type == "regularized",], FUN = min))
merge(df.plot.2, aggregate(log.evidence ~ dataset, data = df.plot.2, FUN = max))

######################################### END OF TASK 4 #########################################