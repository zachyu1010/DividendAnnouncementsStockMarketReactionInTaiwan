# commands
# .libPaths() : Get library locations containing R packages
# library() : Get the list of all the packages installed
# search() : Get all packages currently loaded in the R environment
# install.packages("Package Name") : Install directly from CRAN

rm(list = ls())
library(zoo)
library(eventstudies)

# parameters, number of the firms
nfirms = 50
date.shift = 1

# parmeters, time frames
time.t1 = -5
time.t2 = 5
time.t0 = -125
time.t3 = 125

# Read the prices of stocks from CSV file
data.dir <- "/Users/McInnis/MEGA/Exeter/Courses/Term2/Financial Econometrics/Final Project/R study event/"
data.filename <- "stock prices.csv"
data.filepath <- paste(data.dir, data.filename, sep="")
obser.s.prices <- read.table(data.filepath, header=TRUE, sep=",")

# Caculate approximate returns with LN(P1/P2) -> returns
obser.s.returns <- log(obser.s.prices[-nrow(obser.s.prices),-1]/obser.s.prices[-1,-1])
obser.s.returns <- data.frame(obser.s.returns)
obser.s.returns <- cbind(obser.s.prices[,1][-nrow(obser.s.prices)], obser.s.returns)
names(obser.s.returns) <- c("Date", names(obser.s.returns)[-1])

# Time series of stock's returns -> obser.s.return
obser.s.returns <- zoo(obser.s.returns[-1], 
                  as.Date(sapply(obser.s.returns[1], as.character)[,1], "%b-%d-%Y"))

# Read the event dates from CSV file
data.dir <- "/Users/McInnis/MEGA/Exeter/Courses/Term2/Financial Econometrics/Final Project/R study event/"
data.filename <- "announcement dates.csv"
data.filepath <- paste(data.dir, data.filename, sep="")
obser.s.dates <- read.table(data.filepath, header=TRUE, sep=",")

# Prepare announcement data and categorize the announcement
t_sep = rep(c(TRUE,FALSE), ncol(obser.s.dates)/2)
t_div = obser.s.dates[,!t_sep]
t_news = rbind(t_div[-nrow(t_div),] - t_div[-1,], NA)
t_when = obser.s.dates[,t_sep][!is.na(t_news)]
t_unit = t(matrix(rep(names(t_news), nrow(t_news)), ncol(t_news), nrow(t_news)))[!is.na(t_news)]

obser.s.dates = data.frame(unit = t_unit, when=t_when)
obser.s.dates$unit <- sapply(obser.s.dates$unit, as.character)
# add 1 for delay repsonse to the announcement because of the publishing time after market close.
obser.s.dates$when <- as.Date(sapply(obser.s.dates$when, as.character), "%d/%m/%Y") + date.shift

goodBadNews <- function(x) {
  if(x == 0) {
    ret = "no"
  } else {
    if(x > 0){
      ret = "good"
    } else{
      ret = "bad"
    }
  }
  
  ret
}

obser.ret.category = sapply(t_news[!is.na(t_news)], goodBadNews)

# Read the prices of market from CSV file
data.dir <- "/Users/McInnis/MEGA/Exeter/Courses/Term2/Financial Econometrics/Final Project/R study event/verification_data/"
data.filename <- "TAIEX.csv"
data.filepath <- paste(data.dir, data.filename, sep="")
obser.m.prices <- read.table(data.filepath, header=TRUE, sep=",")

obser.m.returns <- log(obser.m.prices[-nrow(obser.m.prices),-1]/obser.m.prices[-1,-1])
obser.m.returns <- data.frame(obser.m.returns)
obser.m.returns <- cbind(sapply(obser.m.prices[,1], as.character)[-nrow(obser.m.prices)], obser.m.returns)
names(obser.m.returns) <- c("Date", names(obser.m.prices)[-1])

# Time series of market's returns -> obser.m.return
obser.m.returns <- zoo(obser.m.returns[-1], 
                      as.Date(sapply(obser.m.returns[1], as.character)[,1], "%b-%d-%Y"))

obser.m.dates <- obser.s.dates
obser.m.dates$unit <- rep(names(obser.m.prices)[2], nrow(obser.m.dates))


# estimation window of stocks
obser.s.es <-  phys2eventtime(z=obser.s.returns, events=obser.s.dates, width=time.t3)
obser.s.es.w.est <- window(obser.s.es$z.e, start=time.t0, end=time.t1-1)

# estimation window of market
obser.m.es <-  phys2eventtime(z=obser.m.returns, events=obser.m.dates, width=time.t3)
obser.m.es.w.est <- window(obser.m.es$z.e, start=time.t0, end=time.t1-1)

# regress stocks's return on market return
matreg <- function(x){
  z = lm(obser.s.es.w.est[,x] ~ obser.m.es.w.est[,x])
  names(z["coefficients"][[1]]) <- c("alpha", "beta")
  z
}

obser.estimates = sapply(1:ncol(obser.s.es.w.est), matreg)

# caculate abnormal returns
obser.s.es.w.tline <- window(obser.s.es$z.e, start=time.t0, end=time.t3)
obser.m.es.w.tline <- window(obser.m.es$z.e, start=time.t0, end=time.t3)

abret <- function(x){
  beta = obser.estimates["coefficients", x][[1]]["beta"]
  alpha = obser.estimates["coefficients", x][[1]]["alpha"]
  
  obser.s.es.w.tline[,x] - alpha - beta*obser.m.es.w.tline[,x]
}

abret = sapply(1:ncol(obser.estimates), abret)
abret <- zoo(abret, c(time.t0:time.t3))

isolate_news_and_average <- function(s, x, x.all, type, idx){
  t <- s[, obser.ret.category == type]
  t.all <- if(is.null(dim(t))){
                                dim(t) <- c(length(idx), 1)
                                t
                              } else {
                                w = zoo(apply(t, 1, mean), idx)
                                dim(w) <- c(length(idx), 1)
                                w
                              }
  
  assign(x, t, envir=baseenv())
  assign(x.all, t.all, envir=baseenv())
}

abret.w.estimate <- window(abret, start=time.t0, end=time.t1-1)
isolate_news_and_average(abret.w.estimate, "abret.w.estimate.new.good", "abret.w.estimate.new.good.all", "good", c(time.t0:(time.t1-1)))
isolate_news_and_average(abret.w.estimate, "abret.w.estimate.new.bad", "abret.w.estimate.new.bad.all", "bad", c(time.t0:(time.t1-1)))
isolate_news_and_average(abret.w.estimate, "abret.w.estimate.new.no", "abret.w.estimate.new.no.all", "no", c(time.t0:(time.t1-1)))

abret.w.event = window(abret, start=time.t1, end=time.t2)
isolate_news_and_average(abret.w.event, "abret.w.event.new.good", "abret.w.event.new.good.all", "good", c(time.t1:time.t2))
isolate_news_and_average(abret.w.event, "abret.w.event.new.bad", "abret.w.event.new.bad.all", "bad", c(time.t1:time.t2))
isolate_news_and_average(abret.w.event, "abret.w.event.new.no", "abret.w.event.new.no.all", "no", c(time.t1:time.t2))

abret.w.event.news <- merge(abret.w.event.new.good.all, 
                            abret.w.event.new.bad.all,
                            abret.w.event.new.no.all,
                            all = TRUE,
                            suffixes = c("Good News", "Bad News", "No News"))

abret.w.event.all <- zoo(apply(abret.w.event, 1, mean), c(time.t1:time.t2))
# abret.w.event.cum <- remap.cumsum(abret.w.event, is.pc=FALSE, base = 0)
# abret.w.event.cum.all <- zoo(apply(abret.w.event.cum, 1, mean), c(time.t1:time.t2))
abret.w.event.new.good.all.cum <- remap.cumsum(abret.w.event.new.good.all, is.pc=FALSE, base = 0)
abret.w.event.new.bad.all.cum <- remap.cumsum(abret.w.event.new.bad.all, is.pc=FALSE, base = 0)
abret.w.event.new.no.all.cum <- remap.cumsum(abret.w.event.new.no.all, is.pc=FALSE, base = 0)

abret.w.event.cum.news <- merge(abret.w.event.new.good.all.cum, 
                               abret.w.event.new.bad.all.cum,
                               abret.w.event.new.no.all.cum,
                               all = TRUE,
                               suffixes = c("Good News", "Bad News", "No News"))


# Plot the event windows associated good, bad and no news
plot(abret.w.event.cum.news, type="l", plot.type = "single", lty = 1:3, col = 4:2)
legend("topleft", legend = c("Good News", "Bad News", "No News"), col=4:2, pch=1)
#plot(abret.w.event.cum.news, type="l", plot.type = "multiple")

# test statistic cacluation

# CAAR
test.car.good = tail(abret.w.event.new.good.all.cum, 1)
test.car.bad = tail(abret.w.event.new.bad.all.cum, 1)
test.car.no = tail(abret.w.event.new.no.all.cum, 1)

# Variance of average abnormal return
test.var.good = var(abret.w.estimate.new.good.all)
test.var.bad = var(abret.w.estimate.new.bad.all)
test.var.no = var(abret.w.estimate.new.no.all)

# var(CAAR)
test.var.good.car = test.var.good * (time.t2-time.t1+1)
test.var.bad.car = test.var.bad * (time.t2-time.t1+1)
test.var.no.car = test.var.no * (time.t2-time.t1+1)

# statistic 1, CAAR over (var(CAAR))^(1/2) 
test.stat1.good = test.car.good / sqrt(test.var.good.car)
test.stat1.bad = test.car.bad / sqrt(test.var.bad.car)
test.stat1.no = test.car.no / sqrt(test.var.no.car)

# pvalue for statistic 1

pvalue <- function(stat){
  pt(-abs(stat), 10000)*2
}

test.pvalue1.good = pvalue(test.stat1.good)
test.pvalue1.bad = pvalue(test.stat1.bad)
test.pvalue1.no = pvalue(test.stat1.no)

# statistc 2, sign test
signtest <- function(e.all, x, x.all){
  e.mat.all = as.matrix(e.all)
  x.mat.all = as.matrix(x.all)
  ps = length(e.mat.all[e.all < 0])/(time.t1-time.t0)
  p = length(x.mat.all[x.all < 0])/(time.t2-time.t1+1)
  N = ifelse(is.null(dim(x)),1 , ncol(x))
  
  (p-ps)/sqrt(ps*(1-ps)/N)  
}

test.stat2.good = signtest(abret.w.estimate.new.good.all, abret.w.event.new.good, abret.w.event.new.good.all)
test.stat2.bad = signtest(abret.w.estimate.new.bad.all, abret.w.event.new.bad, abret.w.event.new.bad.all)
test.stat2.no = signtest(abret.w.estimate.new.no.all, abret.w.event.new.no, abret.w.event.new.no.all)

# pvalue for statistic 2
test.pvalue2.good = pvalue(test.stat2.good)
test.pvalue2.bad = pvalue(test.stat2.bad)
test.pvalue2.no = pvalue(test.stat2.no)

# statistc 3, rank test
rankevent <- function(x){
  ceiling((x-min(x))/((max(x)-min(x))/(time.t2-time.t1+1)))
}

# skfunc <- function(x){
#   N = ifelse(is.null(dim(x)),1 , ncol(x))
#   L = time.t2 - time.t1 + 1
#   s1 = 0
#   
#   for(i in 1:L){
#     s2 = 0
#     for(j in 1:N){
#       s2 = s2 + x[i, j][[1]]-(L+1)/2
#     }
#     
#     s1 = s1 + (s2/N)^2
#   }
#   
#   sqrt(s1/L)
# }
# 
# ranktest <- function(x){
#   N = ifelse(is.null(dim(x)),1 , ncol(x))
#   L = time.t2 - time.t1 + 1
#   x.rank = rankevent(x)
#   sk = skfunc(x.rank)
#   i0 = (L+1)/2
#   
#   stat = 0
#   for(j in 1:N){
#     stat = stat + x.rank[i0, j][[1]] - (L+1)/2
#   }
#   
#   stat = (stat/N)/sk
# }

skfunc2 <- function(x){
  N = ifelse(is.null(dim(x)),1 , ncol(x))
  L = time.t2 - time.t1 + 1
  s1 = 0
  
  sqrt(sum((rowSums((x-(L+1)/2)/N))^2)/L)
}

ranktest2 <- function(x){
  N = ifelse(is.null(dim(x)),1 , ncol(x))
  L = time.t2 - time.t1 + 1
  x.rank = rankevent(x)
  sk = skfunc2(x.rank)
  i0 = (L+1)/2

  (1/N)*(sum((x.rank[i0, ] - (L+1)/2)))/sk
}


test.stat3.good = ranktest2(abret.w.event.new.good)
test.stat3.bad = ranktest2(abret.w.event.new.bad)
test.stat3.no = ranktest2(abret.w.event.new.no)

test.pvalue3.good = pvalue(test.stat3.good)
test.pvalue3.bad = pvalue(test.stat3.bad)
test.pvalue3.no = pvalue(test.stat3.no)


