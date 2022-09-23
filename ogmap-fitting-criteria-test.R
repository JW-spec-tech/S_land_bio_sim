set.seed(1)
dat= tibble(x = seq(0,1,length=100)) %>%
mutate(y1= rnorm(100, x*(1-x), 0.01),
       y2 = rgamma(100,shape = 2,scale = 1))

plot(y1~x, data= dat)

y1_ecdf = ecdf(dat$y1)

y2_ecdf = ecdf(dat$y2)

dat$y1_p = y1_ecdf(dat$y1)
dat$y2_p = y2_ecdf(dat$y2)

#hist(dat$y1_p)
#hist(dat$y2_p)

y1_unif_test = ks.test(x = dat$y1_p,y = punif)$p.value
y2_unif_test = ks.test(x = dat$y2_p,y = punif)$p.value


y1_unif_test 

p1_diff  = outer(dat$y1_p, dat$y1_p,FUN = "-")
p1_diff = abs(p1_diff)
p1_diff = diag(p1_diff[-1,])



p2_diff  = outer(dat$y2_p, dat$y2_p,FUN = "-")
p2_diff = abs(p2_diff)
p2_diff = diag(p2_diff[-1,])

ptri = function(x) x*(2-x)

y1_indep_test = ks.test(x = p1_diff, y = ptri)$p.value
y2_indep_test = ks.test(x = p2_diff, y = ptri)$p.value


print(y1_indep_test) 
print(y2_indep_test) 

#hist(p1_diff)
#hist(p2_diff)

