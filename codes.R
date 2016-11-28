#####################################

library(MLSeq)
data(cervical)

set.seed(12349)
ratio=0.7
conditions = factor(rep(c("N","T"), c(29,29)))
ind = sample(58, ceiling(58*ratio), FALSE)

train = cervical[,ind]
test = cervical[,-ind]

tr.cond = conditions[ind]
ts.cond = conditions[-ind]

tmmT = voomDDA.train(counts = train, conditions = tr.cond, normalization = "TMM", TRUE)
tmmF = voomDDA.train(counts = train, conditions = tr.cond, normalization = "TMM", FALSE)

quanT = voomDDA.train(counts = train, conditions = tr.cond, normalization = "quan", TRUE)
quanF = voomDDA.train(counts = train, conditions = tr.cond, normalization = "quan", FALSE)

noneT = voomDDA.train(counts = train, conditions = tr.cond, normalization = "none", TRUE)
noneF = voomDDA.train(counts = train, conditions = tr.cond, normalization = "none", FALSE)

tmmNSC = voomNSC.train(counts = train, conditions = tr.cond, normalization = "TMM")
quanNSC = voomNSC.train(counts = train, conditions = tr.cond, normalization = "quan")
  
table(ts.cond, predict.voomDDA(tmmT, test))
table(ts.cond, predict.voomDDA(tmmF, test))
table(ts.cond, predict.voomDDA(quanT, test))
table(ts.cond, predict.voomDDA(quanF, test))
table(ts.cond, predict.voomDDA(noneT, test))
table(ts.cond, predict.voomDDA(noneF, test))
table(ts.cond, predict.voomNSC(tmmNSC, test))
table(ts.cond, predict.voomNSC(quanNSC, test))
