data(cervical)
data = cervical
conditions = rep(c("N","T"),c(29,29))

n = ncol(data)  # number of samples
p = nrow(data)  # number of features

nTest = ceiling(n*0.3)  # number of samples for test set (20% test, 80% train).
ind = sample(n,nTest,FALSE)

# train set
data.train = data[,-ind]
data.train = as.matrix(data.train + 1)
classtr = conditions[-ind]

data.test = data[,ind]
data.test = as.matrix(data.test + 1)
classts = conditions[ind]
