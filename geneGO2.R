library(topGO)
library(limma)
library(edgeR)
library(org.Hs.eg.db)

data = read.table("lung.txt", header = T, sep = "\t")
data2 = data[data$X != "SLC35E2",]
rownames(data2) = data2[,1]
data = data2[,-1]
condition = factor(rep(c("T1", "T2"), c(576, 552)))
design = model.matrix(~condition)

dge = DGEList(counts=as.matrix(data), group = condition)
dge = calcNormFactors(dge, method = "TMM") #webtool'da var. RLE dediği deseq, TMM dediği TMM, none dediği none
v = voom(dge,design,plot=F) 
fit = lmFit(v,design)
fit = eBayes(fit)
res = topTable(fit,coef=ncol(design),number = dim(data)[1])
geneList = res$adj.P.Val #voomNSC tarafından seçilen genler
names(geneList) = rownames(res) #voomNSC tarafından seçilen genler
biomarkers = rownames(res[10:50,]) #voomNSC tarafından seçilen genler
#truefalse = is.element(names(geneList),biomarkers)
#selection = function(x) TRUE 

#allGO2genes = annFUN.org(whichOnto='BP', feasibleGenes = NULL,
#                         mapping="org.Hs.eg.db", ID = "genename")

#GOdata = new("topGOdata", description = "Simple session", ontology = "BP",
# allGenes = geneList, geneSel = topDiffGenes, nodeSize = 10,
# annot = annFUN.db, affyLib = affyLib)
truefalse = function(allScore) {  
  truefalse = is.element(names(geneList),biomarkers)
  return(truefalse)
}

allGO2genes = annFUN.org(whichOnto='BP', feasibleGenes = NULL,
                         mapping="org.Hs.eg.db", ID = "symbol")
GOdata =  new('topGOdata', ontology = 'BP', allGenes = geneList, 
              annot = annFUN.GO2genes, GO2genes = allGO2genes, 
              geneSel = truefalse, nodeSize=10)

results.ks = runTest(GOdata, algorithm = "classic", statistic = "ks")
results.ks
allRes = GenTable(GOdata, KS = results.ks, orderBy = "KS", topNodes = 20)
allRes[,c('GO.ID','Term','KS')]
showSigOfNodes(GOdata, score(results.ks), firstSigNodes = 5, useInfo = 'all')

