shinyServer(function(input, output, session) {
    options(shiny.maxRequestSize=500*1024^2) 
    library(edgeR)
    library(limma)
    library(DESeq2)
    library(caret)
    library(kernlab)
    library(randomForest)
    library(sSeq)
    library(plyr)
    library(pamr)
    library(sfsmisc)
    library(foreach)
    library(digest)
    library(RColorBrewer)
    library(gplots)
    library(d3heatmap)
    library(igraph)
    library(DT)
    library(miRNAtap)
    library(topGO)
    library(org.Hs.eg.db)

    
    source("allFunctions.R")
    
    ## Source codes of PLDA (Witten et. al.)
    source("Classify.cv.R")
    source("Classify.R")
    source("CountDataSet.R")
    source("FindBestTransform.R")
    source("Functions.R")
    source("NullModel.R")
    source("NullModelTest.R")
    source("PoissonDistance.R")
    
    ## Source codes of NBLDA
    source("classnb.R")
    source("CountDataSet1.R")
    source("GetDnb.R")
    
    ## Source codes of voom classifiers
    source("voomGSD.R")
    source("weighted.stats.R")
    source("voomNSC.train.R")
    source("predict.voomNSC.R")
    source("voomDDA.train.R")
    source("predict.voomDDA.R")
    
    observe({
        if (input$varF){
            dat = t(dataM())
            nZ = nearZeroVar(dat,saveMetrics = FALSE)
            updateTextInput(session, inputId = "maxVar", label = "Number of genes with maximum variance", value = as.character(ncol(dat) - length(nZ)))
        } 
    })
    
    
    dataM <- reactive({  ## Data input.
        if(input$dataInput==1){  ## Load example data.
            
            if(input$sampleData==1){
                data <- read.table("cervical_train.txt", header=TRUE)
            }
            
            else if(input$sampleData==2){
                data <- read.table("lung_train.txt", header=TRUE)
            }
        } 
        
        else if(input$dataInput==2){  ## Upload data. Train
            
            inFile <- input$uploadTrain
            mySep <- switch(input$fileSepDF, '1'=",",'2'="\t",'3'=";", '4'="")
            
            if (is.null(input$uploadTrain))  {return(NULL)}
            
            if (file.info(inFile$datapath)$size <= 10485800){
                data <- read.table(inFile$datapath, sep=mySep, header=TRUE, fill=TRUE, na.strings = c("", "NA","."))
            }
            
            else print("File is bigger than 10MB and will not be uploaded.")
        } 
        return(data)
    })
    
    dataC <- reactive({  ## Data input. Conditions
        if(input$dataInput==1){  ## Load example data.
            
            if(input$sampleData==1){
                data <- read.table("cervical_cond.txt", header=FALSE, sep="\t")
            }
            
            else if(input$sampleData==2){
                data <- read.table("lung_cond.txt", header=FALSE, sep="\t")
            }
        }
        
        if(input$dataInput==2){  ## Upload data.
            
            inFile <- input$uploadCond
            mySep <- switch(input$fileSepDF, '1'=",",'2'="\t",'3'=";", '4'="")
            
            if (is.null(input$uploadCond))  {return(NULL)}
            
            if (file.info(inFile$datapath)$size <= 10485800){
                data <- read.table(inFile$datapath, sep=mySep, header=FALSE, fill=TRUE, na.strings = c("", "NA","."), 
                                   stringsAsFactors = TRUE, colClasses = "character")
            }
            
            else print("File is bigger than 10MB and will not be uploaded.")
        } 
        return(data)
    })

    dataT <- reactive({  ## Data input.
        
        if(input$dataInput==1){  ## Load example data.
            
            if(input$sampleData==1){
                data <- read.table("cervical_test.txt", header=TRUE)
            }
            
            else if(input$sampleData==2){
                data <- read.table("lung_test.txt", header=TRUE)
            }
        } 
        
        if(input$dataInput==2){  ## Upload data.
            
            inFile <- input$uploadTest
            mySep <- switch(input$fileSepDF, '1'=",",'2'="\t",'3'=";", '4'="")
            
            if (is.null(input$uploadTest))  {return(NULL)}
            
            if (file.info(inFile$datapath)$size <= 10485800){
                data <- read.table(inFile$datapath, sep=mySep, header=TRUE, fill=TRUE, na.strings = c("", "NA","."))
            }
            
            else print("File is bigger than 10MB and will not be uploaded.")
        } 
        return(data)
    })

    
    output$RawDataTrain <- DT::renderDataTable({
        dataTrain <- dataM()
        conditions <- dataC()[,1]
        dataTrain <- dataTrain[order(rownames(dataTrain)),]
        
        col_dots = rep("...", 10)
        row_dots = rep("...", 15)
        
        ## Buraya data wiev için bir fonsiyon yazılacak.
        rawTrain = rbind(data.frame(dataTrain[1:8, c(1:12)], col_dots = rep("...", 8), 
                         dataTrain[1:8, c(dim(dataTrain)[2]-1, dim(dataTrain)[2])]), 
              row_dots, cbind(dataTrain[c(dim(dataTrain)[1]-1,dim(dataTrain)[1]), 1:12], col_dots = c("...","..."), 
                              dataTrain[c(dim(dataTrain)[1]-1,dim(dataTrain)[1]), c(dim(dataTrain)[2]-1, 
                                                                                    dim(dataTrain)[2])]))
        return(rawTrain)
        
        
        
    })
    
    
    output$RawDataTest <- DT::renderDataTable({
        dataTest <- dataT()
        dataTest <- dataTest[order(rownames(dataTest)),]
        
        return(dataTest)
        
        
        
    })
    
    
    output$trainConsole <- renderPrint({
        
        if(input$runVoomDDA){ 
        dataTrain <- dataM()
        conditions <- as.factor(dataC()[,1])
        dataTest <- dataT()
        dataTrain <- dataTrain[order(rownames(dataTrain)),]
        dataTest <- dataTest[order(rownames(dataTest)),]
        
        ## Control steps
        trGenes <- sort(rownames(dataTrain))
        tsGenes <- sort(rownames(dataTest))
        
        if (!identical(trGenes, tsGenes)) stop(warning("Gene names should be identical for both training and test data sets."))
        if (ncol(dataTrain) != length(conditions)) stop(warning("Number of conditions and sample size do not match."))
        
        ## Filtering
        
        if (input$nearZeroF){
            nZ = nearZeroVar(t(dataTrain),saveMetrics = FALSE)
            if(length(nZ) != 0){
                dataTrain = dataTrain[-nZ,]
                dataTest = dataTest[-nZ,]      
            }
        }
        
        if (input$varF){
            ## Variance Filtering
            vars = apply(log(dataTrain+1),1,var)
            idx = order(vars, decreasing = TRUE)
            dataTrain <- dataTrain[idx[1:input$maxVar],]
            dataTest <- dataTest[rownames(dataTrain),]
        }
        
        ###
        if (input$models == "vDLDA") model = voomDDA.train(counts = dataTrain, conditions = conditions, normalization = input$normMeth, TRUE)
        if (input$models == "vDQDA") model = voomDDA.train(counts = dataTrain, conditions = conditions, normalization = input$normMeth, FALSE)
        if (input$models == "vNSC") model = voomNSC.train(counts = dataTrain, conditions = conditions, normalization = input$normMeth)
        
        cat("Model Summary:", "\n")
        cat("-----------------------------------------","\n")
        cat(paste("Raw Data"),"\n")
        cat(paste("   Data includes the read counts of ", dim(dataM())[1], " genes belong to ", dim(dataM())[2], " observations.", sep=""),"\n\n")
        
        if (input$nearZeroF){
            cat(paste("Near-zero filtering"), "\n")
            cat(paste("   ", length(nZ), " out of ", dim(dataM())[1]," genes are filtered.", sep=""), "\n\n")
        }
        
        if (input$varF){
            cat(paste("Variance filtering"), "\n")
            cat(paste("   ", input$maxVar, " genes are selected based on their maximum variance.",sep=""), "\n")
        }
        cat("-----------------------------------------","\n\n")
        
        
        if ("trSummary" %in% input$advOpts){
            if (input$models != "vNSC") trainSummary <- caret::confusionMatrix(table(Actual = conditions, Predicted = predict.voomDDA(model, dataTrain)))
            if (input$models == "vNSC") trainSummary <- caret::confusionMatrix(table(Actual = conditions, Predicted = predict.voomNSC(model, dataTrain)))
            
            cat("Training Summary:","\n")
            cat("-----------------------------------------","\n")
            caret:::print.confusionMatrix(trainSummary)
            cat("-----------------------------------------","\n\n")
        }

        if (input$models != "vNSC") Predicted = predict.voomDDA(model, dataTest)
        if (input$models == "vNSC") Predicted = predict.voomNSC(model, dataTest)
        
        cat("Predictions:","\n")
        cat("-----------------------------------------","\n")
        cat(paste(c("  ", as.character(Predicted)), sep=""), "\n\n")
        
        if ("selGenes" %in% input$advOpts){
        if (input$models == "vNSC"){
            selectedGenes <- model$SelectedGenes[[which(model$threshold == model$opt.threshold)]]
            cat(paste("Selected Genes (voomNSC): ", length(selectedGenes), " out of ", dim(dataTrain)[1], " genes are selected", sep=""),"\n")
            cat("-----------------------------------------","\n")
            cat(paste(selectedGenes,"\n"))
        }}


        heatMapData <- reactive({


                if (input$models == "vNSC"){
                    nSelFeat <- length(selectedGenes)
                    if (nSelFeat < 2) stop(warning("At least 2 features should be selected in the model. Heatmap can not be drawn with 1 feature."))
                }

                dataTrain_HM <- dataM()
                HM_data <- t(log2(t(dataTrain_HM + 0.5)/(apply(dataTrain_HM, 2, sum) + 1) * 1e+06))
                
                if (input$models == "vNSC") HM_data <- HM_data[selectedGenes,]
                if (input$models != "vNSC") HM_data <- HM_data[rownames(dataTrain),]
                
                if (input$centerHeat) HM_data <- scale(as.matrix(HM_data), scale = FALSE)
                hmcol = colorRampPalette(brewer.pal(8, "GnBu"))(250)

            
                    return(HM_data)
        })
        




        output$heatMap <- renderD3heatmap({
            if ((input$runHeatMap)){
				
                #heatmap.2(HM_data, col = hmcol, trace="none", margin=c(10, 6))
                if(input$darkTheme){
                    d3heatmap(heatMapData(), theme = "dark", color = input$colorsHM, width = "800px", height = "1800px")
                }else{

                    d3heatmap(heatMapData(), theme = "", color = input$colorsHM, width = "800px", height = "1800px")

                }
            }
            
        })#,  height = 800, width = 800)  


        output$newtwork <- renderPlot({

            if ((input$runNetwork)){
                cor_mat = cor = cor(t(heatMapData()))
                diag(cor_mat)<-0
                graph<-graph.adjacency(cor_mat,weighted=TRUE,mode="upper")
                E(graph)[ weight>0.7 ]$color <- "green" 
                E(graph)[ weight < -0.7 ]$color <- "red" 
                E(graph)[ weight>0.6 & weight < 0.7 ]$color <- "black" 
                E(graph)[ weight< -0.6 & weight > -0.7 ]$color <- "black" 
                plot(graph)
            }

        }, height = 650, width = 650)



       if(input$models == "vNSC"){
            observe({
                     updateSelectInput(session, "showResults", choices = selectedGenes[3:length(selectedGenes)], selected = selectedGenes[3:length(selectedGenes)][1])
            })
        }

            observe({
                if(input$goAlgorithm == "classic" || input$goAlgorithm == "elim" || input$goAlgorithm == "weight01" || input$goAlgorithm == "lea"){

                     updateSelectInput(session, "goStatistic", choices = c("ks"="ks", "fisher"="fisher", "t"="t", "ks.ties"="ks.ties"), selected = "ks")
                }else {

                     updateSelectInput(session, "goStatistic", choices = c("fisher"="fisher"), selected = "fisher")
                }
            })


            observe({
                if(input$goGeneAlgorithm == "classic" || input$goGeneAlgorithm == "elim" || input$goGeneAlgorithm == "weight01" || input$goGeneAlgorithm == "lea"){

                     updateSelectInput(session, "goGeneStatistic", choices = c("ks"="ks", "fisher"="fisher", "t"="t", "ks.ties"="ks.ties"), selected = "ks")
                }else {

                     updateSelectInput(session, "goGeneStatistic", choices = c("fisher"="fisher"), selected = "fisher")
                }
            })

        geneOntology <- reactive({

            if(input$runRNAGO && input$miRNAorGene == 1){  
                        mir = input$showResults
                        predictions = getPredictedTargets(mir, species = 'hsa', method = 'geom', min_src = 2)
                        rankedGenes = predictions[,'rank_product']
                        selection = function(x) TRUE 
                        allGO2genes = annFUN.org(whichOnto=input$ontology, feasibleGenes = NULL, mapping="org.Hs.eg.db", ID = "entrez")
                        GOdata =  new('topGOdata', ontology = input$ontology, allGenes = rankedGenes, 
                                                annot = annFUN.GO2genes, GO2genes = allGO2genes, 
                                                geneSel = selection, nodeSize=10)
                        results = runTest(GOdata, algorithm = input$goAlgorithm, statistic = input$goStatistic)
                        allRes = GenTable(GOdata, statistic = results, orderBy = "statistic", topNodes = input$topRNAs)
                        allRes = allRes[,c('GO.ID','Term','statistic')]
                        colnames(allRes)[1] = "GO ID"
                        colnames(allRes)[3] = input$goStatistic
                        return(allRes)


              }
                if(input$runGeneGO && input$miRNAorGene == 2){  
                    data = dataM()
                    condition = as.factor(dataC()[,1])
                    design = model.matrix(~condition)

                    dge = DGEList(counts=as.matrix(data), group = condition)
                    dge = calcNormFactors(dge, method = input$normMeth) #webtool'da var. RLE dediği deseq, TMM dediği TMM, none dediği none
                    v = voom(dge,design, plot=F) 
                    fit = lmFit(v, design)
                    fit = eBayes(fit)
                    res = topTable(fit, coef=ncol(design), number = dim(data)[1])
                    geneList = res$adj.P.Val #voomNSC tarafından seçilen genler
                    names(geneList) = rownames(res) #voomNSC tarafından seçilen genler
                    biomarkers = selectedGenes #voomNSC tarafından seçilen genler
   
                    truefalse = function(allScore) {  
                      truefalse = is.element(names(geneList), biomarkers)
                      return(truefalse)
                    }

                    allGO2genes = annFUN.org(whichOnto=input$ontologyGene, feasibleGenes = NULL,
                                             mapping="org.Hs.eg.db", ID = "symbol")

                    GOdata =  new('topGOdata', ontology = input$ontologyGene, allGenes = geneList, 
                                  annot = annFUN.GO2genes, GO2genes = allGO2genes, 
                                  geneSel = truefalse, nodeSize=10)

                    results = runTest(GOdata, algorithm = input$goGeneAlgorithm, statistic = input$goGeneStatistic)

                    allRes = GenTable(GOdata, statistic = results, orderBy = "statistic", topNodes = input$topGenes)
                    allRes = allRes[,c('GO.ID','Term','statistic')]
                    colnames(allRes)[1] = "GO ID"
                    colnames(allRes)[3] = input$goGeneStatistic
                    return(allRes)

              }
                else{
                    return(allRes = NULL)
                }
            
        })



        geneOntologyPlot <- reactive({

            if(input$runVoomDDA){
                if(input$miRNAorGene == 1){  
                        mir = input$showResults
                        predictions = getPredictedTargets(mir, species = 'hsa', method = 'geom', min_src = 2)
                        rankedGenes = predictions[,'rank_product']
                        selection = function(x) TRUE 
                        allGO2genes = annFUN.org(whichOnto=input$ontology, feasibleGenes = NULL, mapping="org.Hs.eg.db", ID = "entrez")
                        GOdata =  new('topGOdata', ontology = input$ontology, allGenes = rankedGenes, 
                                                annot = annFUN.GO2genes, GO2genes = allGO2genes, 
                                                geneSel = selection, nodeSize=10)
                        results = runTest(GOdata, algorithm = input$goAlgorithm, statistic = input$goStatistic)
                        showSigOfNodes(GOdata, score(results), firstSigNodes = 5, useInfo = 'all')
                       


              }
                if(input$miRNAorGene == 2){  
                    data = dataM()
                    condition = as.factor(dataC()[,1])
                    design = model.matrix(~condition)

                    dge = DGEList(counts=as.matrix(data), group = condition)
                    dge = calcNormFactors(dge, method = input$normMeth) #webtool'da var. RLE dediği deseq, TMM dediği TMM, none dediği none
                    v = voom(dge,design, plot=F) 
                    fit = lmFit(v, design)
                    fit = eBayes(fit)
                    res = topTable(fit, coef=ncol(design), number = dim(data)[1])
                    geneList = res$adj.P.Val #voomNSC tarafından seçilen genler
                    names(geneList) = rownames(res) #voomNSC tarafından seçilen genler
                    biomarkers = selectedGenes #voomNSC tarafından seçilen genler
   
                    truefalse = function(allScore) {  
                      truefalse = is.element(names(geneList), biomarkers)
                      return(truefalse)
                    }

                    allGO2genes = annFUN.org(whichOnto=input$ontologyGene, feasibleGenes = NULL,
                                             mapping="org.Hs.eg.db", ID = "symbol")

                    GOdata =  new('topGOdata', ontology = input$ontologyGene, allGenes = geneList, 
                                  annot = annFUN.GO2genes, GO2genes = allGO2genes, 
                                  geneSel = truefalse, nodeSize=10)

                    results = runTest(GOdata, algorithm = input$goGeneAlgorithm, statistic = input$goGeneStatistic)
                    showSigOfNodes(GOdata, score(results), firstSigNodes = 5, useInfo = 'all')

              }

          }
            
        })


        output$geneOntologyTable <- DT::renderDataTable({

            if(input$runRNAGO || input$runGeneGO){


                  datatable(geneOntology(), extensions = c('Buttons','KeyTable', 'Responsive'), options = list(
                  dom = 'Bfrtip',
                  buttons = c('copy', 'csv', 'excel', 'pdf', 'print'), keys = FALSE, pageLength = 20
                  ))


            }else{return(NULL)}

        }) 

        #output$geneOntologyPlot <- renderPlot({

        #if(input$includePlot){
        #    geneOntologyPlot()         
        #       
        #   }
        #})


        output$downloadGoPlot <- downloadHandler(
    
            filename <- function() { paste('GOplot.pdf') },
            content <- function(file) {
                pdf(file)
                if(input$runVoomDDA == 0){stop('First, run voomDDA')}
                else{
                    geneOntologyPlot()
            }
                dev.off()
            },
            contentType = 'application/pdf'
        )
        
        }   
    })

})





