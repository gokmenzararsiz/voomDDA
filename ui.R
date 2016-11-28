shinyUI(fluidPage(


	titlePanel(title="voomDDA: Discovery of Diagnostic Biomarkers and Classification of RNA-Seq Data  (ver. 1.5)", windowTitle = "voomDDA for RNA-Seq"),

	sidebarPanel(width = 3,
		conditionalPanel(condition="input.tabs1=='Introduction'",
            HTML('<p align="center"><img src="voomdda.jpg" width=200 height=200></p>'),
            tags$head(includeScript("google-analytics.js"))

        ),

        conditionalPanel(condition="input.tabs1=='Manual'",
            HTML('<p align="center"><img src="manual.png" width=200 height=200></p>')

        ),

        conditionalPanel(condition="input.tabs1=='News'",
            HTML('<p align="center"><img src="news.png" width=200 height=200></p>')

        ),

        conditionalPanel(condition="input.tabs1=='Citation'",
            HTML('<p align="center"><img src="cite.png" width=200 height=200></p>')
        ),

        conditionalPanel(condition="input.tabs1=='voomDDA'",
            h4("1. Pre-processing"),
            h5("a) Filtering"),
            fluidRow(column(1),
                     column(10,
                            checkboxInput(inputId = "nearZeroF", label = "Near-zero variance filtering", value = TRUE),
                            checkboxInput(inputId = "varF", label = "Variance filtering", value = FALSE))
            ),
            HTML('<br>'),
            conditionalPanel(condition = "input.varF",
                textInput(inputId = "maxVar", label = "Number of genes with maximum variance", value = "")
            ),
            
            h5("b) Normalization"),
            fluidRow(column(1),
                     column(10,
                            radioButtons(inputId = "normMeth", label = "", choices = list("None" = "none", "DESeq median ratio" = "deseq", "TMM" = "TMM"), 
                                         selected = "none"))
            ),
            #HTML('<br>'),
            h4("2. Model building"),
            fluidRow(column(1),
                     column(11,
                            radioButtons(inputId = "models", label = "", selected = "vNSC",
                                               choices = c("voomNSC" = "vNSC",
                                                            "voomDLDA" = "vDLDA",
                                                            "voomDQDA" = "vDQDA"))
            )),
            HTML('<br>'),
            checkboxInput(inputId = "advanced", label = "Advanced options", value = FALSE),

            actionButton("runVoomDDA", "Run"),
            fluidRow(
                column(1),
                column(11, 
                       conditionalPanel(condition = "input.advanced",
                                 checkboxGroupInput(inputId = "advOpts", label = "", selected = c("trSummary","selGenes","selHeat"),
                                                    choices = c("Display training summary" = "trSummary",
                                                                "Display selected genes (voomNSC)" = "selGenes"))
                ))    
            )
            
        ),


        conditionalPanel(condition="input.tabs1=='Heatmap'",

            h4("Heatmap options"),
            actionButton(inputId = "runHeatMap", label = "Create heatmap"),
            checkboxInput(inputId = "centerHeat", label = "Center heatmap data", value = TRUE),
            checkboxInput(inputId = "darkTheme", label = "Dark theme", value = FALSE),
            selectInput("colorsHM", "Select a color scheme", choices = c("Blues" = "Blues", "Reds" = "Reds", "Greens" = "Greens"),
                              selected = "Blues")




            #sliderInput(inputId = "widthHM" , label = "Width", min = 200, max = 1000, value = 600, step = 50),
            #sliderInput(inputId = "heightHM" , label = "Height", min = 200, max = 1000, value = 600, step = 50)


        ),

        conditionalPanel(condition="input.tabs1=='Network'",

            actionButton(inputId = "runNetwork", label = "Create network plot")

        ),


        conditionalPanel(condition="input.tabs1=='GO'",

            radioButtons(inputId = "miRNAorGene", "", list("miRNA" = 1, "Gene" = 2), selected=1, inline = TRUE),

                conditionalPanel(condition="input.miRNAorGene=='1'",

                    numericInput(inputId = "topRNAs", "Select top genes", value = 20),  
                    
                    selectInput(inputId="goAlgorithm", "Select an algorithm", choices = c("classic"="classic" , "elim"="elim", "weight"="weight", "weight01"="weight01", "lea"="lea", "parentchild"="parentchild"), 
                                      selected = "classic"),

                    selectInput(inputId="goStatistic", "Select a statistic", choices = c("ks"="ks", "fisher"="fisher", "t"="t", "globaltest"="globaltest", "sum"="sum", "ks.ties"="ks.ties"), 
                                      selected = "ks"),

                    selectInput(inputId="showResults", "Show results for", choices = ""),

                    selectInput(inputId="ontology", "Ontology", choices = c("Biological process" = "BP", "Molecular function" = "MF", "Cellular components" = "CC"), 
                                      selected = "BP"),
                    
                    actionButton(inputId = "runRNAGO", label = "GO get it!")
              ),  


                conditionalPanel(condition="input.miRNAorGene=='2'",

                    numericInput(inputId = "topGenes", "Select top genes", value = 20),   

                    selectInput(inputId="goGeneAlgorithm", "Select an algorithm", choices = c("classic"="classic" , "elim"="elim", "weight"="weight", "weight01"="weight01", "lea"="lea", "parentchild"="parentchild"), 
                                      selected = "classic"),

                    selectInput(inputId="goGeneStatistic", "Select a statistic", choices = c("ks"="ks", "fisher"="fisher", "t"="t", "globaltest"="globaltest", "sum"="sum", "ks.ties"="ks.ties"), 
                                      selected = "ks"),

                    selectInput(inputId="ontologyGene", "Ontology", choices = c("Biological process" = "BP", "Molecular function" = "MF", "Cellular components" = "CC"), 
                                      selected = "BP"),
                    
                    actionButton(inputId = "runGeneGO", label = "GO get it!")
              ),
              br(),
              downloadButton("downloadGoPlot", "Download GO Plot")        
        ),
            


        conditionalPanel(condition="input.tabs1=='Authors & News'",
			HTML('<p align="center"> <a href="https://www.hacettepe.edu.tr/english/" target="_blank"><img src="hulogo.JPEG" width=150 height=150></a> </p>')
        ),

		conditionalPanel(condition="input.tabs1=='Data upload'",
			h4("Input data"),
			radioButtons("dataInput", "", list("Load example data"=1,"Upload a file"=2), selected=1),

			conditionalPanel(condition="input.dataInput=='1'",
				h5("Load example data:"),
				
				radioButtons("sampleData", "", list("Cervical cancer (miRNA expression data)"=1, "Lung cancer (Gene expression data)"=2),selected=1)
				
                #HTML('<p>n: number of observations</p>'),
				#HTML('<p>p: number of genes</p>')
			),
			
			conditionalPanel(condition="input.dataInput=='2'",
				h5("Upload train set: "),
				fileInput("uploadTrain", "", multiple = FALSE),

				h5("Upload class labels: "),
				fileInput("uploadCond", "", multiple = FALSE),
				
				h5("Upload test set: "),
				fileInput("uploadTest", "", multiple = FALSE),
				
				radioButtons("fileSepDF", "Delimiter:", list("Comma"=1,"Tab"=2,"Semicolon"=3,"Space"=4),selected=2),  

				HTML('<p>You can upload your data as separated by comma, tab, semicolon or space.</p>'),
                HTML('<p>Note: Samples are in the column and genes are in the rows.</p>')
			)
		)
),


	mainPanel(
		tabsetPanel(
			tabPanel("Introduction",
			        HTML('<br>'),
			        HTML('<p>Identifying the relevant genes (or other genomic features such as transcripts, miRNAs, lncRNAs, etc.) across the conditions (e.g. tumor and non-tumor tissue samples) is a common research interest in gene-expression studies. In this gene selection, researchers are often interested in detecting a small set of genes for diagnostic purpose in medicine that involves identification of the minimal subset of genes that achieves maximal predictive performance. biomarker discovery and classification problem.</p>'),
                    HTML('<p>VoomDDA is a decision support tool developed for RNA-Sequencing datasets to assist researchers in their decisions for diagnostic biomarker discovery and classification problem. VoomDDA consists both sparse and non-sparse statistical learning classifiers adapted with voom method. Voom is a recent method that estimates the mean and variance relationship of the log-counts of RNA-Seq data (log counts per million, log-cpm) at observational level. It also provides precision weights for each observation that can be incorporated with the log-cpm values for further analysis. Algorithms in our tool incorporates the log-cpm values and the corresponding precision weights into biomarker discovery and classification problem. For this purpose, these algorithms use weighted statistics in estimating the discriminating functions of the used statistical learning algorithms.</p>'),
			        HTML('<p>VoomNSC is a sparse classifier that is developed to bring together two powerful methods for RNA-Seq classification:</p>'),
			         
                    h5("1.    to extend voom method for RNA-Seq classification studies,"),
			        h5("2.    to make nearest shrunken centroids (NSC) algorithm available for RNA-Seq technology."),
                    HTML('<p>VoomNSC both provides fast, accurate and sparser classification results for RNA-Seq data. More details can be found in the research paper. This tool also includes RNA-Seq extensions of diagonal linear and diagonal quadratic discriminant classifiers: (i) voomDLDA and (ii) voomDQDA.</p>'),

          
				h5("References"),
                HTML('<p> [1] Zararsiz, G., Goksuluk, G., Korkmaz, S., et al. (2015). VoomDDA: Discovery of Diagnostic Biomarkers and Classification of RNA-Seq Data.</p>'),

                HTML('<p> [2] Law, C.W., Chen, Y., Shi, W., et al. (2014). <a href="http://www.genomebiology.com/2014/15/2/R29">voom: Precision weights unlock linear model analysis tools for RNA-Seq read counts.</a> Genome Biology; 15:R29.</p>'),

				HTML('<p> [3] Tibshirani, R., Hastie, T., Narasimhan, B., et al. (2002). <a href="http://www.pnas.org/content/99/10/6567.abstract" target="_blank">Diagnosis of multiple cancer types by shrunken centroids of gene expression.</a> PNAS; 99(10): 6567-72. </p>'),

                HTML('<p> [4] Dudoit, S., Fridlyand, J. and Speed, T.P. (2002). <a href="http://amstat.tandfonline.com/doi/abs/10.1198/016214502753479248" target="_blank">Comparison of Discrimination Methods for the Classification of Tumors Using Gene Expression Data.</a> Journal of the American Statistical Association; 97(457): 77-87.</p>')
				  
			),

			tabPanel("Data upload",
                navbarPage(
                               title = '',
tabPanel('Train',                DT::dataTableOutput('RawDataTrain')),
                               tabPanel('Test',              DT::dataTableOutput('RawDataTest'))
                 )
            ),


            tabPanel("voomDDA",

                #downloadButton("downloadNormTest", "Download univariate results"),
                #downloadButton("downloadNormPlot", "Download univariate plots"),

                verbatimTextOutput("trainConsole"),
                HTML('<br>'),
                HTML('<br>'),
                tagList(
                    tags$head(
                    tags$link(rel="stylesheet", type="text/css",href="style.css"),
                    tags$script(type="text/javascript", src = "busy.js")
                  )
                ),

                div(class = "busy",
                    p("Getting voomDDA results... this may take a while...."),
                    img(src="loading.gif")
                )

                #plotOutput("heatMap")
            ),



            tabPanel("Heatmap",
		          HTML('<br>'),
		          HTML('<br>'),
              d3heatmap::d3heatmapOutput("heatMap", width = "100%", height = "600px"),
                tagList(
                    tags$head(
                    tags$link(rel="stylesheet", type="text/css",href="style.css"),
                    tags$script(type="text/javascript", src = "busy.js")
                  )
                ),

                div(class = "busy",
                    p("Creating heatmap... this may take a while...."),
                    img(src="loading.gif")
                )



            ),


            tabPanel("Network",

                plotOutput('newtwork'),
                tagList(
                    tags$head(
                    tags$link(rel="stylesheet", type="text/css",href="style.css"),
                    tags$script(type="text/javascript", src = "busy.js")
                  )
                ),

                div(class = "busy",
                    p("Creating network plot... this may take a while...."),
                    img(src="loading.gif")
                )


            ),

            tabPanel("GO",
                h4('Gene ontology results'),
                #verbatimTextOutput('selectedGenes'),
                DT::dataTableOutput('geneOntologyTable'),
                plotOutput('geneOntologyPlot'),
                tagList(
                    tags$head(
                    tags$link(rel="stylesheet", type="text/css",href="style.css"),
                    tags$script(type="text/javascript", src = "busy.js")
                  )
                ),

                div(class = "busy",
                    p("Getting ontology results... this may take a while...."),
                    img(src="loading.gif")
                )


                #verbatimTextOutput("geneResult")
            ),



            tabPanel("Manual",

            h5("Tutorial"),
            HTML('<p><b>1.Uploading the data</b></p>'),

            HTML('<p>Two example datasets are available in voomDDA web application. Cervical cancer is a miRNA, lung cancer is a gene expression dataset. For GO analysis, users should select the necessary option (miRNA or gene) to obtain the related analysis results.</p>'),
            HTML('<p>VoomDDA application requires three inputs from the user. Train and test sets should be text files (.txt) that contain the raw mapped read counts in a matrix form, where rows correspond to genomic features (for simplicity of language, let’s say genes) and the columns correspond to observations (or samples). This type of count data can be obtained from feature counting softwares such as HTSeq [1] or featureCounts [2]. Note that this type of count data should contain the raw number of mapped reads, should not be normalized or contain RPKM values. Class labels should also be in a text file (.txt) and should contain each sample condition. Note that each row should contain only one label of a sample. Example datasets for Witten et al. cervical dataset are given as below:</p>'),
            HTML('<p><a href="cervical_train.txt" download>Training set of cervical data</a></p>'),
            HTML('<p><a href="cervical_test.txt" download>Test set of cervical data</a></p>'),
            HTML('<p><a href="cervical_cond.txt" download>Class labels of cervical data</a></p>'),
            HTML('<p>If the purpose is the prediction of the class labels of new test observations, users should upload all three necessary files. However, test set is not required, when the purpose is just the identification of the diagnostic biomarkers.</p>'),
            HTML('<p>After uploading the data, make sure that the data is displayed in the screen.</p>'),
            HTML('<br>'),
            HTML('<p><b>2. Pre-processing the Data</b></p>'),
            HTML('<p><b>2.1. Filtering</b></p>'),
            HTML('<p>VoomDDA classifiers (VoomNSC, VoomDLDA and VoomDQDA) introduced in this application have the same assumptions with voom+limma pipeline [3], that is to filter out the rows with zero or very very low counts. In RNA-Seq data, we often meet with count data that contains rows with single unique values (mostly zero). This type of data may lead to unreliable estimation of the mean and variance relationship of the data and unstable model fitting for the introduced classifiers. Three possible filtering criteria are available: (i) DESeq2 outlier and independent filtering, (ii) near-zero variance filtering, (iii) variance filtering.</p>'),
            HTML('<p>DESeq2 package [4] contains a filtering criteria based on outlier detection and independent filtering.  Outliers are detected based on the Cook’s distance and independent filtering is applied based on the gene-wise mean normalized counts. More details can be obtained in the vignette of DESeq2 package [5].</p>'),
            HTML('<p>Near-zero variance filtering is described in caret package of R [6]. This package applies filtering based on two criteria: (i) the frequency of the most frequent value to the most frequent second value is higher than 19 (95/5), (ii) the number of unique values divided by the sample size is less than 10%.</p>'),
            HTML('<p>Variance filtering is another option to filter out the non-informative genes. This option may also be selected to decrease the computational cost of the model building process for very large datasets. After selecting this option, users can enter the number of genes desired to be included to the classification models.</p>'),
            HTML('<p>After selecting one or multiple filtering criteria, filtering statistics are demonstrated in the screen.</p>'),
            HTML('<p><b>2.2. Normalization</b></p>'),
            HTML('<p>Library sizes for each observation are dependent on the experimental design and may lead to the existence of technical biases. These biases can have significant effect on the classification results and should be corrected before starting to classification model building. In our experiments, we found that normalization has a significant effect on the classification results for datasets that have very large library size differences across samples. Two normalization approaches are available in the application: (i) DESeq median ratio [7], (ii) trimmed mean of M values (TMM) [8]. More details about this approaches can be found in referenced papers.</p>'),
            HTML('<br>'),
            HTML('<p><b>3. Model Building for Classification</b></p>'),
            HTML('<p>After data processing, users can build classification models with three introduced algorithms: (i) voomNSC, (ii) voomDLDA, (iii) voomDQDA. VoomNSC is a sparse classifier that brings together two powerful methods, voom method [3] and nearest shrunken centroids algorithm [9], for the classification of RNA-Seq data. VoomDLDA and voomDQDA are non-sparse classifiers which are the extensions of diagonal discriminant classifiers [10]. Details of these classifiers are given in the referenced paper [11].</p>'),
            HTML('<p>After selecting any of the three classifiers, a summary of the fitting process is displayed in the screen. A confusion matrix and several statistical diagnostic measures are given to examine how successful the classifier fit to the given data. Furthermore, a heatmap plot is constructed to display the expression levels of genes and the gene-wise and sample-wise relationships. Heatmap is displayed for the entire unfiltered genes for non-sparse classifiers, while displayed for the selected gene subset for sparse voomNSC classifier.</p>'),
            HTML('<br>'),
            HTML('<p><b>4. Identification of Diagnostic Biomarkers</b></p>'),
            HTML('<p>If VoomNSC is the selected classifier, the subset of genes, that are most relevant with the class condition, are identified and the gene names are displayed in the screen. Several plots are also given. First plot demonstrates the selection of the threshold parameter. The parameter which fits the most accurate and sparsest model is identified as optimal. Second plot displays the distribution of selected genes in each class. Third plot displays the shrunken differences of the selected genes. Final plot is the heatmap plot discussed in the previous section.</p>'),
            HTML('<br>'),
            HTML('<p><b>5. Prediction</b></p>'),
            HTML('<p>Based on the selected classifier, predictions appear on the screen for each test observation. Note that the test observations should be processed as same as the training observations. Same experimental and computational procedures should be applied before obtaining the raw count data. Data should be in the same format as the training data to obtain the predictions. It should contain the raw mapped read counts, and the gene names should match with the training data.</p>'),
            HTML('<p>VoomDDA application filters and normalizes the test data based on the information obtained from the training data. Thus, the estimated parameters from the training data are used for the test data. This guarantees that both sets are on the same scale and homoscedastic each other.</p>'),
            HTML('<br>'),
            HTML('<p><b>6. Downstream Analysis</b></p>'),
            HTML('<p>After detecting diagnostic biomarkers via voomNSC algorithm, it may be useful to visualize the results to see the interactions or go further analysis, such as GO analysis. For this purpose, several downstream analysis tools are also available in this web application. These tools include heatmaps, network analysis and gene ontology analysis. Detailed information about gene ontology analysis can be found in topGO BIOCONDUCTOR package.</p>'),
            HTML('<br>'),
            HTML('<p><b>References</b></p>'),
            HTML('<p>[1] Anders, S., Pyl, P.T., and Huber, W. (2015) HTSeq - a Python framework to work with high-throughput sequencing data. Bioinformatics; 31(2):166-9.</p>'),
            HTML('<p>[2] Liao, Y., Smyth, G.K., and Shi, W. (2013). featureCounts: an efficient general-purpose program for assigning sequence reads to genomic features. Bioinformatics. doi: 10.1093/bioinformatics/btt656.</p>'),
            HTML('<p>[3] Law, C.W., Chen, Y., Shi, W. and Smyth, G.K. (2014). voom: Precision weights unlock linear model analysis tools for RNA-Seq read counts. Genome Biology; 15:R29.</p>'),
            HTML('<p>[4] Love, M.I., Huber, W. and Anders, S. (2015). Moderated estimation of fold change and dispersion for RNA-seq data with DESeq2. Genome Biology; 15(550).  doi:10.1186/s13059-014-0550-8 .</p>'),
            HTML('<p>[5] Love, M.I., Huber, W. and Anders, S. (2015). Differential analysis of count data – the DESeq2 package. http://www.bioconductor.org/packages/release/bioc/vignettes/DESeq2/inst/doc/DESeq2.pdf (19.06.2015).</p>'),
            HTML('<p>[6] Kuhn, M. (2008). Building Predictive Models in R Using the caret Package. Journal of Statistical Software; 28(5).</p>'),
            HTML('<p>[7] Anders, S. and Huber, W. (2010). Differential expression analysis for sequence count data. Genome Biology; 11(R106): doi:10.1186/gb-2010-11-10-r106 .</p>'),
            HTML('<p>[8] Robinson, M.D., and Oshlack, A. (2010). A scaling normalization method for differential expression analysis of RNA-seq data. Genome Biology; 11(R25).</p>'),
            HTML('<p>[9] Tibshirani, R., Hastie, T., Narasimhan, B. and Chu, G. (2002). Diagnosis of multiple cancer types by shrunken centroids of gene expression. PNAS; 99(10): 6567–72.</p>'),
            HTML('<p>[10] Dudoit, S., Fridlyand, J. and Speed, T.P. (2002). Comparison of Discrimination Methods for the Classification of Tumors Using Gene Expression Data. Journal of the American Statistical Association; 97(457): 77-87.</p>'),
            HTML('<p>[11] Zararsiz, G., Goksuluk, D, Korkmaz, S., et al. (2015). VoomDDA: Discovery of Diagnostic Biomarkers and Classification of RNA-Seq Data.</p>')
            ),

			tabPanel("Authors & News",
                h4("Authors"),
			    HTML('<p><a href="http://bit.do/gokmenzararsiz" target="_blank"> <b>Gokmen Zararsiz, PhD</b></a><p>'),
			    HTML('<p>Hacettepe University Faculty of Medicine <a href="http://www.biostatistics.hacettepe.edu.tr" target="_blank"> Department of Biostatistics</a><p>'),
			    HTML('<p><a href="mailto:gokmen.zararsiz@hacettepe.edu.tr" target="_blank">gokmen.zararsiz@hacettepe.edu.tr</a><p>'),

                br(),

			    HTML('<p><a href="http://www.biostatistics.hacettepe.edu.tr/cv/Dincer_Goksuluk_CV_Eng.pdf" target="_blank"> <b>Dincer Goksuluk, PhD</b></a><p>'),
			    HTML('<p>Hacettepe University Faculty of Medicine <a href="http://www.biostatistics.hacettepe.edu.tr" target="_blank"> Department of Biostatistics</a><p>'),
			    HTML('<p><a href="mailto:dincer.goksuluk@hacettepe.edu.tr" target="_blank">dincer.goksuluk@hacettepe.edu.tr</a><p>'),
                br(),

                HTML('<p><a href="http://www-huber.embl.de/users/klaus/" target="_blank"> <b>Bernd Klaus, PhD</b></a><p>'),
                HTML('<p>EMBL Heidelberg<p>'),
                HTML('<p><a href="bernd.klaus@embl.de" target="_blank">bernd.klaus@embl.de</a><p>'),
                br(),

			    HTML('<p><a href="http://yunus.hacettepe.edu.tr/~selcuk.korkmaz/" target="_blank"> <b>Selcuk Korkmaz, PhD</b></a><p>'),
                HTML('<p>Hacettepe University Faculty of Medicine <a href="http://www.biostatistics.hacettepe.edu.tr" target="_blank"> Department of Biostatistics</a><p>'),
                HTML('<p><a href="mailto:selcuk.korkmaz@hacettepe.edu.tr" target="_blank">selcuk.korkmaz@hacettepe.edu.tr</a><p>'),
                br(),

			    HTML('<p><a href="http://aves.istanbul.edu.tr/vahap.eldem/" target="_blank"> <b>Vahap Eldem, PhD</b></a><p>'),
			    HTML('<p>Istanbul University Faculty of Science <a href="http://fen.istanbul.edu.tr/biyoloji/#" target="_blank"> Department of Biology</a><p>'),
			    HTML('<p><a href="mailto:vahap.eldem@istanbul.edu.tr" target="_blank">vahap.eldem@istanbul.edu.tr</a><p>'),
                br(),

			    HTML('<p><a href="http://unverlab.karatekin.edu.tr" target="_blank"> <b>Turgay Unver, PhD</b></a><p>'),
			    HTML('<p>Cankiri Karatekin University Faculty of Science <a href="http://biyoloji.karatekin.edu.tr/default.aspx" target="_blank"> Department of Biology</a><p>'),
			    HTML('<p><a href="mailto:turgayunver@gmail.com" target="_blank">turgayunver@gmail.com</a><p>'),
                br(),

			    HTML('<p><a href="http://www.biyoistatistik.hacettepe.edu.tr/erdem.html" target="_blank"> <b>Erdem Karabulut, PhD</b></a><p>'),
			    HTML('<p>Hacettepe University Faculty of Medicine <a href="http://www.biostatistics.hacettepe.edu.tr" target="_blank"> Department of Biostatistics</a><p>'),
			    HTML('<p><a href="mailto:ekarabul@hacettepe.edu.tr" target="_blank">ekarabul@hacettepe.edu.tr</a><p>'),
                br(),

			    HTML('<p><a href="http://aves.erciyes.edu.tr/ahmetozturk/" target="_blank"> <b>Ahmet Ozturk, PhD</b></a><p>'),
			    HTML('<p>Erciyes University Faculty of Medicine <a href="http://biyoistatistik.erciyes.edu.tr" target="_blank"> Department of Biostatistics</a><p>'),
			    HTML('<p><a href="mailto:ahmets67@hotmail.com" target="_blank">ahmets67@hotmail.com</a><p>'),


                HTML('<br>'),


                h4("News"),
                h5("Version 1.5 (November 25, 2016)"),
                HTML('<p>(2) Lung cancer data added as an example dataset <p>'),
                HTML('<p>(2) Bug fixes and improvements <p>'),

                h5("Version 1.4 (November 14, 2016)"),
                HTML('<p>(1) Bug fixes and improvements <p>'),

                h5("Version 1.3 (November 5, 2016)"),
                HTML('<p>(1) Gene ontology results added <p>'),
                HTML('<p>(2) Gene ontology plot added <p>'),
                HTML('<p>(3) Bug fixes and improvements <p>'),

                h5("Version 1.2 (August 20, 2016)"),
                HTML('<p>(1) Heatmap added <p>'),
                HTML('<p>(2) Network plot added <p>'),
                HTML('<p>(3) Bug fixes and improvements <p>'),

                h5("Version 1.1 (July 18, 2016)"),
                HTML('<p>(1) Upgraded to shiny version 0.14 <p>'),

                h5("Version 1.0 (June 18, 2015)"),
                HTML('<p>(1) VoomDDA web application has been released. <p>'),








                HTML('<br>'),

                h5("Other Tools"),

                HTML('<p><a href="http://www.bioconductor.org/packages/release/bioc/html/MLSeq.html" target="_blank"> <b>MLSeq: Machine learning interface for RNA-Seq data </b></a><p>'),
                HTML('<p><a href="http://www.biosoft.hacettepe.edu.tr/MLViS/" target="_blank"> <b>MLViS: machine learning-based virtual screening tool </b></a><p>'),
                HTML('<p><a href="http://www.biosoft.hacettepe.edu.tr/easyROC" target="_blank"> <b>easyROC: a web-tool for ROC curve analysis </b></a><p>'),
                HTML('<p><a href="http://www.biosoft.hacettepe.edu.tr/MVN" target="_blank"> <b>MVN: a web-tool for assessing multivariate normality </b></a><p>'),
                HTML('<p><a href="http://www.biosoft.hacettepe.edu.tr/DDNAA/" target="_blank"> <b>DDNAA: Decision support system for differential diagnosis of nontraumatic acute abdomen </b></a><p>'),
                HTML('<br>'),



                h6("Please feel free to send us bugs and feature requests.")

            ),

            tabPanel("Citation",
             HTML('<br>'),
             HTML('<p>If you use this application, please cite it as below:</p>'),
             HTML('<p><b>Zararsiz, G., Goksuluk, D, Korkmaz, S., et al. (2015). VoomDDA: Discovery of Diagnostic Biomarkers and Classification of RNA-Seq Data. Available at: http://www.biosoft.hacettepe.edu.tr/voomDDA .</b></p>')
            ),

            id="tabs1", type = "pills"
		),

        tags$head(tags$style(type="text/css", "label.radio { display: inline-block; }", ".radio input[type=\"radio\"] { float: none; }"),
        tags$style(type="text/css", "select { max-width: 200px; }"),
        tags$style(type="text/css", "textarea { max-width: 185px; }"),
        tags$style(type="text/css", ".jslider { max-width: 200px; }"),
        tags$style(type='text/css', ".well { max-width: 330px; }"),
        tags$style(type='text/css', ".span4 { max-width: 330px; }"))

     #   tags$head(
     #   tags$link(rel = "shortcut icon", href = "favicon-2.ico"))

 )
))




