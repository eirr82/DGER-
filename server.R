#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#
#DGR Ebayes Analyzer v1.1
#Developed by Javier Bertol Chorro


##Global options
options(spinner.color="#FF8C00")
options(shiny.maxRequestSize = 70*1024^2)
options(expressions = 5000)
options(repos = BiocInstaller::biocinstallRepos())
getOption("repos")
 

#libraries 

library(shiny)
library(shinydashboard)
library(pd.clariom.s.human)
library(limma)
library(annotate)
library(clariomshumantranscriptcluster.db)
library(GOfuncR)
library(knitr)
library(oligo)
library(DT)
library(Heatplus)
library(shinycssloaders)
library(Homo.sapiens)
library(affycoretools)
library(org.Hs.eg.db)
library(ReactomePA)

#header

  
header <- dashboardHeader(
  title = "DGE Ebayes Analyzer"
)

#sidebar

sidebar <- dashboardSidebar(
  #sidebarmenu
  sidebarMenu(id ="tabs",
              menuItem("Introduction", tabName = "intro", icon = icon("edit") ),
              menuItem("Data Input", tabName = "data", icon = icon("server")),
              menuItem("Quality Contol", tabName = "qc", icon = icon("braille"),
                       menuSubItem("Before Normalization", tabName = "qcbn"),
                       menuSubItem("After Normalization", tabName = "qcan" )
              ),
              menuItem("Filtering", tabName = "filter", icon = icon("filter")),
              menuItem("DGE", tabName = "DGE", icon = icon("hourglass"),
                       menuSubItem("One Group", tabName = "DGE1"),
                       menuSubItem("Multiple Groups", tabName = "DGE2")
              ),
              menuItem("Results", tabName = "res", icon = icon("eye"),
                       menuSubItem("One Group", tabName = "res1"),
                       menuSubItem("Multiple Groups", tabName = "res2")
              ),
              menuItem("GO Analysis", tabName = "topgo", icon = icon("spinner")
                       ),
              menuItem("Gsea Enrichemnt", tabName = "gsea")
              
  )
)

#body

body <- dashboardBody(
  
  #tabitems
  tabItems(
    
    ##################Intro tab#################
    tabItem(tabName = "intro",
            h1("Introduction"),
            br(),
            p("DGE Ebayes Analyzer is an online app to perform an analysis to a 
              differential expression genes microarrays experiment"),
            p("This app uses the", tags$a(href="http://www.statsci.org/smyth/pubs/ebayes.pdf",
                                          "Empirical Bayes Method"), "from", 
              tags$a(href="https://bioconductor.org/help/course-materials/2005/BioC2005/labs/lab01/printing/6.appendix.pdf",
                     "Limma"),"(Linear Models for Microarray Data), R package to calculate the expression of the genes"),
            br(),
            
            h2("How to use it?"),
            p("This section provides a short guide to run the expression experiment analysis"),
            h3("Data Input section"),
            p("The Data input section, provides all the data tha the app need tu run the analysis"),
            p("In the CEL FILES section, the user have to provide the path to his experiment files.(.cel),
              the directory must contain only CEL files, otherwise the user will get and error"),
            p("In the Targets section, the user have to provide the customs targets file in csv format"),
            p("The custom targets file must be as follows:"),
            tags$ol(
              tags$li("FileNames column with the CEL files name"),
              tags$li("Group column with the group of each sample of the experiment"),
              tags$li("Name column with an alias for CEL files(for better results understanding"),
              tags$li("PlotColor column indicating wich color the user wants to represent the groups 
                      for plotting purposes")
              
              ),
            p(tags$a(href="http://www.stat.columbia.edu/~tzheng/files/Rcolor.pdf", "Here"),
              "is a guide for all colors R can accept"),
            p("The input data section ,after the user provides the files and use submit button, shows
              the CEL files loaded, the groups of the experiment, and the targets custom file info"),
            h3("Quality Control"),
            p("The quality control section runs a quality control on the microarray data, on the raw data,
              and the normalized data:"),
            p("Normalization of the samples is done with",
              tags$a(href="https://www.ncbi.nlm.nih.gov/pubmed/12925520)","RMA"),"method"),
            p("The quality control includes the following plots:"),
            tags$ol(
              tags$li("Histogram of arrays intensities(only for rawdata)"),
              tags$li("Boxplot"),
              tags$li("Dendogram"),
              tags$li("PCA Analysis")
            ),
            h3("Filtering section(NOT IMPLEMENTED YET)"),
            h3("DGE section"),
            p("DGE section is divided in two sections, one for one group comparion, and another for an 
              experiments with more than a group to compare(max 6 groups for three comparations)"),
            p("In the DGE section in the tab 'select groups'the user must indicate on wich groups
              run the analysis, after loading the groups the app shows the matrix and contrasts designed
              to perform the ebayes method"),
            h3("Results section"),
            p("This section shows the results of the DGE analysis, user must navigate to the proper tab(one or
              multiple groups)"),
            p("The One Group Results Tab includes the following outputs:"),
            tags$ol(
              tags$li("Ebayes results in table format with the expression of the genes"),
              tags$li("Selected Genes,the table with the selected genes to perform the GO analysis.
                      The User can control the Adjust.P.Value to filter results and the number of genes"),
              tags$li("Volcano plot"),
              tags$li("Heatmap")
              ),
            p("The multiple groups tab includes the following outputs:"),
            tags$ol(
              tags$li("Ebayes results in table format. The usar can control wich of the comparations see"),
              tags$li("Selected Genes,the table with the selected genes to perform the GO analysis.
                      The User can control the Adjust.P.Value to filter results and the number of genes"),
              tags$li("Volcano plot"),
              tags$li("Heatmap"),
              tags$li("Summary of decide tests function to select the up and down regulated genes"),
              tags$li("Venn Diagram with the result of decide tests"),
              tags$li("A list with the names of the genes selected by decision tests")
              ),
            h3("GO Analysis(NOT IMPLEMENTED YET)")
            
            
            
            
            ),
    ################End Intro tab##############
    
    ################Inputdata tab##############
    tabItem(tabName = "data",
            fluidRow(
              box(title = "Input Data",
                  solidHeader = T, status = "warning", width = 12,
                  
                  fluidRow(
                    ###Input CEL files###
                    box(title ="CEL FILES",
                        solidHeader = T, status = "info", width = 6,
                        
                        fileInput("file", "Load CEL files(.CEL)",
                                  multiple = TRUE,
                                  accept = ".CEL",
                                  placeholder = "Please, insert the datapath of  CEL files"),
                        submitButton("Submit")
                    ),
                    
                    ###Input target files###
                    box(title = "TARGETS",
                        solidHeader = T, status = "info", width = 6,
                        
                        fileInput("file2", "Load custom targets file",
                                  accept = ".csv",
                                  placeholder = "Please, load custom targets file"),
                        submitButton("Submit")
                    )
                  )
                  
              )
            ),
            
            #####Data summary########
            
            fluidRow(
              box(title = "Input data Summary",
                  solidHeader = T, status = "warning", width = 12,
                  
                  fluidRow(
                    box(title = "Samples Loaded",
                        solidHeader = T, status = "info", width = 6,
                        h4("The following samples are loaded"),
                        verbatimTextOutput("expression")),
                    
                    box(title = "Groups Levels",
                        h4("These are the following group Levels"),
                        solidHeader = T, status = "info", "width" =6,
                        verbatimTextOutput("groups"))
                  )
                  
                  
              )
              
            ),
            
            fluidRow(
              box(title = "Targets info",
                  solidHeader = T, status = "warning", width = 12,
                  
                  fluidRow(
                    box(title = "Targets info",
                        solidHeader = T, status = "info",width = 12,
                        withSpinner(dataTableOutput("targets")))
                    
                  ))
            ),
            actionButton("btn_switchtab", "Go to QC", class ="btn-info")
    ),
    
    ###########End datainput tab#########
    
    ###########QC tab################
    
    ############QC before normalization###########
    tabItem(tabName = "qcbn",
            h2("Quality Control Before Normalization"),
            fluidRow(
              box(title = "Samples intensity Distribution",
                  solidHeader = T, status = "warning", width = 12,
                  withSpinner(plotOutput("histogram")),
                  downloadButton("plot1", label = "Download")
              ),
              
              fluidRow(
                box(title = "Boxplot",
                    solidHeader = T, status = "info", width = 6,
                    withSpinner(plotOutput("boxplotb")),
                    downloadButton("plot2", label = "Download")),
                box(title = "Dendogram",
                    solidHeader = T, status = "info", width = 6,
                    withSpinner(plotOutput("dendogramb")),
                    downloadButton("plot3", label = "Download")),
                box(title= "PCA Analysis", 
                    solidHeader = T, status = "info", width = 12,
                    withSpinner(plotOutput("PCAB")),
                    downloadButton("plot4", label = "Download"))
              )
            )),
    ########End QC before normalization tab########
    
    ########QC after normalization tab#############
    tabItem(tabName = "qcan",
            h2("Quality Control After Normalization"),
            fluidRow(
              box(title = "Boxplot",
                  solidHeader = T, status = "info", width = 6,
                  withSpinner(plotOutput("boxplotb2")),
                  downloadButton("plot5", label = "Download")),
              box(title = "Dendogram",
                  solidHeader = T, status = "info", width = 6,
                  withSpinner(plotOutput("dendogramb2")),
                  downloadButton("plot6", label = "Download")),
              box(title= "PCA Analysis", 
                  solidHeader = T, status = "info", width = 12,
                  withSpinner(plotOutput("PCAB2")),
                  downloadButton("plot7", label = "Download"))
              
            ),
            h4("Press the button to dowload normalized expression data
               for furhter analysis"),
            downloadButton("datanorm", "Download", class = "butt1"),
            tags$head(tags$style(".butt1{background-color:orange;} 
                                  .butt1{color: black;} 
                                  .butt1{font-family: Courier New}"))
           
    ),
    
    #####End QC after normalization tab######
    
    #####Filter data tab(not implemented)#####
    tabItem(tabName = "filter",
            h2("Filter the data(NOT IMPLEMENTED YET"),
            fluidRow(
              box(title = "Filtering",
                  solidHeader = T, status = "warning", width = 12,
                  fluidRow(
                    box(title = "Filter Results",
                        solidHeader = T, status = "info", width = 6,
                        textOutput("filter")),
                    box(title = "Select Filter Options",
                        solidHeader = T, status = "info", width = 6,
                        #input cutoff
                        sliderInput("cutoff", label = h3("Select Cut-Off:"),
                                    min = 0, max = 1,
                                    value = 0.10),
                        #Require Entrez
                        radioButtons("ReqEntrez","Require Entrez",
                                     choices = c(Yes = TRUE,
                                                 No = FALSE),
                                     selected = TRUE),
                        #Remove Duplicates
                        radioButtons("RemDup", "Remove Duplicates",
                                     choices = c(Yes = TRUE, 
                                                 No = FALSE),
                                     selected = TRUE),
                        submitButton("Submit")
                    )
                  ))
            )),
    
    #####End filtering tab##############################
    
    #####DGE one group tab #################################
    tabItem(tabName = "DGE1",
            h2("DGE Analysis"),
            fluidRow(
              box(title = "Selecting Groups",
                  solidHeader = T, status = "warning", width = 12,
                  box(title = "Select Groups",
                      solidHeader = T, status = "info", width = 12,
                      textOutput("groupsava"),
                      tags$style("#groupsava{ font-weight: bold; color: #0033FF; 
                                 padding-top: .3cm; padding-bottom: .3cm;}"),
                      
                      fluidRow(
                        column(6,
                               textInput(inputId = "gr1",
                                         label="Group 1",
                                         value = "")
                        ),
                        column(6,
                               textInput(inputId = "gr2",
                                         label = "Group 2",
                                         value = "")
                        )),
                      
                      submitButton("Submit")
                      
                      
                      ))
            ),
            fluidRow(
              box(title = "Experiment Design",
                  solidHeader = T, status = "warning", width = 12,
                  fluidRow(
                    box(title = "Matrix Design",
                        solidHeader = T, status = "info", width = 6,
                        dataTableOutput("design")),
                    box(title = "Contrasts Experiment",
                        solidHeader = T, status = "info", width = 6,
                        dataTableOutput("cont")
                        
                        
                    )
                    
                  )
              )
              
            )),
    
    ######End DGE one group tab####################
    
    ######DGE more groups tab#####################
    tabItem("DGE2",
            fluidRow(
              box(title = "Selecting Groups",
                  solidHeader = T, status = "warning", width = 12,
                  fluidRow(
                    box(title = "Select Groups",
                        solidHeader = T, status = "info", width = 12,
                        textOutput("groupsava2"),
                        tags$style("#groupsava2{ font-weight: bold; color: #0033FF; 
                                   padding-top: .3cm; padding-bottom: .3cm;}"),
                        fluidRow(
                          column(6,
                                 textInput(inputId = "gr3",
                                           label="Group 1",
                                           value = "")
                          ),
                          column(6,
                                 textInput(inputId = "gr4",
                                           label = "Group 2",
                                           value = "")
                                 
                          ),
                          column(6,
                                 textInput(inputId = "gr5",
                                           label = "Group 3",
                                           value = "")
                          ),
                          column(6,
                                 textInput(inputId = "gr6",
                                           label = "Group 4",
                                           value = "")
                          ),
                          column(6,
                                 textInput(inputId = "gr7",
                                           label = "Group 5",
                                           value = "")
                          ),
                          
                          column(6,
                                 textInput(inputId = "gr8",
                                           label = "Group 6",
                                           value = "")
                          )
                        ),
                        submitButton("Submit")
                        )
              )),
              fluidRow(
                box(title = "Experiment Design",
                    solidHeader = T, status = "warning", width = 12,
                    box(title = "Matrix Design",
                        solidHeader = T, status = "info", width = 7,
                        dataTableOutput("desm")),
                    box(title = "Contrasts Design",
                        solidHeader = T, status = "info", width = 5,
                        dataTableOutput("contm"))
                    
                    
                )
              ))
            ),
    
    ######End DGE more groups tab#######################
    
    #####Results one group tab#########################
    
    tabItem("res",
            h2("Summary of Results")
    ),
    
    tabItem("res1",
            h2("Results one group comparation"),
            fluidRow(
              box(title = "Toptable",
                  solidHeader = TRUE, status = "warning", width = 12,
                  fluidRow(
                    box(title = "Ebayes Results",
                        solidHeader = T, status = "info", width = 12,
                        withSpinner(dataTableOutput("toptableresults")),
                        downloadButton("table1", label = "Download")),
                    box(title = "Selected Genes",
                        solidHeader = T, status = "info", width = 2,
                        numericInput("pvalue", "Select Adjust P Value:",value = 0.05, min = 0.001, max= 0.7),
                        numericInput("ngenes","Select Number of Genes", value =10,min = 1, max = 20),
                        submitButton("Submit")),
                    box(title = "Selected Genes",
                        solidHeader = T, status = "info", width = 10,
                        withSpinner(dataTableOutput("selectedtable")),
                        downloadButton("table2", label = "Download"))
                  ))
            ),
            fluidRow(
              box(title = "Volcano Plot",
                  solidHeader = T, status = "warning", width = 12,
                  fluidRow(
                    box(title = "Volcano plot",
                        solidHeader = T, status = "info", width =12,
                        withSpinner(plotOutput("volcano")),
                        sliderInput("volcano", "number of genes:",
                                    min = 1, max = 10,
                                    value = 2, step = 1),
                        submitButton("Submit"),
                        br(),
                        downloadButton("plot8", label= "Download")),
                    box(title = "HeatMap",
                        solidHeader = T, status = "info", width =12,
                        withSpinner(plotOutput("heatmap")),
                        downloadButton("plot9", label = "Download"))
                    
                  ))
            )
            
            
    ),
    
    ######End Results one group tab#####################
    
    ######Results more groups tab#######################
    tabItem(tabName = "res2",
            h2("Results Multiple Comparation"),
            fluidRow(
              box("Toptable",
                  solidHeader = T, status = "warning", width = 12,
                  fluidRow(
                    box(title = "Ebayes Results",
                        solidHeader = T, status = "info", width = 12,
                        numericInput("ngroups",label = "Select Group:",
                                     value = 1, min = 1,
                                     max = 3),
                        submitButton("Submit"),
                        br(),
                        withSpinner(dataTableOutput("toptablemresults")),
                        downloadButton("table3", "Download")),
                    
                    box(title = "Selected Genes",
                        solidHeader = T, status = "info", width = 12,
                        numericInput("pvaluem", "Select Adjust P Value:",
                                     value = 0.05, min = 0.01, max = 0.5),
                        numericInput("ngenesm","Selec number of Genes",
                                     value = 10, min = 1, max = 20),
                        submitButton("Submit"),
                        br(),
                        withSpinner(dataTableOutput("selectedtablem")),
                        downloadButton("table4", label = "Download"))
                  ))
            ),
            
            fluidRow(
              box(title = "Volcano Plot",
                  solidHeader = T, status = "warning", width = 12,
                  fluidRow(
                    box(title = "Volcano plot",
                        solidHeader = T, status = "info", width =12,
                        sliderInput("volcanom", "number of genes:",
                                    min = 1, max = 10,
                                    value = 2,step = 1),
                        submitButton("Submit"),
                        (plotOutput("volcanom")),
                        downloadButton("plot10", label= "Download")),
                    box(title = "Heatmap",
                        solidHeader = T, status = "info", width = 12,
                        withSpinner(plotOutput("heatmapm")),
                        downloadButton("plot11", label = "Download"))
                    
                  ))
            ),
            
            
            fluidRow(
              box(title = "Decide Tests",
                  solidHeader = T, status = "warning", width = 12,
                  fluidRow(
                    box(title = "Decide tests summary",
                        solidHeader = T, status = "info", width = 6,
                        numericInput("dectest",label ="Select P.Value", value = 0.01,
                                     min = 0.01, max = 0.5),
                        submitButton("submit"),
                        br(),
                        withSpinner(dataTableOutput("sumdectes")),
                        downloadButton("table5", "Download")),
                    box(title = "Venn Diagram",
                        solidHeader = T, status = "info", width = 6,
                        withSpinner(plotOutput("venn")),
                        downloadButton("plot12", "Download")),
                    box(title = "List of UP and Down  Genes",
                        solidHeader = T, status = "info", width = 12,
                        withSpinner(dataTableOutput("list")),
                        downloadButton("table6","Download"))
                    
                  )
                  
              )
            )
            
    ),
    
    #####End results more groups tab##########
    
    #####TOP GO analysis tab##################
    
    tabItem(tabName = "topgo",
            fluidRow(
              box(title = "Gene Ontology Analysis",
                  solidHeader = T, status = "warning", width = 12,
                  fluidRow(
                    radioButtons("radio", label = "Select Experiment",
                                 choices = list("One Group experiment" =1,
                                                "Multiple groups experiment"=2),
                                 selected=1),
                    submitButton("Submit"),
                    box(title = "Go Table",
                        solidHeader = T, status = "info", width = 12,
                        dataTableOutput("topgo")),
                    box(title = "Go Graph",
                        solidHeader = T, status = "info",width = 12,
                        plotOutput("goplot"))
                  ))
            )),
    tabItem(tabName = "gsea",
            fluidRow(
              box(title = "Enrichment analysis",
                  solidHeader = T, status = "warning", width =12,
                  fluidRow(
                    box(title = "Bar plot",
                        solidHeader = T, status = "info", width = 12,
                        plotOutput("bar")),
                    box(title = "Dot plot",
                        solidHeader = T, status = "info", width = 12,
                        plotOutput("dot")),
                    box(title = "Enrichment map",
                        solidHeader = T, status = "info", width = 12,
                        plotOutput("map")),
                    box(title = "Complex associations",
                        solidHeader = T, status = "info", width = 12,
                        plotOutput("mapc"))
                    
                  ))
            ))
    
    
) 
)



#######End UI################################  

ui <- dashboardPage(header, sidebar, body, skin = "yellow")

#######Server side###############

server <- function(input, output){
  
  ####Load CEL Files function###########
  rawData <- reactive({
    validate(
      need(input$file != "", "Please Load CEL files")
    )
    
    read.celfiles(input$file$datapath)
    
  })
  
  ##
  
  #####Load targets file#################
  
  targets <- reactive({
    validate(
      need(input$file2 !="", "Please load the custom targets file")
    )
    
    targets1 = input$file2
    data1 = read.csv(targets1$datapath,sep = ";")
    return(data1)
  })
  
  
  
  
  #####Show files loades########
  
  output$expression <- renderText({
    data <- targets()$Name
    data1 <- rawData()
    
    samples <- length(colnames(data1))
    paste(data,sep = "")
    
  })
  
  ####Show experiment groups#####
  
  output$groups <- renderText({
    data <- targets()$Group
    
    paste(unique(data), sep = "")
  })
  
  #####Custom targets information#####
  
  output$targets <- renderDataTable({
    targets()
  })
  
  
  #####Histogram plot##############
  output$histogram <-renderPlot({
    withProgress(message = 'Creating Plots',  value = 0.1, {
    Sys.sleep(0.25)
    for (i in 1:15){
      incProgress(1/15)
    }  
    data <- rawData()
    info <- targets()
    hist(data, main="Signal distribution", names = info[,3])
    })
  })
  
  #####Download histogram###############
  
  output$plot1 <- downloadHandler(
    filename = function(){
      paste("plot1", "png", sep = ".")
    },
    content<- function(file){
      png(file)
      
      data <- rawData()
      info <- targets()
      
      plot <- hist(data, main="Signal distribution", names = info[,3], col=as.character(info[,4]))
      
      print(plot)
      
      dev.off()
    },
    contentType = "image/png"
  )
  
  ########BOXPLOT before normalization########
  
  output$boxplotb <- renderPlot({
    data <- rawData()
    info <- targets()
    boxplot(data, cex.axis=0.6, col=as.character(info[,4]), las=2, names=info[,3], main = "Before Normalization")
    
  })
  
  ###########Download boxplot##################
  
  output$plot2 <- downloadHandler(
    filename = function(){
      paste("plot1", "png", sep = ".")
    },
    content<- function(file){
      png(file)
      
      data <- rawData()
      info <- targets()
      
      plot <- boxplot(data, cex.axis=0.6, col=as.character(info[,4]), las=2, names=info[,3],
                      main = "Before Normalization")
      
      print(plot)
      
      dev.off()
    },
    contentType = "image/png"
  )
  
  #######Dendogram before normalization#############
  
  output$dendogramb <- renderPlot({
    data <- rawData()
    info <- targets()
    clust.euclid.average <- hclust(dist(t(exprs(data))),method="average")
    plot(clust.euclid.average, labels=info[,3], main="Hierarchical clustering of samples",  hang=-1)
  })
  
  ##########Download dendogram####################
  
  output$plot3 <- downloadHandler(
    filename = function(){
      paste("plot1", "png", sep = ".")
    },
    content<- function(file){
      png(file)
      
      data <- rawData()
      info <- targets()
      clust.euclid.average <- hclust(dist(t(exprs(data))),method="average")
      plot <- plot(clust.euclid.average, labels=info[,3], main="Hierarchical clustering of samples",  hang=-1)
      print(plot)
      
      dev.off()
    },
    contentType = "image/png"
  )
  
  #######PCA before normalization############
  
  output$PCAB <- renderPlot({
    data <- rawData()
    info <- targets()
    
    plotPCA(exprs(data), groupnames = info[,3])
  })
  
  ####Download PCA########################
  
  output$plot4 <- downloadHandler(
    filename = function(){
      paste("plot1", "png", sep = ".")
    },
    content<- function(file){
      png(file)
      
      data <- rawData()
      info <- targets()
      plot <- plotPCA(exprs(data), groupnames = info[,3])
      print(plot)
      
      dev.off()
    },
    contentType = "image/png"
  )
  
  
  
  ####Normalization##########
  
  eset_rma <- reactive({
    
    withProgress(message = 'Normalizing Raw Data',  value = 0.1, {
    
    incProgress(0.3, "Background correcting")
    incProgress(0.6, "Normalizing")
    incProgress(0.94, "Calculating Expression")
    
    eset_rma <- rawData()
    
    rma(eset_rma)
    
    
    
    })
  })
  
  
  ###Download normalized data
  
  normalized <- reactive({
    norma <- exprs(eset_rma())
    df <- as.data.frame(norma)
    df1 <- cbind(rownames(df), data.frame(df,row.names = NULL))
    return(df1)
  })
  
  output$datanorm <- downloadHandler(
    filename = function(){
      "normalizeddata.csv"},
    content <- function(file){
      write.csv(normalized(),file)
    }
  )
  
  
  ########Boxplot after normalization########
  
  output$boxplotb2 <- renderPlot({
    data <- eset_rma()
    info <- targets()
    boxplot(data, cex.axis=0.6, col=as.character(info[,4]), las=2, names=info[,3], main = "After Normalization")
  })
  
  
  #########Download Boxplot###################
  
  output$plot5 <- downloadHandler(
    filename = function(){
      paste("plot1", "png", sep = ".")
    },
    content<- function(file){
      png(file)
      
      data <- eset_rma()
      info <- targets()
      plot <- boxplot(data, cex.axis=0.6, col=as.character(info[,4]), las=2, names=info[,3], main = "After Normalization")
      print(plot)
      
      dev.off()
    },
    contentType = "image/png"
  )
  
  #######Dendogram after normalization#############
  
  output$dendogramb2 <- renderPlot({
    data <- eset_rma()
    info <- targets()
    clust.euclid.average <- hclust(dist(t(exprs(data))),method="average")
    plot(clust.euclid.average, labels=info[,3], main="Hierarchical clustering of samples",  hang=-1)
  })
  
  ####Download dendogram#####################
  
  output$plot6 <- downloadHandler(
    filename = function(){
      paste("plot1", "png", sep = ".")
    },
    content<- function(file){
      png(file)
      
      data <- eset_rma()
      info <- targets()
      clust.euclid.average <- hclust(dist(t(exprs(data))),method="average")
      plot <- plot(clust.euclid.average, labels=info[,3], main="Hierarchical clustering of samples",  hang=-1)
      print <- plot
      
      dev.off()
    },
    contentType = "image/png"
  )
  
  #######PCA after normalization#############
  
  output$PCAB2 <- renderPlot({
    data <- eset_rma()
    info <- targets()
    
    plotPCA(exprs(data), groupnames = info[,3])
  })
  
  ####Download PCA########################
  
  output$plot7 <- downloadHandler(
    filename = function(){
      paste("plot1", "png", sep = ".")
    },
    content<- function(file){
      png(file)
      data <- eset_rma()
      info <- targets()
      plot <-plotPCA(exprs(data), groupnames = info[,3])
      print <- plot
      dev.off()
    },
    contentType = "image/png"
  )
  
  
  ######Show experiment groups for selection####
  
  output$groupsava <- renderText({
    data <- targets()$Group
    paste("The Avalaible groups level are :", paste(as.character(levels(data)), collapse = ", "),".", sep = "")
  })
  
  ###Show experimetn gropus for multianalysis##########
  output$groupsava2 <- renderText({
    data <- targets()$Group
    paste("The Avalaible groups level are :", paste(as.character(levels(data)), collapse = ", "),".", sep = "")
  })
  
  
  ###Design Matrix######
  
  design <- reactive({
    groups <- as.character(targets()$Group)
    groups1 <- as.factor(groups)
    lev <- factor(groups1, levels = unique(groups1))
    design <- model.matrix(~ 0 + lev)
    colnames(design) <- levels(lev)
    rownames(design) <- targets()$Name
    return(design)
  })
  
  #Make Contrasts
  cont.matrix <- reactive({
    validate(
      need(input$gr2 !="", "Please select the groups")
    )
    
    contrasts <- makeContrasts(
      contrasts =paste(input$gr1, input$gr2, sep="-"),
      levels = design())
  })
  
  cont.matrix.m <- reactive({
    validate(
      need(input$gr3 !="", "Please select the groups")
    )
    x <-c(paste(input$gr3, input$gr4, sep = "-"),paste(input$gr5, input$gr6, sep = "-"),
          paste(input$gr7, input$gr8, sep = "-"))
    
    contrasts <- makeContrasts(
      
      contrasts = x,
      
      levels = design())
    return(contrasts)
  })
  
  
  
  ####DGE objetct##########
  
  fit <- reactive({
    data <- eset_rma()
    des <- design()
    fit1 <- lmFit(exprs(data), des)
    return(fit1)
  })
  
  #####DGE multianalysis########
  
  fit.m <- reactive({
    data <- eset_rma()
    des <- design()
    fit1 <- lmFit(exprs(data), des)
    return(fit1)
  })
  
  
   
  ###Ebayes limma####
  
  fit.main <- reactive({
    fit.main <- contrasts.fit(fit(), cont.matrix())
    return(fit.main)
  })
  
  
  fit.ebayes <- reactive({
    withProgress(message = 'Calculating DGE',  value = 0.1, {
    Sys.sleep(0.25)
    fit.ebayes <- eBayes(fit.main())
    return(fit.ebayes)
    })
  })
  
  
  ###Ebayes limma multianalysis#######
  
  fit.main.m <- reactive({
    fit.main <- contrasts.fit(fit.m(), cont.matrix.m())
    return(fit.main)
  })
  
  
  fit.ebayes.m <- reactive({
    fit.ebayes <- eBayes(fit.main.m())
    return(fit.ebayes)
  })
  
  
  ######Toptable###########
  toptable<- reactive({
    options(digits = 6)
    genes <- topTable(fit.ebayes(), coef = 1, n=1000)
  })
  
  #####Toptable multianalysis###########
  
  toptable.m<- reactive({
    options(digits = 6)
    genes <- topTable(fit.ebayes.m(), coef = input$ngroups, n=1000)
  })
  
  
  #####Anotate toptable##############
  genesymbols <- reactive({
    gensym <- getSYMBOL(rownames(toptable()), "clariomshumantranscriptcluster.db")
  })
  
  ####Anotate toptable multinalysis###
  
  genesymbols.m <- reactive({
    gensym <- getSYMBOL(rownames(toptable.m()), "clariomshumantranscriptcluster.db")
  })
  
  
  #######Toptable anotated########
  
  results <- reactive({
    results1 <- cbind(genesymbols(), toptable())
    results1 <- na.omit(results1)
  })
  
  
  output$toptableresults <- renderDataTable({
    toptable <- results()
    datatable(toptable, rownames = FALSE)
  })
  
  ######Toptable anotated multianalysis###########
  
  results.m <- reactive({
    results1 <- cbind(genesymbols.m(), toptable.m())
    results1 <- na.omit(results1)
  })
  
  
  output$toptablemresults <- renderDataTable({
    toptable <- results.m()
    datatable(toptable[,c(1:3,5:6)], rownames = FALSE)
  })
  
  ###Download toptable############
  
  output$table1 <- downloadHandler(
    filename = function(){
      "toptable.csv"},
    content <- function(file){
      write.csv(results(),file)
    }
  )
  
  ###Download toptable multianalysis
  output$table3 <- downloadHandler(
    filename = function(){
      "toptable.csv"},
    content <- function(file){
      write.csv(results.m(),file)
    }
  )
  
  
  #####Selected genes table######
  
  selected <- reactive({
    
    selected <- results()[which(results()$adj.P.Val <= input$pvalue),]
    selected <- selected[1:input$ngenes,]
  })
  
  
  #####Selected genes table multianalysis###########
  
  selected.m <- reactive({
    
    selected <- results.m()[which(results.m()$adj.P.Val <= input$pvaluem),]
    selected <- selected[1:input$ngenesm,]
  })
  
  #####Output selected genes table#######
  
  output$selectedtable <- renderDataTable({
    data <- selected()[,c(1:6)]
    
    
    datatable(data, rownames = FALSE)
  })
  
  #####Output selected genes table multianalysis#######
  
  output$selectedtablem <- renderDataTable({
    data <- selected.m()[,c(1:3,5:6)]
    
    
    datatable(data, rownames = FALSE)
  })
  
  
  ####Download selected genes table#####
  
  output$table2 <- downloadHandler(
    filename = function(){
      "selectedtable.csv"},
    content <- function(file){
      write.csv(selected(),file)
    }
  )
  
  ####Download selected genes table multianalysis#####
  
  output$table4 <- downloadHandler(
    filename = function(){
      "toptable.csv"},
    content <- function(file){
      write.csv(selected.m(),file)
    }
  )
  
  ######Volcanoplot#################
  
  output$volcano <- renderPlot({
    volcanoplot(fit.ebayes(), coef = 1, highlight = input$volcano, 
                names = names(fit.ebayes()$coefficients[,1]), 
                main = paste("Differentially expressed genes", sep = "\n"))
    abline(v=c(-1,1))
  })
  
  ###Download volcano plot#####################
  
  output$plot8 <- downloadHandler(
    filename = function(){
      paste("plot8", "png", sep = ".")
    },
    content<- function(file){
      png(file)
      
      data <- rawData()
      info <- targets()
      
      plot <- volcanoplot(fit.ebayes(), coef = 1, highlight = input$volcano, 
                          names = names(fit.ebayes()$coefficients[,1]), 
                          main = paste("Differentially expressed genes", sep = "\n"))
      abline(v=c(-1,1))
      print(plot)
      
      dev.off()
    },
    contentType = "image/png"
  )
  
  #########Volcano plot multiple analysis##############
  
  output$volcanom <- renderPlot({
    volcanoplot(fit.ebayes.m(), coef = input$ngroups, highlight = input$volcanom, 
                names = names(fit.ebayes.m()$coefficients[,input$ngroups]), 
                main = paste("Differentially expressed genes", sep = "\n"))
    abline(v=c(-1,1))
  })
  
  ####Download volcano plot###########
  
  output$plot10 <- downloadHandler(
    filename = function(){
      paste("plot10", "png", sep = ".")
    },
    content<- function(file){
      png(file)
      
      data <- rawData()
      info <- targets()
      
      plot <- volcanoplot(fit.ebayes.m(), coef = input$ngroups, highlight = input$volcanom, 
                          names = names(fit.ebayes.m()$coefficients[,1]), 
                          main = paste("Differentially expressed genes", sep = "\n"))
      abline(v=c(-1,1))
      print(plot)
      
      dev.off()
    },
    contentType = "image/png"
  )
  
  ######Heatmap###################
  
  heatmap <- reactive({
    geneID <- rownames(selected())
    exdat <- eset_rma()[geneID]
    colnames(exdat)<- targets()$Name
    gene <- getSYMBOL(rownames(exdat),"clariomshumantranscriptcluster.db")
    rownames(exdat) <-gene
    return(exdat)
  })
  
  exdat <- reactive({
    exdat1 <- exprs(heatmap())
  })
  
  
  output$heatmap <- renderPlot({
    reg2 <- regHeatmap(exdat(), legend = 2, col = heat.colors,
                       breaks = -3:3)
    plot(reg2)
  })
  
  ###Download heatmap##############
  output$plot9 <- downloadHandler(
    filename = function(){
      paste("plot9", "png", sep = ".")
    },
    content<- function(file){
      png(file)
      
      reg2 <- regHeatmap(exdat(), legend = 2, col = heat.colors,
                         breaks = -3:3)
      plot(reg2)
      
      dev.off()
    },
    contentType = "image/png"
  )
  
  #####Heatmap multianalysis######
  
  heatmap.m <- reactive({
    geneID <- rownames(selected.m())
    exdat <- eset_rma()[geneID]
    colnames(exdat)<- targets()$Name
    gene <- getSYMBOL(rownames(exdat),"clariomshumantranscriptcluster.db")
    rownames(exdat) <-gene
    return(exdat)
    
  })
  
  exdat.m <- reactive({
    exdat1 <- exprs(heatmap.m())
  })
  
  output$heatmapm <- renderPlot({
    reg2 <- regHeatmap(exdat.m(), legend = 2, col = heat.colors,
                       breaks = -3:3)
    plot(reg2)
  })
  
  ###Download heatmap###########
  
  output$plot11 <- downloadHandler(
    filename = function(){
      paste("plot11", "png", sep = ".")
    },
    content<- function(file){
      png(file)
      
      reg2 <- regHeatmap.m(exdat(), legend = 2, col = heat.colors,
                           breaks = -3:3)
      plot(reg2)
      
      dev.off()
    },
    contentType = "image/png"
  )
  
  #########Decide tests###################
  
  dectest <- reactive({
    data <- decideTests(fit.ebayes.m(), method = "separate",adjust.method = "fdr", 
                        p.value = input$dectest)
    
  })
  
  
  output$sumdectes <- renderDataTable({
    datatable(summary(dectest()))
  })
  
  
  output$table5 <- downloadHandler(
    filename = function(){
      "toptable.csv"},
    content <- function(file){
      write.csv(dectest(),file)
    }
  )
  
  ####Venn diagram##########
  
  venn <- reactive({
    sum.res.rows <- apply(abs(dectest()),1,sum)
    res.selected <- dectest()[sum.res.rows !=0,]
    return(res.selected)
  })
  
  output$venn <- renderPlot({
    data <- venn()
    vennDiagram(data[,1:3], main="Genes in common", cex =0.9,
                circle.col = c("red","blue","green"), include = c("up","down"))
  })
  
  ###Download venn diagram############
  output$plot12 <- downloadHandler(
    filename = function(){
      paste("plot12", "png", sep = ".")
    },
    content<- function(file){
      png(file)
      
      data <- venn()
      vennDiagram(data[,1:3], main="Genes in common", cex =0.9,
                  circle.col = c("red","blue","green"), include = c("up","down"))
      
      dev.off()
    },
    contentType = "image/png"
  )
  
  #######List of up and down genes###########
  
  res.selected1 <- reactive({
    genesymbols <- getSYMBOL(rownames(venn()),"clariomshumantranscriptcluster.db")
    GeneSymbol <- as.data.frame(genesymbols)
    df <- as.data.frame(venn())
    df1 <- cbind(GeneSymbol, df)
    df1 <- na.omit(df1)
    return(df1)
  })
  
  
  output$list <- renderDataTable({
    datatable(res.selected1(), rownames =FALSE)
  })
  
  ###Download list of up and down genes##########
  
  output$table6 <- downloadHandler(
    filename = function(){
      "toptable.csv"},
    content <- function(file){
      write.csv(res.selected1(),file)
    }
  )
  
  ########Go Analysis#############
  geneID <- reactive ({
    validate(
      need(input$radio !="", "Please choose experiment")
    )
    if(input$radio==1){
    geneid <- as.character(selected()[,1])
    } else {
    geneid <- as.character(selected.m()[,1])
        }
    return(geneid)
  })
 
 input_hyper <- reactive({
   input_hyper = data.frame(geneID(), is_candidate=1)
 })
  
  
  res_hyper <- reactive({
    res_hyper = go_enrich(input_hyper(), n_randset=100)
  })
  
  top_gos_hyper <- reactive({
    res_hyper()[[1]][1:20,"node_id"]
  })
  
  annogenes <- reactive({
    gos <- as.character(top_gos_hyper())
    if(input$radio==1){
    genes <- as.character(selected()[,1])
    }else {
    genes <- as.character(selected.m()[,1])
    }
    anno_genes <- get_anno_genes(go_ids=gos, genes=genes)
    final <- cbind(anno_genes, get_names(anno_genes$go_id)[,2:3])
    return(final)
  })
  
  #####Output#######
  
  output$topgo <- renderDataTable({
    
    annogenes()
    
  })
  
  
  top_gos_hyperplot <- reactive({
    res_hyper()[[1]][1:20,"node_id"]
  })
  
  output$goplot <- renderPlot({
    plot_anno_scores(res_hyper(), top_gos_hyperplot())
  })
  
  ###GSEA#########################
  
  de <- reactive({
    hs <- org.Hs.eg.db
    my.symbols <- as.character(geneID())
    entrez <-select(hs, 
                    keys = my.symbols,
                    columns = c("ENTREZID", "SYMBOL"),
                    keytype = "SYMBOL")
    de <- as.character(entrez$ENTREZID)
    return(de)
  })
  
  ####enrichpathway####
  
  x <- reactive({
    x <- enrichPathway(gene=de(), readable = T, pvalueCutoff = 0.05)
  })
  
  
  ###plots###########
  
  output$bar <- renderPlot({
    barplot(x(), showCategory=8)
  })
  
  output$dot <- renderPlot({
    dotplot(x(), showCategory=15)
  })
  
  output$map <- renderPlot({
    emapplot(x())
  })
  
  output$mapc <- renderPlot({
    cnetplot(x(), categorySize="pvalue", foldChange=targets()[,3])
  })
  
  
   
  
  #######output matrix######
  
  designres <- reactive({
    data <- design()
    data1 <- as.data.frame(data)
    return(data1)
  })
  
  
  output$design <- renderDataTable({
    datatable(designres(), rownames = TRUE)
  })
  
  ####output matrix multianalysis########
  
  designresm <- reactive({
    data <- design()
    data1 <- as.data.frame(data)
    return(data1)
  })
  
  output$desm <- renderDataTable({
    datatable(designresm(), rownames = TRUE)
  })
  
  ####output contrasts#######
  contres <- reactive({
    data <- cont.matrix()
    data1 <- as.data.frame(data)
    return(data1)
  })
  
  output$cont <- renderDataTable({
    datatable(contres(), rownames = TRUE)
  })
  
  ####Output contrasts multianalysis#######
  
  contres.m <- reactive({
    data <- cont.matrix.m()
    data1 <- as.data.frame(data)
    return(data1)
  })
  
  output$contm <- renderDataTable({
    datatable(contres.m(), rownames = TRUE)
  })
  
  
 
##Nav Buttons###

  
  
}



# Run the application 
shinyApp(ui = ui, server = server)

