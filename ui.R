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
