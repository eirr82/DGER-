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

