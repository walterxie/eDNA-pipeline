# server.R

library(shiny)
library(DT)

source("TaxonomyAssignment.r")

shinyServer(function(input, output, session) {
  # community matrix of taxa in groups
  # 1st col is taxa, last col is group
  data <- reactive({
    inFile <- input$cmTaxaFile
    
    if (is.null(inFile))
      return(NULL)
    
    taxaAssg <- read.csv(inFile$datapath, header = input$header, sep = input$sep, quote = input$quote, check.names=F)
    
    # set NA (empty cell) to 0
    taxaAssg[is.na(taxaAssg)] <- 0
    
    # rm rows rowSums==0 
#    taxaAssg <- taxaAssg[rowSums(taxaAssg[,-c(1,ncol(taxaAssg))]) > 0,]
    
    taxaAssg
  })
  
  output$cmTaxaTable <- DT::renderDataTable(data())
  
  # taxonomy assignment bar chart
  ta <- reactive({
    percThr = input$percThr / 100
    TaxonomyAssignment(taxaAssg = data(), percThr=percThr, legend_nrow=input$legend_nrow, barValueType=input$bar_vt)
  }) 
  
  # not easy to use pdf in web UI
  output$imageTA <- renderImage({
    if (!is.null(data())) {
      maxLabelLen = ta()$maxlablen
      width = 20 + maxLabelLen * 6 + (ta()$ncol-2) * 120
      height = 20 + input$legend_nrow * 30 + ta()$nrow * 10
      
      outfile <- tempfile(fileext = '.png')
      png(outfile, width = width, height = height)
      print(ta()$plot)
      invisible(dev.off())
      
      width = width * input$zoom / 100
      height = height * input$zoom / 100
      
      # Return a list containing the filename
      list(
        src = outfile,
        contentType = 'image/png',
        width = width,
        height = height,
        alt = "Taxonomy assignment bar chart"
      )
    }
  }, deleteFile = TRUE)
  
  # console output TODO: not reacted properly with input$percThr
#  values <- reactiveValues()
#  output$ta_info <- renderPrint({
#    values[["log"]] <- invisible(capture.output(data <- p()))
#    return( print(unlist(lapply(values[["log"]], paste, collapse=" ")), quote = F) ) 
#  })
  
  # TODO: size chaos
  output$downloadTA <- downloadHandler(
    filename = 'taxa.pdf',
    content = function(con) {
      temp <- tempfile()
      on.exit(unlink(temp))
      
      if (is.null(data()))
        return(NULL)
      maxLabelLen = ta()$maxlablen
      pdfWidth = 0.1 + maxLabelLen / 10 + (ta()$ncol-2) * 1.5
      pdfHeight = 1 + input$legend_nrow * 1 + ta()$nrow * 0.12
      
      pdf(temp, width = pdfWidth, height = pdfHeight)
      print(ta()$plot)
      invisible(dev.off())
      
      bytes <- readBin(temp, "raw", file.info(temp)$size)
      writeBin(bytes, con)
    }
  )
  
  output$summary <- renderPrint({
    summary(data())
  })
  
})