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
    taxaAssg <- taxaAssg[rowSums(taxaAssg[,-c(1,ncol(taxaAssg))]) > 0,]
    
    taxaAssg
  })
  
  output$cmTaxaTable <- DT::renderDataTable(data())
  
  # taxonomy assignment bar chart
  p <- reactive({
    percThr = input$percThr / 100
    TaxonomyAssignment(taxaAssg = data(), percThr, "community matrix", input$legend_nrow)
  }) 
  
  # not easy to use pdf in web UI
  output$imageTA <- renderImage({
    if (!is.null(data())) {
      height = 10 + input$legend_nrow * 20 + nrow(data()) * 4
      maxLabelLen = max(nchar(data()[1]))
      width = maxLabelLen/2 + (ncol(data()) - 2) * 80
      
      outfile <- tempfile(fileext = '.png')
      png(outfile, width = width, height = height)
      print(p())
      invisible(dev.off())
      
      width = width * input$zoom / 100
      height = height * input$zoom / 100
      
      # Return a list containing the filename
      list(
        src = outfile,
        contentType = 'image/png',
        width = width,
        height = height,
        alt = "This is alternate text"
      )
    }
  }, deleteFile = TRUE)
  
  # console output TODO: not reacted properly with input$percThr
  values <- reactiveValues()
  output$ta_info <- renderPrint({
    values[["log"]] <- invisible(capture.output(data <- p()))
    return( print(unlist(lapply(values[["log"]], paste, collapse=" ")), quote = F) ) 
  })
  
  # TODO: size chaos
  output$downloadTA <- downloadHandler(
    filename = 'taxa.pdf',
    content = function(con) {
      temp <- tempfile()
      on.exit(unlink(temp))
      
      if (is.null(data()))
        return(NULL)
      
      pdfHeight = 2 + nrow(data()) * 0.12
      maxLabelLen = max(nchar(data()[1]))
      pdfWidth = 0.1 + maxLabelLen / 10 + (ncol(data()) - 2) * 1.5
      
      pdf(temp, width = pdfWidth, height = pdfHeight)
      print(p())
      invisible(dev.off())
      
      bytes <- readBin(temp, "raw", file.info(temp)$size)
      writeBin(bytes, con)
    }
  )
  
  output$summary <- renderPrint({
    summary(data())
  })
  
})