#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

# Load libraries
library(shiny)
library(ggplot2)
library(gridExtra)
library(shinyBS)
library(magrittr)
library(dplyr)
library(cowplot)
library(stringr)
library(igraph)
library(grid)
library(ggraph)
library(tidyr)
library(networkD3)

# Change data load maximum
options(shiny.maxRequestSize = 500*1024^2)

#Source functions
source("./source.R")

# Define UI for application
ui <- shinyUI(fluidPage(
                        # Set style for vertical alignment
                        tags$head(tags$style(
                          HTML('
                               #vert { 
                               display: flex;
                               align-items: center;
                               margin-top: 50px;
                               }
                               .tooltip .tooltip-inner {
                               max-width: 100%;
                               }
                               '))),
                        shinyjs::useShinyjs(),
                        
                        # Application title
                        titlePanel(textOutput("title")),
                        fluidRow(tags$hr(style="border-color: black;")),
                        
                        # Create Input
                        fluidRow(column(2, wellPanel(
                        fileInput("smeans", 'Choose significant_means.txt',
                                  accept = c('significant_means.txt')),
                        fileInput("pval", 'Choose pvalues.txt',
                                  accept = c('pvalues.txt')),
                        fileInput("meta", 'Choose metadata file',
                                  accept = c('.txt')),
                        #Set load button
                        actionButton("load", "Generate report"),
                        
                        fluidRow(tags$hr(style="border-color: black;")),
    
                        # Select plot to download
                        selectInput("plot",
                                    "Select plot to download",
                                    c("Senders and Recievers", "Network", "Sankey", "Sankey for gene of interest")),
                        #Download
                        downloadButton("download", "Download Plot"))),
                        column(5, align="center", id="vert", plotOutput("sr", width="90%", height="800px")), 
                        column(5, id="vert", align="center", plotOutput("network", width="90%", height="800px"))),
                        
                        # Sankey
                        fluidRow(tags$hr(style="border-color: black;")),
                        fluidRow(column(10, align="center", offset=2, id="vert", h3("Sankey network of interacting celltypes"))),
                        fluidRow(column(10, align="center", offset=2, id="vert", sankeyNetworkOutput("sankey", width="80%", height="700px"))),
                        fluidRow(tags$hr(style="border-color: black;")),
                        fluidRow(column(10, align="center", offset=2, id="vert", h3("Sankey network of interacting celltypes over gene of interest"))),
                        fluidRow(column(2, textInput("gene", "Show feature:", ""), actionButton("search", "Search")), column(10, align="center", id="vert", sankeyNetworkOutput("sk.goi", width="80%", height="700px")))
                        ))

# Define server logic
server <- shinyServer(function(input, output, session) {
  
  #Set text outputs
  output$title <- renderText("Visualisation of CellPhoneDB results")
  
  observeEvent(input$load,{
    pval.file <- isolate({input$pval})
    smeans.file <- isolate({input$smeans})
    metadata.file <- isolate({input$meta})
    
    cp.out <- cellphonedb(pval=pval.file$datapath, smeans=smeans.file$datapath, metadata=metadata.file$datapath)
    
    #Create plots on input of gene
    featurePlotter <- eventReactive(input$load, {
      
      sender <- cp.out$plot.senders
      reciever <- cp.out$plot.recievers
      network <- cp.out$network
      sankey <- cp.out$sankey
      
      list(sender, reciever, network, sankey)
    })
    
    #Set senders and recievers plot
    output$sr <- renderPlot({grid.arrange(featurePlotter()[[1]], featurePlotter()[[2]], nrow=1, 
                                          top = textGrob("Number of interactions from senders and recievers",gp=gpar(fontsize=20)))})
    
    #Set network plot
    output$network <- renderPlot({featurePlotter()[[3]] + 
        ggtitle("Network of interacting celltypes") +
        theme(plot.title = element_text(size = 20))})
    
    #Set sankey plot
    output$sankey <- renderSankeyNetwork({featurePlotter()[[4]]})
    

    #Create plots needed for saving
    observeEvent(input$search,{
      if(input$gene != "") {
        goi.plt.out <- sk.goi(smeans=smeans.file$datapath, metadata=metadata.file$datapath, goi=input$gene)
        } else {
      
      #Set sankey goi plot
        goi.plt.out <- validate(need(exists("gene"), "No feature supplied"))
        }
      #Output sankey
      output$sk.goi <- renderSankeyNetwork({goi.plt.out})
      
      #Create download to include sankey
      savePlot <- function(){
        pnames <- list("Senders and Recievers"=sr, "Network"=network, "Sankey"=sankey, "Sankey for gene of interest"=goi.plt.out)
        pnames[input$plot][[1]]
      }
      
      #Create action for download button
      output$download <- downloadHandler(
        filename = function() {
          fnames <- c("Senders and Recievers"="senders_recievers.pdf", "Network"="network.pdf", "Sankey"="sankey.html", "Sankey for gene of interest"="sankey.html")
          ifelse(input$plot == "Sankey for gene of interest", paste(input$gene, fnames[input$plot][[1]], sep = "_"), fnames[input$plot][[1]])
        },
        content = function(file) {
          
          if (input$plot == "Sankey for gene of interest" | input$plot == "Sankey") {
            saveNetwork(savePlot(), file, selfcontained = TRUE)
            
          } else if (input$plot == "Senders and Recievers"){
            pdf(file)
            grid.draw(savePlot())
            dev.off()
          } else if (input$plot == "Network"){
            pdf(file)
            print(savePlot())
            dev.off()
          }
        })
    })
    
    sr <- grid.arrange(cp.out$plot.senders, cp.out$plot.recievers, nrow=1, 
                       top = textGrob("Number of interactions from senders and recievers",gp=gpar(fontsize=20)))
    
    network <- cp.out$network + 
        ggtitle("Network of interacting celltypes") +
        theme(plot.title = element_text(size = 20))
    
    sankey <- cp.out$sankey
    
    #Create function to save plots as pdfs
    savePlot <- function(){
      pnames <- list("Senders and Recievers"=sr, "Network"=network, "Sankey"=sankey, "Sankey for gene of interest"=NA)
      pnames[input$plot][[1]]
      }
    
    #Create action for download button
    output$download <- downloadHandler(
      filename = function() {
        fnames <- c("Senders and Recievers"="senders_recievers.pdf", "Network"="network.pdf", "Sankey"="sankey.html", "Sankey for gene of interest"="sankey.html")
        ifelse(input$plot == "Sankey for gene of interest", paste(input$gene, fnames[input$plot][[1]], sep = "_"), fnames[input$plot][[1]])
      },
      content = function(file) {
        
        if (input$plot == "Sankey for gene of interest" | input$plot == "Sankey") {
          saveNetwork(savePlot(), file, selfcontained = TRUE)
          
        } else if (input$plot == "Senders and Recievers"){
          pdf(file)
          grid.draw(savePlot())
          dev.off()
        } else if (input$plot == "Network"){
        pdf(file)
        print(savePlot())
        dev.off()
        }
      }
    )
  
  })
})

# Run the application 
shinyApp(ui = ui, server = server)

