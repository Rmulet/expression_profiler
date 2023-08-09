library(shiny)
library(data.table)
library(ggplot2)
library(plotly)
library(cowplot)
library(ggrepel)
library(DT)
options(bitmapType = "cairo")

# setwd("/srv/shiny-server/expression_profiler/")
CCLE = fread("data/CCLE.txt")
LL100 = fread("data/LL100.txt")
AMLRNA = fread("data/AMLRNA.txt")[,-"Id"]
AML.clinical = fread("data/metadata_AML.txt")
theme_set(theme_cowplot())

ui <- fluidPage(
  # Application title
  titlePanel("Expression Profiler"),
  hr(),
  sidebarPanel(
    
    radioButtons("database", "Select the database:",
                 c("LL100" = "LL100","CCLE" = "CCLE","AML RNA-seq" = "AMLRNA")),
    selectizeInput('genes', label='Gene(s) of interest', choices = NULL,
                   multiple=T, options = list(maxOptions=10,placeholder='Type the name of a gene')
    ),
   selectizeInput('lines', label='Cell line(s) or patient(s) of interest', choices = NULL,
                  multiple=T, options = list(placeholder='Type a cell line name or patient ID')
    ),
    actionButton(inputId = "click",
                 label = "Plot"),
    checkboxInput(inputId = "log",
                   label = "Log2 scale"),
  ),

    mainPanel(
      # Output: Tabset w/ plot, summary, and table ----
      tabsetPanel(type = "tabs",
                  tabPanel("Barplot", plotOutput("barplot")),
                  tabPanel("Waterfall plot", 
                    checkboxInput(inputId = "interactive",
                                label = "Interactive"),
                    conditionalPanel(condition = "input.interactive == 1",
                                       plotlyOutput("waterfall_int")),
                   conditionalPanel(condition = "input.interactive == 0",
                                    plotOutput("waterfall")),
                    ),
                  tabPanel("Table", 
          			dataTableOutput("table"),
              			downloadButton("downloadData","Download")		
	                ),
	                tabPanel("Patient data",conditionalPanel(condition = "input.database == 'AMLRNA'",
	                  h5(paste("Here you can view clinical characteristics of in-house patients and"),
	                     strong("select them"),style="color:red"),
	                  verbatimTextOutput('test'),
	                  dataTableOutput("clinical")))
      )
    
    
))

server <- function(input, output,session) {
  data <- reactiveVal()
  observeEvent(input$database, {
    data(get(input$database))
  })
  observeEvent(input$log,{
    if (input$log == TRUE) {
      data(cbind(data()[,"Symbol"],log2(data()[,-"Symbol"]+1)))
    } else {
      data(get(input$database))
    }
    
  })

  observeEvent(input$database,{
    choices.lines <- sort(c(colnames(data())[colnames(data()) != "Symbol"],"Myeloid","Lymphoid"))
    updateSelectizeInput(session,inputId='genes', choices = sort(data()$Symbol), server = TRUE)
    updateSelectizeInput(session,inputId='lines', choices = choices.lines, server = TRUE)
  }) 
  observeEvent(input$click,{ # We can remove this if we want data to be updated dynamically
      d.barplot = data()[Symbol %in% input$genes,c(input$lines,"Symbol"),with=F]
      output$barplot <- renderPlot({
        if(!is.null(input$genes) & !is.null(input$lines)) {
        molten = melt(d.barplot,measure.vars=colnames(d.barplot)[colnames(d.barplot)!="Symbol"])
        ggplot(molten,aes(x=variable,y=value)) + geom_bar(stat="identity",fill="firebrick") +
          labs(x="Cell lines",y="Expression (TPM)") + facet_wrap(~Symbol) # Use facets for multiple genes?s
      }
    })
      d.waterfall = data()[Symbol %in% input$genes,] # Use facets for multiple genes?
      molten = melt(d.waterfall,measure.vars=colnames(d.waterfall)[colnames(d.waterfall)!="Symbol"])
      molten$value <- as.numeric(molten$value)
      molten = molten[order(Symbol,value)]
      molten$term = factor(molten[,paste0(Symbol,"__",variable)],
                           levels=molten[,paste0(Symbol,"__",variable)])
      molten$color <- ifelse(molten$variable  %in% input$lines,"red","black")
      g <- ggplot(molten,aes(x=term,y=value,color=color,label=variable)) + geom_point() + theme(axis.text.x=element_blank()) +
        labs(x="Cell lines",y="Expression (TPM)") + scale_colour_identity() +
        scale_x_discrete(labels = function(x) gsub("__", "", x)) + facet_wrap(~Symbol,scales="free") + 
        theme(legend.position = "none") 
      #molten$color <- sample(c("red","black"),398,repl=T)
      output$waterfall <- renderPlot({
        if(!is.null(input$genes) & !is.null(input$lines)) {
        g + geom_text_repel(data=molten[color == "red"],aes(x=term,y=value,label=variable))
        }
          })
      output$waterfall_int <- renderPlotly({
        if(!is.null(input$genes) & !is.null(input$lines)) {
        ggplotly(g,tooltip=c("label","value"))
        }
          })

    })
      
    dataset <- reactive({data()[Symbol %in% input$genes,c("Symbol",input$lines),with=F]})
    output$table <- renderDataTable(dataset())
    
    output$downloadData <- downloadHandler(
      filename = function() {
        paste(input$database, "_selection.csv", sep = "")
      },
      content = function(file) {
        write.csv(dataset(), file, row.names = FALSE)
      }
    )
    
    output$clinical <- DT::renderDataTable(AML.clinical)
    
    output$test = renderPrint({
      rows = input$clinical_rows_selected
      if (length(rows)) {
        cat('These patients were selected:\n\n')
        cat(unlist(AML.clinical[rows,label]), sep = ', ')
      }
    })
    observe({
    updateSelectizeInput(session,"lines", label = NULL, selected = unlist(
      AML.clinical[input$clinical_rows_selected,label]))
    })
}

# In one table histogram, in another cascade plot
# Create pre-selection of myeloid, lymphoid cell lines?

shinyApp(ui = ui, server = server)

### HELP PAGES ###

# https://stackoverflow.com/questions/52214071/how-to-order-data-by-value-within-ggplot-facets
# https://stackoverflow.com/questions/13904212/r-shiny-plots-do-not-show-in-the-browser

# showModal(modalDialog(
#   title = "my message",
#   "This should only appear when the checkbox is TRUE, given that the_user is FALSE"
# )

# data <- reactiveValues()
# observeEvent(input$database, {
#   data$raw <- get(input$database)
# })
# observeEvent(input$log == TRUE,{
#   data$value <- data$raw
#   for (j in colnames(data$value)[colnames(data$value) != "Symbol"]) set(data$value, j = j, value = log2(data$value[[j]]+0.01))
# })
# observeEvent(input$log == FALSE,{
#   data$value <- data$raw
# })
