###
# Shiny app purely for testing action buttons
###

library(mgcv)
library(shiny)
library(shinyWidgets)
library(dplyr)
library(ggplot2)
library(scales)


# Code that should run at the start of the app
# load("../data/gam_object_Figure5.rda") # loads an object called fit
load("GRASE_FIT_OBJECT.rda")

# Set graphical theme
theme_set(theme_minimal())

# Note on transformations:
# Model is trained on log(bp) (i.e., natural log), - bp must be genome size
# fraction^6 - fraction must be prob
# log(relative_abundance) - relative_abundance must be target

# define ranges of each variable
# Run this code first; it is static
range.prob <- rev(c(round(10^(-1*seq(0,4,length.out = 1000)),5)))#range(exp(fit$model$probability)) # untransformed range is -4.6 to 0 - exp(probability) seems to make sense
grid.prob <- seq(min(range.prob), max(range.prob), length.out = 100)

range.gs <- c(seq(0.5,20,by=0.1))#range(exp(fit$model$genome_size))
grid.gs <- seq(min(range.gs), max(range.gs), length.out = 100) # Transformed range is 500kb to 20 Mb which makes sense

range.target <- c(round(seq(0.5,1,length.out=1000),3))#range(fit$model$target^{1/6}) # 
grid.target <- seq(min(range.target), max(range.target), length.out = 100)






# Define UI for application that draws a histogram
ui <- fluidPage(
  wellPanel(
    titlePanel("Genome Relative Abundance to Sequencing Effort (GRASE)"),
    numericInput(inputId="read_len",
                 label = "Sequence Read Length (>25)",
                 value=100,
                 min=25,
                 max=Inf),
    numericInput(inputId="prob",
                 label = "Relative Abundance of Genome (fraction from 0.0001 to 1)",
                 value=0.0001,
                 min=0.0001,
                 max=1),
    numericInput(inputId="gs",
                 label = "Genome Size (Mbp; from 0.5 to 20)",
                 value=5.15,
                 min=0.5,
                 max=20),
    numericInput(inputId="target",
                 label = "Desired Fraction of Genome to Sequence (fraction from 0.5 to 1)",
                 value=0.5,
                 min=0.5,
                 max=1),
    # shinyWidgets::sliderTextInput(inputId="prob",
                 # label = "Relative Abundance of Genome (fraction from 0.0001 to 1)",
                 # choices=range.prob,
                           # selected=0.0001,
                           # grid=F),
    # shinyWidgets::sliderTextInput(inputId="gs",
    #                               label = "Genome Size (Mbp)",
    #                               choices=range.gs,
    #                               selected=0.5,
    #                               grid=F),
    # shinyWidgets::sliderTextInput(inputId="target",
    #                               label = "Desired Fraction of Genome to Sequence",
    #                               choices=range.target,
    #                               selected=0.5,
    #                               grid=F),

    actionButton("est.bp", "Estimate Sequences")
  ),

  tags$h3("Estimated Number of Sequences (in millions):"),
  tags$p(),
  tags$h1(textOutput("txt")),
  tags$h3("If you kept the same..."),
  tags$p(),
  tags$p(),
  tags$h5("genome size and desired fraction but different genome relative abundance."),
  plotOutput("p1"),
  tags$p(),
  tags$h5("genome relative abundance and desired fraction but different genome size."),
  # dataTableOutput("p2"),
  plotOutput("p2"),
  tags$p(),
  tags$h5("genome gize and genome relative abundance but different desired fraction of genome to sequence."),
  plotOutput("p3")

)

# Define server logic required to draw a histogram
server <- function(input, output,session) {

  output$readout1 <- reactive({
    paste0("Selected value (numbers): ", input$log_slider, " = ", 10^input$log_slider)
  })

   preds_list  <- eventReactive(input$est.bp,
                         {
                           # Can't think of how to do this any way except by just making 4 prediction vectors
                           # 1, just the predicted value
                           # 2. the predicted value as a function of probability, given the assumed values of genome_size and target
                           # 3. the predicted value as a function of genome_size, given the assumed values of probability and target
                           # 4. the predicted value as a function of (target, given the assumed values of probability and genome_size)
                           preds_list <- list()
                           main_grid_pred <- data.frame(probability = log(input$prob), genome_size = log(input$gs*1e6), target = input$target^6)
                           preds_list$main.pred <- predict(fit, newdata = main_grid_pred)
                           text.pred <- scientific(exp(preds_list$main.pred)/input$read_len/1e6)

                         })



   # Predicted number of base pairs
   output$txt <- renderText({
     preds_list()
   })

   # Plot 1 (what is this?)
   plot1 <- eventReactive(input$est.bp, { #expected bases as a function of probability"
     d1 <- data.frame(genome_size = log(input$gs*1e6), target = input$target^6, probability = log(range.prob))
     d1 <- cbind(d1, data.frame(sequences = exp(predict.gam(fit, newdata = d1))/input$read_len, row.names = NULL))
     ggplot(d1, aes(x = exp(probability), y = sequences/1e6)) +
       geom_line()+
       labs(x="Genome Relative Abundance",y="Sequences")+
       scale_y_log10()+
       scale_x_log10(breaks=c(0.0001,0.001,0.01,0.1,1),labels=c("0.0001","0.001","0.01","0.1","1"))
   })

   # Plot 2
   plot2 <- eventReactive(input$est.bp, {
     d2 <- data.frame(genome_size = log(range.gs*1e6), target = input$target^6, probability = log(input$prob))
     d2 <- cbind(d2, data.frame(sequences = exp(predict.gam(fit, newdata = d2))/input$read_len, row.names = NULL))
     ggplot(d2, aes(x = exp(genome_size)/1e6, y = sequences/1e6)) +
       geom_line()+
       labs(x="Genome Size (Mbp)",y="Sequences")+
       scale_y_log10()+
       scale_x_continuous(breaks=seq(0.5,20,by=0.5))
   })

   # Plot 3
   plot3 <- eventReactive(input$est.bp, {
     d3 <- data.frame(genome_size = log(input$gs*1e6), target = range.target^6, probability = log(input$prob))
     d3 <- cbind(d3, data.frame(sequences = exp(predict.gam(fit, newdata = d3))/input$read_len, row.names = NULL))
     ggplot(d3, aes(x = target^(1/6), y = sequences/1e6)) + 
       geom_line()+
       labs(x="Fractional Genome",y="Sequences")+
       scale_y_log10()+
       scale_x_continuous(breaks=seq(0.5,1,0.05))
   })


   # # Text output
   # output$txt2 <- eventReactive(input$est.bp, {
   #   as.character(scientific(rnorm(1)))
   # })

   # plot1
   output$p1 <- renderPlot({
     plot1()
   })

   # plot2
   output$p2 <- renderPlot({
     plot2()
   })

   # plot3
   output$p3 <- renderPlot({
     plot3()
   })

}

# Run the application
shinyApp(ui = ui, server = server)

