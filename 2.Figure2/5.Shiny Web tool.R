library(shiny)
res <- readRDS("all_res.rds")
ui <- fluidPage(
  titlePanel("Calculate score for guiding ICI therapy"),
  
  sidebarLayout(
    sidebarPanel(
      helpText("If the score is greater than 0.966, it may benefit from ICI therapy."),
      
      tabsetPanel(
        type = "tabs",
        tabPanel("Comparison", 
                 fluidRow(
                   column(6, 
                          selectInput("checkBox1", label = "BMP2|SELE", choices = list("BMP2>SELE"=1, "BMP2<=SELE"=0)),
                          selectInput("checkBox2", label = "CD274|SH3TC1", choices = list("CD274>SH3TC1"=1, "CD274<=SH3TC1"=0)),
                          selectInput("checkBox3", label = "CHST15|LAG3", choices = list("CHST15>LAG3"=1, "CHST15<=LAG3"=0)),
                          selectInput("checkBox4", label = "CKLF|TLR7", choices = list("CKLF>TLR7"=1, "CKLF<=TLR7"=0)),
                          selectInput("checkBox5", label = "ESCO2|RXRA", choices = list("ESCO2>RXRA"=1, "ESCO2<=RXRA"=0))
                   )
                 ),
                 actionButton("submit", "Submit")
        )
      )
    ),
    
    mainPanel(
      h1(textOutput("selected_var"))
    )
  )
)

server <- function(input, output) {
  observeEvent(input$submit, {
    dat <- data.frame(`BMP2|SELE` = as.numeric(input$checkBox1),`CD274|SH3TC1`= as.numeric(input$checkBox2),
                      `CHST15|LAG3`= as.numeric(input$checkBox3),`CKLF|TLR7`=as.numeric(input$checkBox4),
                      `ESCO2|RXRA`=as.numeric(input$checkBox5),check.names = F)
    v1 <- which(dat[,1]==res[,1])
    v2 <- which(dat[,2]==res[,2])
    v3 <- which(dat[,3]==res[,3])
    v4 <- which(dat[,4]==res[,4])
    v5 <- which(dat[,5]==res[,5])
    jj <- table(c(v1,v2,v3,v4,v5))
    jj <- as.numeric(names(jj)[jj==5])
    GMS_score <- res[jj,6]
    
    output$selected_var <- renderText({ 
      paste("ICBP_score= ", as.character(GMS_score))
    })
  })
}

shinyApp(ui = ui, server = server)
