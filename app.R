library(shiny)
library(shinythemes)
rm(list=ls())

lambda <- 1
nb.theta <- 100
nb.pbar  <- 100

ui <- fluidPage(
  
  theme = shinytheme("superhero"),
  img(src = "image.jpg"),
  h1("Set the marks of your multiple choice test"),
  h3("Optimized marks derived from a statistical model of scoring."),
  span("Source: ", style="size-font:2.4em;"),
  span(a("A. Direr (2020) \"Efficient Scoring of Multiple-Choice Tests\"", 
    href="https://www.researchgate.net/publication/332222013_Efficient_Scoring_of_Multiple-Choice_Tests"),
    style="size:2.4em;"),
  
  br(),
  br(),
  
  sidebarLayout(
    
    sidebarPanel(
      
      radioButtons(
        "nbOptions", 
        label = "Select the number of options in your items:",
        choices = list(
          "2 options" = 2,
          "3 options" = 3,
          "4 options" = 4),
        selected = 3
      ),
      
      sliderInput(
        "nbItems", 
        label = "Select the number of items in your test:",
        min = 1, max = 100, value = 20
      ),
      
    ),
    
    mainPanel(
      span("Mark for right selection (by convention): ", style="font-size:1.1em;"), 
      br(),
      strong("1", style="color:red;font-size:1.1em;"),
      div("Mark for wrong selection: ", style="font-size:1.1em;"),
      strong(textOutput("resul_theta"), style="color:red; font-size:1.1em;"),
      span("Mark for omission: ", style="font-size:1.1em;"),
      strong(textOutput("resul_gamma"), style="color:red; font-size:1.1em;"),
      span("Targeted proportion of omission (%): ", style="font-size:1.1em;"),
      strong(textOutput("resul_omit"), style="color:red; font-size:1.1em;"),
      span("Total score's mean error: ", style="font-size:1.1em;"),
      span(textOutput("resul_msetot"), style="color:red; font-size:1.1em;"),
      span("Mean error for each item: ", style="font-size:1.1em;"),
      span(textOutput("resul_mse"), style="color:red; font-size:1.1em;")
      
    ),
    
  ),
  
  "F.A.Q.",
  br(), br(),
  em("Why is the mark for wrong answer negative ? "),
  br(),
  "A negative mark corrects for selecting by chance the right answer. It is set so that the score 
  is close to zero if examinees select options at random.",
  br(), br(),
  em("Why is the mark for omission postive ? "),
  br(),
  "The mark gives enough incentives to omit in the case examinees are not confindent 
  enough about which option is right. Omission provides more accurate score than blind guessing.",
  br(), br(),
  em("Why do marks depend on the number of options in each item?"),
  br(),
  "The less options per item, the easier to guess the right option. Wrong selection should be more 
  heavily penalized and omission more encouraged.",
  br(), br(),
  em("Why do marks depend on the number of items in the test?"),
  br(),
  "With more items, examinees' ability is more accurately estimated. Omission should be discouraged
  by lowering the mark for omission.",
  br(), br(),
  em("What is targeted proportion of omission?"),
  br(),
  "It is the proportion of omitted items found efficient and that the marks tries to induce.", 
  br(), br(),
  em("What is total score's mean error?"),
  br(),
  "It measures by how much total score deviates on average from true examinees' ability. 
  For instance, with 20 items and 3 options per item, scores distribute over a 0-20 scale and 
  deviate from true ability by +/-2.67 on average.",  
  br(), br(),
  em("What is mean error for each item?"),
  br(),
  "It measures by how much expected score deviates on average from true examinees' ability for one item.", 
  br(), br(),
  em("What are the main assumptions of the estimation model?"),
  br(),
  "The marks jointly minimize mean square error between scores and abilities averaged over 
  all examinees. Abilities are uniformly distributed. Examinees are risk neutral. See 
  the paper for a full description of the model.",
  br(), br(),
  em("Aside from marks, what other factors are important for scores to be informative?"),
  br(),
  "The quality of questions and answers is of first importance for scores' accuracy.
  Items should be well written, without obvious answers, traps, or ambiguous formulations. 
  Options should be  correctly randomized within each item. Enough time should be granted for 
  all questions to be answered.",
 br(), br(), 
)

server <- function(input, output) {

  reactive_results <- eventReactive(
    c(input$nbItems, input$nbOptions),
    {
      n <- input$nbItems
      m <- as.numeric(input$nbOptions)
      p0 <- 1/m
      tstar <- -1/(m-1)
      theta.inf <- tstar -.1
      theta.sup <- .01
      pas.theta <- (theta.sup - theta.inf)/nb.theta
      v.theta   <- seq(from=theta.inf, to=theta.sup, by=pas.theta)
      pbar.inf <- p0
      pbar.sup <- .85
      pas.pbar <- (pbar.sup - pbar.inf)/nb.pbar
      v.pbar <- seq(from=pbar.inf, to=pbar.sup, by=pas.pbar)
      msemin <- 100
      
      for (pbar in v.pbar)
      {
        (pbar)
        for (theta in v.theta)
        {
          gamma <- pbar + (1-pbar) * lambda * theta
          
          a <- tstar - gamma
          b <- 1 - tstar
          
          sb1 <- (a^2)*p0   + a*b*(p0^2)   + (1/3)*(b^2)*(p0^3)
          sb2 <- (a^2)*pbar + a*b*(pbar^2) + (1/3)*(b^2)*(pbar^3)
          
          A <- (1/n) * ((1-theta)^2)
          B <- (tstar-theta)^2
          
          ms1 <- B*pbar + (1/2)*(A-2*B)*(pbar^2) + (1/3)*(B-A)*(pbar^3)
          ms2 <- B      + (1/2)*(A-2*B)          + (1/3)*(B-A)
          
          mse <- sb2 - sb1 + ms2 - ms1
          
          if ( mse < msemin )
          {
            msemin <- mse
            thetachap <- theta
            pbarchap <- pbar
          }
        }
      }
    
      gammachap <- pbarchap + (1-pbarchap) * lambda * thetachap
      prop.omit <- 100*(pbarchap-p0)/(1-p0)
      msemin <- (msemin/(1-p0))^(1/2)
      msetot <- n*msemin
      
      Espo <- 0.5*(1-tstar)*( pbarchap^2 - p0^2 ) + tstar*(pbarchap - p0)
      if ( pbarchap == p0 ) { bias <- 0
      } else { bias <- 100*(gammachap - Espo/(pbarchap-p0)) }
      
      c(gammachap, thetachap, msemin, prop.omit, msetot)
      
    }
  )
  
  server <- function(input, output) {
    output$slider <- renderUI({
      sliderInput("slider", "Slider", min = 0,
                  max = input$num, value = 0)
    })
  }
  
  output$text_items <- renderText({ 
    paste("Your test includes ", input$nbItems, 
          "items. Each item has ", 
          input$nbOptions, " options.")
  })
  
  output$resul_theta <- renderText({
    results <- reactive_results()
    print(round(results[2], digits=2))
  })
  
  output$resul_gamma <- renderText({
    results <- reactive_results()
    print(round(results[1], digits=2))
  })
  
  output$resul_mse <- renderText({
    results <- reactive_results()
    paste("+/-", round(results[3], digits=2))
  })
  
  output$resul_msetot <- renderText({
    results <- reactive_results()
    paste("+/-", round(results[5], digits=2))
  })
  
  output$resul_omit <- renderText({
    results <- reactive_results()
    print(round(results[4], digits=0))
  })
  
}
shinyApp(ui = ui, server = server)

