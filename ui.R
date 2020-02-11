# if (!require("pacman")) install.packages("pacman")
# pacman::p_load(shiny, DT, shinythemes, shinydashboard, shinycssloaders, plotly, ggplot2, rsconnect, ggsci, shinyWidgets)

library(shiny)
library(plotly)
library(DT)
library(shinythemes)
library(shinydashboard)
library(shinycssloaders)
library(ggplot2)
library(rsconnect)
library(ggsci)
library(shinyWidgets)
library(plogr)
source("funciones.R")

# START OF UI 
shinyUI(tagList(
  dashboardPage(
    
  
    
    
    header =   dashboardHeader(title="Meta-DiSc",
                               titleWidth = 350),
    
    sidebar =   dashboardSidebar(
      width = 350,
      sidebarMenu( id = "tabs",
                   
                   #INTRO
                   menuItem(
                     "Introduction", tabName = "intro", icon = icon("info-circle")
                   ),
                   
                   #MENU FILE
                   menuItem(
                     "File upload", tabName ="fidata", icon = icon("folder-open")
                   ),
                   
                   # MENU STATS
                   menuItem(
                     "Statistics", tabName="stats", icon = icon("table"),
                     menuSubItem("Study-level outcome", tabName = "stats1"),
                     menuSubItem("Meta-analysis", tabName = "glm"),
                     menuSubItem("Heterogeneity", tabName = "stats2"),
                     menuSubItem("Analysis of subgroups", tabName = "subgroups"),
                     menuSubItem("Meta-regression", tabName = "metaregression"),
                     menuSubItem("Sensibility", tabName = "stats3")
                   ),
                   
                   #REFERENCES
                   menuItem(
                     "References", tabName = "ref", icon = icon("search-plus")
                   )
      )),
    
    
    body =   dashboardBody(
      
      tags$head(
        tags$link(rel="stylesheet", type = "text/css", href="custom.css")
      ),
      
      
      tabItems(
        
        # MAIN DASHBOARD INTRO
        tabItem( tabName = "intro",
                 h1("Meta-DiSc"), br(),
                 p(
                   "Meta-DiSc",tags$sub("[1]"), "is freeware software to perform Meta-analysis of studies of evaluations of Diagnostic and Screening tests."
                 ),
                 p(
                   "Data must be imported from text files in '.csv' format."
                 ), 
                 p("Meta-DiSc"),
                 
                 tags$ol(type = "a",
                         tags$li("allows exploration of heterogeneity, with a variety of statistics including chi-square, I-squared and Spearman correlation tests,"),
                         tags$li("implements meta-regression techniques to explore the relationships between study characteristics and accuracy estimates, "),
                         tags$li("performs statistical pooling of sensitivities, specificities, likelihood ratios and diagnostic odds ratios using fixed and random effects models, both overall and in subgroups and"),
                         tags$li("produces high quality figures, including forest plots and summary receiver operating characteristic curves that can be exported for use in manuscripts for publication. All computational algorithms have been validated through comparison with different statistical tools and published meta-analyses.")
                 ),
                 p("Meta-DiSc has a Graphical User Interface with roll-down menus, dialog boxes, and online help facilities.")
        ),
        
        # MAIN DASHBOARD FILE
        tabItem( tabName="fidata",
                 fluidRow(
                   column(width = 4,
                          box( fileInput(
                            label = "Select file (.csv)",
                            inputId = "file1",
                            multiple = F, 
                            accept = c("text/csv",
                                       "text/comma-separated-values,text/plain",
                                       ".csv")
                            
                          ),
                          p(
                            "File must contains the columns \"ID, TP, FP, FN, TN\". "
                          ),
                          
                          tags$hr(),
                          actionButton(inputId = "Apply", "Apply"), width =  NULL),
                          
                          box("Download template csv",
                              actionButton(inputId="template", label = "Download Template"), width = NULL)
                   ),
                   column(width = 8,
                          box(withSpinner(DT::dataTableOutput("DataImport"), type = 5),title = "Input Data", solidHeader = TRUE, status = "primary", width = NULL)
                   ))
        ),
        
        # MAIN DASHBOARD STATS
        tabItem(tabName = "stats1",
                fluidRow(
                  box(title="Study-level outcome",withSpinner(DT::dataTableOutput("slevelTable")), status="primary", solidHeader = TRUE,
                      downloadButton("downloadData1", label = "Download"), width = 8, collapsible = TRUE),
                  infoBox("Diagnostic test accuracy summary statistics",
                          value = tags$ul(tags$li("id: Study id"),tags$li("tp: true positives"), tags$li("fp: false positives"),tags$li("fn: false negatives"),
                                          tags$li("tn: true negatives"),tags$li("n1: tp + fn"), tags$li("n0: fp + tn"),tags$li("pos: tp + fp"),tags$li("neg: tn + fn"), 
                                          tags$li("sens: tp/n1"),tags$li("spec: tn/n0"))
                          , width = 4, color = "light-blue")
                ),
                
                tags$br(), tags$br(), tags$br(), tags$br(),
                
                fluidRow(
                  
                  box(title = "Forest Plot of sensitivity", status="primary", solidHeader = TRUE,
                      withSpinner(plotlyOutput("CIsen", height = "1000px")),
                      tags$br(), tags$br()
                      ),
                  box(title = "Forest Plot of specificity", status="primary", solidHeader = TRUE,
                      withSpinner(plotlyOutput("CIspe", height = "1000px")),
                      tags$br(), tags$br()
                      )
                ),
                
                tags$br(), tags$br(), tags$br(), tags$br(), tags$br(), tags$br(), tags$br(), tags$br(), tags$br(), tags$br()
        ),
        
        # META-ANALYSIS
        tabItem(tabName = "glm", 
                fluidRow(
                  box(title = "GLM statistic estimates", width = 4, header = T, status = "primary", solidHeader = TRUE, height = 450,
                      withSpinner(tableOutput("glm_table")),
                      downloadButton("downloadDataGlm", "Download")),
                  box(title = "Parameters for Review Manager", width = 3, header=T, status = "primary", solidHeader = TRUE, height = 450,
                      withSpinner(tableOutput("tablaRevMan")),
                      downloadButton("downloadRevMan", "Download"))
                      ),
                
                fluidRow(
                  box(title = "Options for the SROC plot", width = 4, header = T, status = "primary", solidHeader = TRUE,
                      checkboxGroupInput("sroc_options", label = "Check the options to be shown", selected = c("studies", "ellipse"),
                                         choices = c("Studies" = "studies", "Sens/Spec average" = "averages", "Prediction ellipse" = "ellipse"))
                      ),
                  box(title = "SROC plane",  header = T, status = "primary", solidHeader = T, height = 550, 
                    withSpinner(plotlyOutput("SROC")))
                    ),
                
                br(), br(), br(), br()
        ),
        
        # HETEROGENEITY
        tabItem(tabName = "stats2",
                fluidRow(
                  box(title="Heterogeneity statistics estimates", width = 4, header = T, status = "primary", solidHeader = TRUE,
                      withSpinner(tableOutput("tabla_het")),
                      downloadButton("downloadHet_est", "Download")
                      ),
                  
                  box(title="Heterogeneity statistics estimates","Heterogeneity statistics estimates (95% credible interval) by the Bayesian approach. This model does not take in account the covariates." ,
                      header = T, status = "primary", solidHeader = TRUE,
                      tags$br(), tags$br(),
                      actionButton("heter", "Calculate estimates"),
                      tags$br(), tags$br(),
                      conditionalPanel(condition = "input.heter > 0", withSpinner(tableOutput("Table2"), type = 5), 
                                       downloadButton('downloadDataHet', "Download") ) 
                      ),
                  tags$br(), tags$br(), tags$br(), tags$br(), tags$br(), tags$br(), tags$br(), tags$br(), tags$br(), tags$br()
                  
                )
        ),
        
        # ANALYSIS OF SUBGROUPS
        tabItem(tabName = "subgroups",
                fluidRow(
                  box(title = "Analysis of subgroups", status = "primary", header = T, solidHeader = TRUE,
                      p("Separate meta-analysis for each subgroup."),
                      selectInput(inputId = "subgroups_list", label = "Select  a subgroup to fit the new model.",
                                  choices = NULL, selected = NULL ),
                      actionButton("subgroup_button", label = "Fit model")
                      ),
                  conditionalPanel(condition= "input.subgroup_button > 0", 
                                   box(title= "Analysis completed.", 
                                       textOutput("subgroup"), 
                                       width = 2, status = "success"))
                ),
                fluidRow(conditionalPanel(condition= "input.subgroup_button > 0", 
                                          box(title = "Estimates", status = "primary", header = T, solidHeader = TRUE,
                                              h5("GLM statistic estimates for each subgroup"), 
                                              withSpinner(uiOutput("glm_subgroups")),
                                              br(), br()),
                                          box(title = "SROC plane", status = "primary", header = T, solidHeader = TRUE,
                                              withSpinner(plotlyOutput("sroc_subgroups")))
                                          ),
                         br(), br(), br(), br()), 
                fluidRow(
                  conditionalPanel(condition = "input.comparation > 0",
                                   box(title= "Compare test accuracy", status = "primary", header = T, solidHeader = TRUE)
                                   )
                )
        ),
        
        # META-REGRESSION
        tabItem(tabName = "metaregression",
                fluidRow(
                  box(title = "Meta-regression", status = "primary", header = T, solidHeader = TRUE,
                      withSpinner(tableOutput("combinations")), width = 7, 
                      
                      selectInput(inputId = "covariates_list", label = "Select  a covariate with two levels: ",
                                  choices = NULL, selected = NULL ),
                      actionButton("metareg_button", label = "Fit model")
                  )
                ),
                
                fluidRow(
                  
                  
                  conditionalPanel(condition= "input.metareg_button > 0", 
                                   box(title = "Models", status = "primary", header = T, solidHeader = T, width = 10,
                                       #withSpinner(verbatimTextOutput("error_E")),
                                       withSpinner(verbatimTextOutput("models")),
                                       #withSpinner(verbatimTextOutput("model_A")),
                                       withSpinner(tableOutput("model_final")),
                                       withSpinner(verbatimTextOutput("nsub"))
                                   )
                  )
                ),
                
                br(), br(), br(), br()

        ),

        
        # SENSIBILITY
        tabItem(tabName = "stats3",
                fluidRow(
                  box(title="Sensibility", width = 12, header = T, status = "primary", solidHeader = TRUE,
                      p("The bivariate binomial.......
                        Fit a generalized linear mixed model, which incorporates both fixed-effects parameters and random effects in a linear predictor, via maximum likelihood. The linear predictor is related to the conditional mean of the response through the inverse link function defined in the GLM family.")),
                  box(
                    title = "Selection", header = T, status = "primary", solidHeader = TRUE,
                    checkboxGroupInput(inputId = "studies_list", 
                                       label = "Select the studies that you want to exclude from the model:",
                                       choices = NULL), 
                    width = 4, 
                    actionButton("refresh", "Refresh")
                    
                  ),
                  
                  conditionalPanel(condition = "input.refresh > 0", 
                                   
                                   
                                   box(title="GLM estimates for selected studies", header = T, status = "primary", solidHeader = TRUE,
                                       withSpinner(tableOutput("glm_table_sens")),  width = 4,
                                       downloadButton("downloadDataGlm_sens", "Download")
                                   ),
                                   
                                   box(title="SROC plot of selected studies", header = T, status = "primary", solidHeader = TRUE,
                                       withSpinner(plotlyOutput(outputId = "sroc_sens")), width = 6, height = 550)
                  )
                )
        ),
      
      # REFERENCES
      tabItem(tabName = "ref",
              h2("References"),
              tags$ol(type = "1",
                      tags$li("Zamora, J., Abraira, V., Muriel, A., Khan, K., & Coomarasamy, A. (2006). Meta-DiSc: a software for meta-analysis of test accuracy data. BMC medical research methodology, 6, 31. doi:10.1186/1471-2288-6-31"),
                      tags$li("Douglas Bates, Martin Maechler, Ben Bolker, Steve Walker (2015). Fitting Linear Mixed-Effects Models Using lme4. Journal of Statistical Software, 67(1), 1-48.<doi:10.18637/jss.v067.i01>.")
              ), 
              
              br(), br(),
              h2("Packages used"),
              tags$ul(
                      tags$li("shiny"),
                      tags$li("shiydashboard"),
                      tags$li("plotly"),
                      tags$li("DT"),
                      tags$li("shiythemes"),
                      tags$li("shinycssloaders"),
                      tags$li("ggplot2"),
                      tags$li("rsconnect"),
                      tags$li("ggsci"),
                      tags$li("shinyWidgets"),
                      tags$li(""),
                      )
      )
      
      )
  )
    ), #end dashboardPage
  tags$footer("Meta-DiSc 2.0, 2019.", align = "center", style = "
              position:fixed;
              bottom:0;
              width:100%;
              height:5%;
              color: white;
              padding: 10px;
              background-color: #686A6F;
              z-index: 1000;
              ")
  )
)
  
  
  ## 
  ## 
  ## END UI
  