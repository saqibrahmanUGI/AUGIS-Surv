library(shiny)
library(shinyWidgets)
library(shinycssloaders)
library(shinydashboard)
library(scales)
library(shinyBS)
library(plotly)
library(shinythemes)
library(data.table)
library(formattable)
library(dashboardthemes)
library(shinyjqui)
library(rintrojs)
library(tableone)
library(listviewer)
library(summarytools)
library(plotly)
library(qs)

bsModalNoClose <-function(...) {
  b = bsModal(...)
  b[[2]]$`data-backdrop` = "static"
  b[[2]]$`data-keyboard` = "false"
  return(b)
}



header <- dashboardHeader(
  title = "AUGIS-Surv",
  
  tags$li(a(href = 'http://www.augis.org',
            img(src = 'logo1.png',
                title = "AUGIS", width = "75vw"),
            style = "padding-top:0px; padding-bottom:0x;"),
          class = "dropdown"),
  tags$li(a(href = 'https://www.rcseng.ac.uk',
            img(src = 'RCSLogo.png',
                title = "RCS", width = "75vw"),
            style = "padding-top:0px; padding-bottom:0px;"),
          class = "dropdown"),
  tags$li(a(href = 'https://www.nogca.org.uk',
            img(src = 'nogca.png',
                title = "NOGCA", width = "75vw"),
            style = "padding-top:0px; padding-bottom:0px;"),
          class = "dropdown"),
tags$li(a(href = 'https://www.baso.org.uk',
          img(src = 'BASO.jpg',
              title = "BASO", width = "30vw"),
          style = "padding-top:0px; padding-bottom:0px;"),
        class = "dropdown"))

sidebar<-dashboardSidebar(
  width = 325,
  tags$head(
    tags$style(HTML("
                      .sidebar { height: 90vh; overflow-y: auto; }
                      " )
    )
  ),
  sidebarMenu(
    id='tabs',
    HTML("Version 0.3.1 ©2020<br>See Usage Disclaimer for terms <br>and conditions of use."),
    fluidRow(sliderInput("Age", "Age (years)", min = 20, max =100, value = 60,step = 1)),
    fluidRow(
    div(style='display:inline-block;width:45%;text-align:left;margin-bottom:0px !important;',
    selectInput("Histopathology", "Histology", choices = c("Adenocarcinoma","SCC"),selected='Adenocarcinoma')),
    div(style='display:inline-block;width:45%;text-align:left;margin-bottom:0px !important;',
    selectInput("Gender", "Gender", choices = c("Male","Female"),selected='Male'))),
  fluidRow(
  div(style='display:inline-block;width:45%;text-align:left;margin-bottom:0px !important;',
      selectInput("cT", "cT", choices = c('T0','T1','T2','T3','T4'),selected='T3')),
  div(style='display:inline-block;width:45%;text-align:left;margin-bottom:0px !important;',
      selectInput("cN", "cN", choices = c('N0','N1','N2','N3'),selected='N1'))),
  fluidRow(
    div(style='display:inline-block;width:45%;text-align:left;margin-bottom:0px !important;',
        selectInput("Neoadjuvant.Treatment", "Neoadjuvant Treatment", choices = c("None", "Chemoradiotherapy","Chemotherapy"),selected = 'Chemotherapy')),
        div(style='display:inline-block;width:45%;text-align:left;margin-bottom:0px !important;',
            uiOutput('Completion'))),
  fluidRow( 
    div(style='display:inline-block;width:45%;text-align:left;margin-bottom:0px !important;',
  selectInput("pT.ypT", "pT", choices = c('T0','T1','T2','T3','T4'),selected='T3')),
  div(style='display:inline-block;width:45%;text-align:left;margin-bottom:0px !important;',
  selectInput("Grade", "Grade", choices = c('G1','G2','G3/4','GX'),selected='G1'))),
        sliderInput("Total.LN.Positive", "Positive LN", min = 0, max =
                  30, value = 0,step = 1),
  fluidRow(
    div(style='display:inline-block;width:45%;text-align:left;margin-bottom:0px !important;',
    selectInput("Site.of.Tumour", "Tumour Site", choices=c('GOJ','Oesophagus'),selected='Oesophagus')),
    div(style='display:inline-block;width:45%;text-align:left;margin-bottom:0px !important;',
        selectInput("Any.Complication", "Surgical Complications", choices=c('No','Yes'),selected='No'))),
  
    fluidRow(
      div(style='display:inline-block;width:45%;text-align:left;margin-bottom:0px !important;',
      selectInput("Long", "Logitudinal Margin +ve", choices = c("No","Yes"))),
      div(style='display:inline-block;width:45%;text-align:left;margin-bottom:0px !important;',
        selectInput("Circ", "Circumferential Margin +ve", choices = c("No","Yes")))),
    
    tags$li(a(href = 'https://www.southampton.ac.uk',
              img(src = 'uoslogo.png',
                  title = "UoS", width = '100%'),
              style = "padding-top:0px; padding-bottom:0px;"),
            class = "dropdown")
))

body<-dashboardBody(
  bsModalNoClose("window", "Window",
                 title = "AUGIS-Surv: Prediction of survival after treatment for oesophageal cancer",size='large',
                 p(HTML('By using this app you agree to be bound by the following conditions.<br><br>
The model presented is as accurate as possible at the time of first derivation. It has not yet been subjected to external peer review, but will undergo this shortly. It is provided free to use on an ‘as is’ basis for use by healthcare professionals as a reference tool. Although every effort has been made to ensure its accuracy, no guarantee or warranty is provided.<br><br>
The model is not a replacement for clinical judgement by an experienced professional and is not designed for use of the general public.<br><br>
The output from this model is a prediction and not a guarantee of outcome.<br><br>
No liability for death, injury or other adverse outcome as a result of decisions made using this tool is taken by the authors/publishers/hosting organisation and responsibility remains with the user.<br><br>
The tool is valid in the populations described but may not standardise beyond that and may become less accurate with time or changing medical practice.<br><br>
No data entered into this application is collected or stored.<br><br>
               Application for CE certification is intended, but has not yet been completed.'),),
                 footer = tagList(
                   modalButton("I agree"),
                   actionButton("ok", "I disagree")
                 ),
                 
                
                 tags$head(tags$style("#window .modal-footer{display:none}
                                       .modal-header .close{display:none}
                                      "),
                           tags$script("$(document).ready(function(){
                                        $('#window').modal();
                                        })
                                       $(window).on('load', function(){
                                       $('.modal.fade').appendTo('body');
                                       });"),
                           tags$head(tags$style(HTML('
      .modal.in .modal-dialog{
        width:100%;
        height:100%;
        margin:0px;
      }

      .modal-content{
        width:100%;
        height:100%;
      }
    ')))
                           
                 )),
  
  
  tags$style(type="text/css",
             ".shiny-output-error { visibility: hidden; }",
             ".shiny-output-error:before { visibility: hidden; }"),
  fluidRow(fluidPage(
  
    box(titlePanel(HTML("AUGIS-Surv")),
  tags$p("This model was derived from the National Oesophagogastric Cancer Audit, a compulsory audit of patients undergoing treatment for Oesophageal and Gastric Cancer in England and Wales, using a Random Survival Forest methodology. It is designed to estimate prognosis after completing a course of treatment, and is not valid to guide preoperative treatment decisions."),
  tags$p("Patients treated between 2012 and 2019 who underwent a planned curative Oesophagectomy for Adenocarcinoma or SCC and survived to discharge from hospital were included for model training (n = 6399). It is valid only in this patient group."),
  tags$p("The model was validated internally using bootstrap resampling, with a time dependent AUC at 5 years of 84.8%."),width='100%',collapsible = TRUE),
  tabsetPanel(
        tabPanel('Detailed Prediction and Survival Curve',
    box(withSpinner(fillPage(plotlyOutput("plot")),type=6),
    tableOutput("result"),width='90vh'),
    HTML('<br><br><br><br><br><br>')
        ),
    tabPanel("Conditional Survival",
             box(
               plotlyOutput("condsurv"),
               column(6,dateInput('datesurg','Date of Surgery',value=Sys.Date()-30,min=Sys.Date()-1642.5,max=Sys.Date(),format = "dd/mm/yy",startview = 'decade')),
               column(6,dateInput('dateapp','Date of Assessment',value=Sys.Date(),format = "dd/mm/yy")),
             HTML('<br><br><br><br><br><br>'),
             plotOutput('condcalib'),
             HTML('<br><br><br><br><br><br>'),width='100%')
             ),
    tabPanel('Outcome Pictograms',
             box(tags$p("Pictograms illustrating predicted survival at time points after surgery"),
                 column(4,plotOutput('person1yr')),
                 column(4,plotOutput('person2yr')),
                 column(4,plotOutput('person5yr')),width='100%'),
             HTML('<br><br><br><br><br><br>')),
      tabPanel('Dataset Details',
               box(HTML("POPULATION: Patients undergoing a planned curative oesophagectomy for Adenocarcinoma or Squamous Cell Carcinoma of the Oesophagus or GOJ, survived to discharge and had an adequate lymphadenectomy. In England and Wales between 2012-2019, 8980 patients who underwent an oesophagectomy were identified. Of these, 6838 fulfilled the inclusion criteria and were analysed.<br><br>"),
      HTML("VARIABLE SELECTION: A total of 41 plausible variables were considered for inclusion in the models. <br><br>"),
      box(tableOutput('fullvar'),collapsible = TRUE,collapsed=TRUE,width='100%'),
      HTML("VARIABLE SELECTION: A variable selection step based on Random Forest permutation variable importance was performed (details in manuscript). <br><br>"),
      box(tableOutput('usevar'),collapsible = TRUE,collapsed=TRUE,width='100%'),
      HTML("SUMMARY CHARACTERISTICS. <br><br>"),
      box(uiOutput('table1'),collapsible = TRUE,collapsed=TRUE,width='100%'),
      HTML("MISSING DATA: Missing data was handled by multiple imputation by chained equations (MICE).  Height, Weight, Smoking status and tumour length were not considered for analysis due to excessive missing values (>50%)<br>"),
      box(tags$img(
                    img(src = 'missvar.png',width='100%'),
                    style = "padding-top:0px; padding-bottom:0px;"),width='100%',collapsible=TRUE,collapsed = TRUE),
      HTML('<br><br><br><br><br><br>'),width='100%')),
  tabPanel('Model Details',
           box(HTML("VARIABLE IMPORTANCE: Variable influence is best visualised graphically as marginal effect plots (effect adjusted for all other variables)"),
                    box(tags$img(img(src = 'pdp.png',width='100%'),style = "padding-top:0px; padding-bottom:0px;"),width='100%',collapsible=TRUE,collapsed = TRUE),
           downloadButton("downloadData", "Download Models for External Validation"),
           HTML('<br><br><br><br><br><br>'),width='100%')),
                       
    tabPanel("Usage Disclaimer", box(HTML("By using this app you agree to be bound by the following conditions.<br><br> 
The model presented is as accurate as possible at the time of first publication. It is provided free to use on an ‘as is’ basis for use by healthcare professionals as a reference tool. Although every effort has been made to ensure its accuracy, no guarantee or warranty is provided.<br><br> 
                                       The model is not a replacement for clinical judgement by an experienced professional and is not designed for use of the general public.<br><br>
                                       The output from this model is a prediction and not a guarantee of outcome.<br><br> 
                                       No liability for death, injury or other adverse outcome as a result of decisions made using this tool is taken by the authors/publishers/hosting organisation and responsibility remains with the user. <br><br> 
                                       The tool is valid in the populations described but may not standardise beyond that and may become less accurate with time or changing medical practice.<br><br>
                                       No data entered into this application is collected or stored."),width='100%')
    )
      
      ),style = "height:100vh; overflow-y: auto;"
    )
  ))


dashboardPage(header,sidebar,body,skin='black')