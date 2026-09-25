library(shiny)
library(bslib)
library(plotly)
library(shinycssloaders)
library(ggplot2)
library(rhandsontable) # for full manual dependence table
library(tidyr)

# ---- Load method functions ----
if (dir.exists("R")) {
  r_files <- list.files("R", full.names = TRUE, pattern = "\\.R$")
  invisible(lapply(r_files, source))
}

ui <- page_navbar(
  title = "Power & Sample Size Calculator for Ordinal (Composite) Outcomes",
  theme = bs_theme(version = 5, bootswatch = "minty"),
  
  header = tagList(
    tags$style(HTML("
      .navbar {
        background-color: #f2f7f6 !important;
        border-bottom: none;
        border-left:5px solid #78c2ad;
      }

      .navbar-brand {
        color: #34495E !important;
        font-weight:600;
      }
      .navbar-nav {
        margin-top: 10px;
        margin-right: 50px;
      }
      .navbar-nav .nav-link {
        background-color: #eef7f5;
        color: #4F9F8A !important;
        border-radius: 10px;
        margin-left: 8px;
        min-width: 100px;
        padding: 10px 20px;
        font-weight: 600;
        line-height: 1.4;
        text-align: center;
      }

      .navbar-nav .nav-link.active {
        background-color: #66B4A0;
        color: #f7fbfa !important;
        box-shadow: 0 2px 4px rgba(0,0,0,0.08);
      }

      .navbar-nav .nav-link:hover {
        background-color: #dff1ec;
      }
      
      #comp_dist_copy {
    font-family: monospace;
    font-size: 13px;
    background-color: #fbfcfc;
    }
    ")),
    tags$div(
      style = "
        background:#f2f7f6;
        padding:10px 20px;
        margin-top:-16px;
        margin-bottom:15px;
        font-size:17px;
        color:#B85C6B;
        font-weight:600;
        border-left:5px solid #78c2ad;
      ",
      tags$span(style="font-size:15px; margin-right:6px;", "◤"),
      "Supports Individual and Cluster Randomized Trials",
      tags$span(style="font-size:15px; margin-left:6px;", "◢")
    )
  ),
  
  nav_spacer(),   # pushes tabs to the right
  
  nav_panel(
    "Composite",
    sidebarLayout(
      sidebarPanel(
        h4("Ordinal Composite"),
        selectInput("comp_type", "Composition Type", 
                    choices = list("Rank by Severity" = "order", 
                                   "Summation" = "sum", 
                                   "COVID-19 Ordinal Scale Example" = "covid")),
        
        # Only show n_components for Sum and Order
        conditionalPanel(
          condition = "input.comp_type != 'covid'",
          h4("Binary Components"),
          textInput(
            "prev_input",
            "Prevalence (Y1, Y2, ...)",
            value = "0.2, 0.2, 0.2, 0.2, 0.2, 0.2"
          ),
          div(
            style = "margin-top: 6px; margin-bottom: 8px; font-size: 0.95em; color: #7a7a7a;",
            textOutput("n_comp_text", inline = TRUE)
          ),
          helpText("Enter one prevalence per binary component, each between 0 and 1.")
        ),
        
        # Severity Score input only for Order
        conditionalPanel(
          condition = "input.comp_type == 'order'",
          textInput(
            "sev_score",
            "Severity order (least severe to most severe)",
            value = "1, 2, 3, 4, 5, 6"
          ),
          helpText("Example: 1, 2, 3 means component 1 is least severe and component 3 is most severe.")
        ),
        
        conditionalPanel(
          condition = "input.comp_type != 'covid'",
          #numericInput("n_comp", "Number of Components", value = 6, min = 2, max = 10),
          hr(),
          radioButtons(
            "dep_mode",
            h4("Dependence specification"),
            choices = c("Fast & Easy" = "easy", "Full Manual" = "manual"),
            selected = "easy",
            inline = TRUE
          ),
          helpText("In Fast & Easy, same copula family and correlation level are assumed for all trees. Deeper trees can have a fading correlation using decay rate."),
          h5("VC Configuration"),
          selectInput("vine_type", "Vine Copula Type", choices = c("D", "C"))
        ),
        
        conditionalPanel(
          condition = "input.comp_type != 'covid' && input.dep_mode == 'easy'",
          selectInput(
            "cov_fam", "Copula Family",
            choices = copula_choices,
            selected = "Gaussian"
          ),
          h5("Correlation & Decay"),
          selectInput(
            "corr_level", "Correlation Level",
            choices = c(
              "High (τ = 0.6)" = "High",
              "Moderate (τ = 0.4)" = "Moderate",
              "Low (τ = 0.2)" = "Low",
              "Very low (τ = 0.05)" = "Very low",
              "Indep (τ = 0)" = "Indep"
            ),
            selected = "Low"
          ),
          checkboxInput("add_decay", "Add Decay Rate (Deeper Trees)", value = TRUE),
          conditionalPanel(
            condition = "input.add_decay == true",
            sliderInput("decay_val", "Decay Rate", min = 0, max = 1, value = 0.67, step = 0.01)
          )
        ), 
        
        conditionalPanel(
          condition = "input.comp_type != 'covid' && input.dep_mode == 'manual'",
          hr(),
          h4("Manual Bicopula Specification"),
          checkboxInput("manual_use_decay_init", "Initialize tau values using decay", value = TRUE),
          conditionalPanel(
            condition = "input.manual_use_decay_init == true",
            selectInput(
              "manual_init_fam", "Initial Copula Family",
              choices = copula_choices,
              selected = "Gaussian"
            ),
            selectInput(
              "manual_init_corr", "Initial Correlation Level",
              choices = c(
                "High (τ = 0.6)" = "High",
                "Moderate (τ = 0.4)" = "Moderate",
                "Low (τ = 0.2)" = "Low",
                "Very low (τ = 0.05)" = "Very low",
                "Indep (τ = 0)" = "Indep"
              ),
              selected = "Low"
            ),
            sliderInput("manual_decay_val", "Initial Decay Rate", min = 0.01, max = 1, value = 0.67, step = 0.01)
          ),
          actionButton("reset_manual_tbl", "Generate / Reset pair table"),
          br(), br(),
          rHandsontableOutput("paircop_table")
        ),
        
        # COVID-specific prevalence inputs
        conditionalPanel(
          condition = "input.comp_type == 'covid'",
          h5("Prevalence Inputs (%)"),
          numericInput("p_hosp", "Hospitalization", value = 1.4, step = 0.1),
          numericInput("p_vent", "Ventilation/ECMO", value = 0.28, step = 0.1),
          numericInput("p_oxy", "Supp. Oxygen", value = 0.84, step = 0.1),
          numericInput("p_symp", "Symptoms", value = 28.85, step = 1),
          numericInput("p_act", "Activity Limitation", value = 15.13, step = 1),
          
          hr(),
          h5("Subgroup 1: Hospital — Vent/ECMO"),
          selectInput("cov_fam1", "Copula Family", choices = copula_choices, selected = "Gaussian"),
          selectInput("cov_corr1", "Correlation", choices = corr_choices, selected = "Low"),
          
          h5("Subgroup 2: Hospital — Supp. Oxygen"),
          selectInput("cov_fam2", "Copula Family", choices = copula_choices, selected = "Gaussian"),
          selectInput("cov_corr2", "Correlation", choices = corr_choices, selected = "High"),
          
          h5("Subgroup 3: Hospital — Symptoms — Activity"),
          p(tags$i("Pair A: Hospital — Symptoms")),
          selectInput("cov_fam3a", "Copula Family", choices = copula_choices, selected = "Gaussian"),
          selectInput("cov_corr3a", "Correlation", choices = corr_choices, selected = "Low"),
          
          p(tags$i("Pair B: Symptoms — Activity")),
          selectInput("cov_fam3b", "Copula Family", choices = copula_choices, selected = "Gaussian"),
          selectInput("cov_corr3b", "Correlation", choices = corr_choices, selected = "High"),
          
          p(tags$i("Pair C: Hospital — Activity (Cond. on Symptoms)")),
          selectInput("cov_fam3c", "Copula Family", choices = copula_choices, selected = "Indep"),
          selectInput("cov_corr3c", "Correlation", choices = corr_choices, selected = "Indep")
        )
        
        #hr(),
        #actionButton("run_composite", "Generate Distribution", class = "btn-primary")
      ),
      mainPanel(
        navset_card_pill(
          nav_panel("Distribution", plotOutput("distPlot")),
          nav_panel("Vine Structure", plotOutput("vinePlot", height = "900px"))#,
          #nav_panel("Bi-copula Families (Shapes)", plotOutput("bicopPlot", height = "900px"))
        ),
        
        tableOutput("pmfTable"),
        
        tags$hr(),
        
        h5("Copy distribution"),
        
        div(
          style = "
      margin-top: 6px;
      margin-bottom: 6px;
      font-size: 0.95em;
      color: #5D6D7E;
      line-height: 1.45;
    ",
          "Copy and paste this distribution into the corresponding IRT or CRT tab. The copied distribution is ordered from least severe to most severe."
        ),
        
        textAreaInput(
          "comp_dist_copy",
          label = NULL,
          value = "",
          width = "100%",
          height = "70px"
        ),
        
        div(
          style = "
      margin-top: 4px;
      font-size: 0.9em;
      color: #7a7a7a;
      line-height: 1.35;
    ",
          textOutput("comp_dist_copy_note", inline = TRUE)
        )
      )
    )
  ),
  
  nav_panel(
    "IRT",
    
    sidebarLayout(
      sidebarPanel(
        h4("Design inputs"),
        
        h5("Plot set-up"),
        div(
          style = "
    margin-bottom: 12px;
    padding: 8px 10px;
    background: #fbfcfc;
    border-left: 4px solid #78C2AD;
    border-radius: 6px;
    font-size: 0.95em;
    color: #34495E;
  ",
          "Power vs Total sample size"
        ),
        
        sliderInput(
          "irt_n_range",
          "X-axis range: Total sample size (C:T = 1:1)",
          min = 0,
          max = 5000,
          value = c(200, 2000),
          step = 10
        ),
        
        sliderInput(
          "irt_target_power",
          "Target power (%)",
          min = 0,
          max = 100,
          value = 80,
          step = 1
        ),
        
        numericInput(
          "irt_alpha",
          "Type I error",
          value = 0.05,
          min = 0.001,
          max = 0.2,
          step = 0.001
        ),
        
        tags$hr(),
        h5("Ordinal distribution specification"),
        
        radioButtons(
          "irt_dist_source",
          "User specifies:",
          choices = c(
            "Control distribution" = "control",
            "Overall distribution" = "overall"
          ),
          selected = "overall"
        ),
        
        textInput(
          "irt_dist_input",
          "Ordinal distribution (comma-separated)",
          value = "0.25, 0.25, 0.25, 0.25"
        ),
        
        div(
          style = "margin-top: 6px; margin-bottom: 8px; font-size: 0.95em; color: #7a7a7a;",
          "Enter probabilities separated by commas. They must be positive and sum to 1, e.g. 0.30, 0.40, 0.30 or 0.10, 0.20, 0.30, 0.40."
        ),
        
        tags$hr(),
        h5("Treatment effect"),
        
        textInput(
          "irt_OR_text",
          "Odds Ratio",
          value = "1.30"
        )
      ),
      
      mainPanel(
        tags$details(
          open = NA,
          style = "
          margin-bottom: 12px;
          background: #fbfcfc;
          border: 1px solid #dfe7e5;
          border-radius: 8px;
          padding: 8px 12px;
        ",
          tags$summary(
            style = "font-weight:600; cursor:pointer;",
            "Outcome Distributions & Treatment Effect"
          ),
          uiOutput("irt_resolved_inputs")
        ),
        
        tags$br(),
        
        withSpinner(plotlyOutput("irt_power_plot", height = "650px")),
        
        textOutput("irt_plot_summary"),
        
        tags$br(),
        
        downloadButton("download_irt_power_csv", "Download power table (CSV)")
      )
    )
  ),
  
  nav_panel(
    "CRT (3 Cat.)",
    
    sidebarLayout(
      sidebarPanel(
        h4("Design inputs"),
        h5("Plot set-up"),
        selectInput(
          "plot_setup",
          label = NULL,   # or ""
          choices = c(
            "Power vs Total number of clusters" = "power_vs_N",
            "Power vs Cluster size" = "power_vs_J",
            "Total number of clusters vs Cluster size" = "N_vs_J"
          ),
          selected = "power_vs_N"
        ),
        
        conditionalPanel(
          condition = "input.plot_setup == 'power_vs_N'",
          sliderInput(
            "N_range",
            "X-axis range: Total number of clusters (C:T = 1:1)",
            min = 2,
            max = 300,
            value = c(2, 80),
            step = 1
          ),
          sliderInput(
            "J_fixed",
            "Cluster size",
            min = 2,
            max = 500,
            value = 50,
            step = 1
          )
          
        ),
        
        conditionalPanel(
          condition = "input.plot_setup == 'power_vs_J'",
          sliderInput(
            "J_range",
            "X-axis range: Cluster size",
            min = 2,
            max = 500,
            value = c(10, 100),
            step = 1
          ),
          sliderInput(
            "N_fixed",
            "Total number of clusters (C:T = 1:1)",
            min = 2,
            max = 300,
            value = 40,
            step = 1
          )
        ),
        
        conditionalPanel(
          condition = "input.plot_setup == 'N_vs_J'",
          sliderInput(
            "J_range_req",
            "X-axis range: Cluster size",
            min = 0,
            max = 500,
            value = c(50, 200),
            step = 1
          ) #,
          # sliderInput(
          #   "N_search_range",
          #   "Search range: Total number of clusters",
          #   min = 2,
          #   max = 300,
          #   value = c(2, 120),
          #   step = 1
          # )
        ),
        
        sliderInput(
          "target_power",
          "Target power (%)",
          min = 0,
          max = 100,
          value = 80,
          step = 1
        ),
        
        numericInput(
          "alpha", "Type I error",
          value = 0.05, min = 0.001, max = 0.2, step = 0.001
        ),
        
        tags$hr(),
        h5("Ordinal distribution specification"),
        
        radioButtons(
          "dist_source",
          "User specifies:",
          choices = c(
            "Control distribution" = "control",
            "Overall distribution" = "overall"
          ),
          selected = "control"
        ),
        
        textInput(
          "dist_input",
          "Ordinal distribution (comma-separated)",
          value = "0.30, 0.40, 0.30"
        ),
        div(
          style = "margin-top: 6px; margin-bottom: 8px; font-size: 0.95em; color: #7a7a7a;",
          "Enter exactly 3 probabilities separated by commas, e.g. 0.30, 0.40, 0.30"
        ),
        
        tags$hr(),
        h5("Treatment effect"),
        textInput("OR_text", "Odds Ratio", value = "1.30"),
        
        tags$hr(),
        h5("Within-cluster correlation matrix"),
        
        div(
          tags$label(HTML("&rho;<sub>11</sub>")),
          sliderInput(
            "rho11", NULL,
            min = -0.1,
            max = 0.1,
            value = 0.005,
            step = 0.0001
          ),
          div(
            style = "margin-top: -8px; margin-bottom: 10px; font-size: 0.95em;",
            textOutput("rho11_range_text", inline = TRUE)
          )
        ),
        
        div(
          tags$label(HTML("&rho;<sub>12</sub>")),
          sliderInput(
            "rho12", NULL,
            min = -0.1,
            max = 0.1,
            value = -0.001,
            step = 0.0001
          ),
          div(
            style = "margin-top: -8px; margin-bottom: 10px; font-size: 0.95em;",
            textOutput("rho12_range_text", inline = TRUE)
          )
        ),
        
        div(
          tags$label(HTML("&rho;<sub>22</sub>")),
          sliderInput(
            "rho22", NULL,
            min = -0.1,
            max = 0.1,
            value = 0.001,
            step = 0.0001
          ),
          div(
            style = "margin-top: -8px; margin-bottom: 10px; font-size: 0.95em;",
            textOutput("rho22_range_text", inline = TRUE)
          )
        ),
        
        
        div(
          style = "margin-top: 8px; font-size: 0.95em; color: #7a7a7a;",
          uiOutput("rho_validity_text")
        ),
        
        tags$hr(),
        h5("Methods to display"),
        checkboxInput("add_bin", "Include GEE-Binary", value = TRUE),
        checkboxInput("add_WH_DE", "Include Whitehead-Design Effect", value = TRUE),
        conditionalPanel(
          condition = "input.add_WH_DE == true",
          numericInput(
            "WH_DE_ICC",
            "ICC for Whitehead-Design Effect",
            value = 0.005,
            min = -1,
            max = 1,
            step = 0.001
          )
        ),
        checkboxInput("add_indep", "Include GEE-Indep and Whitehead", value = FALSE)
      ),
      
      mainPanel(
        tags$details(
          open = NA,
          style = "
    margin-bottom: 12px;
    background: #fbfcfc;
    border: 1px solid #dfe7e5;
    border-radius: 8px;
    padding: 8px 12px;
  ",
          tags$summary(
            style = "font-weight:600; cursor:pointer;",
            "Outcome Distributions, Treatment Effect & Correlation Structure"
          ),
          uiOutput("resolved_inputs")
        ),
        tags$br(),
        withSpinner(plotlyOutput("main_power_plot", height = "650px")),
        textOutput("main_plot_summary"),
        uiOutput("whitehead_crt_note"),
        tags$br(),
        downloadButton("download_power_csv", "Download power table (CSV)")
      )
    )
  ),
  
  nav_panel(
    "CRT (4+ Cat.)",
    
    sidebarLayout(
      sidebarPanel(
        h4("Design inputs"),
        
        h5("Plot set-up"),
        selectInput(
          "crtK_plot_setup",
          label = NULL,
          choices = c(
            "Power vs Total number of clusters" = "power_vs_N",
            "Power vs Cluster size" = "power_vs_J",
            "Total number of clusters vs Cluster size" = "N_vs_J"
          ),
          selected = "power_vs_N"
        ),
        
        conditionalPanel(
          condition = "input.crtK_plot_setup == 'power_vs_N'",
          sliderInput(
            "crtK_N_range",
            "X-axis range: Total number of clusters (C:T = 1:1)",
            min = 2,
            max = 300,
            value = c(2, 80),
            step = 1
          ),
          sliderInput(
            "crtK_J_fixed",
            "Cluster size",
            min = 2,
            max = 500,
            value = 50,
            step = 1
          )
        ),
        
        conditionalPanel(
          condition = "input.crtK_plot_setup == 'power_vs_J'",
          sliderInput(
            "crtK_J_range",
            "X-axis range: Cluster size",
            min = 2,
            max = 500,
            value = c(10, 100),
            step = 1
          ),
          sliderInput(
            "crtK_N_fixed",
            "Total number of clusters (C:T = 1:1)",
            min = 2,
            max = 300,
            value = 40,
            step = 1
          )
        ),
        
        conditionalPanel(
          condition = "input.crtK_plot_setup == 'N_vs_J'",
          sliderInput(
            "crtK_J_range_req",
            "X-axis range: Cluster size",
            min = 2,
            max = 500,
            value = c(50, 200),
            step = 1
          )
        ),
        
        sliderInput(
          "crtK_target_power",
          "Target power (%)",
          min = 0,
          max = 100,
          value = 80,
          step = 1
        ),
        
        numericInput(
          "crtK_alpha",
          "Type I error",
          value = 0.05,
          min = 0.001,
          max = 0.2,
          step = 0.001
        ),
        
        tags$hr(),
        h5("Original ordinal distribution specification"),
        
        radioButtons(
          "crtK_dist_source",
          "User specifies:",
          choices = c(
            "Control distribution" = "control",
            "Overall distribution" = "overall"
          ),
          selected = "overall"
        ),
        
        textInput(
          "crtK_dist_input",
          "Original ordinal distribution (comma-separated)",
          value = "0.10, 0.20, 0.30, 0.25, 0.15"
        ),
        
        div(
          style = "margin-top: 6px; margin-bottom: 8px; font-size: 0.95em; color: #7a7a7a;",
          "Enter 4–7 probabilities separated by commas. They must be positive and sum to 1."
        ),
        
        tags$hr(),
        h5("Treatment effect"),
        textInput("crtK_OR_text", "Odds Ratio", value = "1.30"),
        
        tags$hr(),
        h5("3-category collapsing specification"),
        
        numericInput(
          "crtK_n_collapse",
          "Number of 3-category collapsing schemes",
          value = 1,
          min = 1,
          max = 5,
          step = 1
        ),
        
        div(
          style = "margin-top: 6px; margin-bottom: 8px; font-size: 0.95em; color: #7a7a7a;",
          HTML(
            "Use semicolons to separate the 3 collapsed categories and commas within each group. 
           Example for 6 categories: <code>1,2; 3,4; 5,6</code>."
          )
        ),
        
        uiOutput("crtK_collapse_inputs"),
        
        tags$hr(),
        h5("Within-cluster correlation matrix for collapsed 3-category outcome"),
        
        div(
          tags$label(HTML("&rho;<sub>11</sub>")),
          sliderInput(
            "crtK_rho11",
            NULL,
            min = -0.1,
            max = 0.1,
            value = 0.005,
            step = 0.0001
          )
        ),
        
        div(
          tags$label(HTML("&rho;<sub>12</sub>")),
          sliderInput(
            "crtK_rho12",
            NULL,
            min = -0.1,
            max = 0.1,
            value = -0.001,
            step = 0.0001
          )
        ),
        
        div(
          tags$label(HTML("&rho;<sub>22</sub>")),
          sliderInput(
            "crtK_rho22",
            NULL,
            min = -0.1,
            max = 0.1,
            value = 0.001,
            step = 0.0001
          )
        ),
        
        div(
          style = "margin-top: 8px; font-size: 0.95em; color: #7a7a7a;",
          "The same collapsed 3-category correlation matrix is applied to all selected 3-category collapsing schemes."
        ),
        
        tags$hr(),
        h5("Methods to display"),
        
        checkboxInput(
          "crtK_add_WH_DE",
          "Include Whitehead-Design Effect",
          value = TRUE
        ),
        
        conditionalPanel(
          condition = "input.crtK_add_WH_DE == true",
          numericInput(
            "crtK_WH_DE_ICC",
            "ICC for Whitehead-Design Effect",
            value = 0.005,
            min = -1,
            max = 1,
            step = 0.001
          )
        ),
        
        checkboxInput(
          "crtK_add_bin",
          "Include GEE-Binary",
          value = FALSE
        ),
        
        conditionalPanel(
          condition = "input.crtK_add_bin == true",
          
          numericInput(
            "crtK_n_binary",
            "Number of binary dichotomizations",
            value = 1,
            min = 1,
            max = 5,
            step = 1
          ),
          
          div(
            style = "margin-top: 6px; margin-bottom: 8px; font-size: 0.95em; color: #7a7a7a;",
            HTML(
              "Use semicolons to separate the two binary groups. 
             Default is first category vs the rest, e.g. <code>1; 2,3,4,5</code>."
            )
          ),
          
          uiOutput("crtK_binary_inputs"),
          
          numericInput(
            "crtK_binary_ICC",
            "ICC for GEE-Binary",
            value = 0.005,
            min = -1,
            max = 1,
            step = 0.001
          )
        )
      ),
      
      mainPanel(
        tags$details(
          open = NA,
          style = "
          margin-bottom: 12px;
          background: #fbfcfc;
          border: 1px solid #dfe7e5;
          border-radius: 8px;
          padding: 8px 12px;
        ",
          tags$summary(
            style = "font-weight:600; cursor:pointer;",
            "Original and Collapsed Outcome Distributions, Treatment Effect & Correlation Structure"
          ),
          uiOutput("crtK_resolved_inputs")
        ),
        
        tags$br(),
        
        withSpinner(plotlyOutput("crtK_power_plot", height = "650px")),
        
        textOutput("crtK_plot_summary"),
        
        tags$br(),
        
        tableOutput("crtK_collapse_table"),
        
        tags$br(),
        
        downloadButton("download_crtK_power_csv", "Download power table (CSV)")
      )
    )
  )
)


server <- function(input, output, session) {
  
  prev_vec <- reactive({
    req(input$comp_type != "covid")
    
    x <- unlist(strsplit(input$prev_input, ","))
    x <- trimws(x)
    x <- x[nzchar(x)]
    x <- as.numeric(x)
    
    shiny::validate(
      need(length(x) >= 2, "Please enter at least 2 prevalence values."),
      need(length(x) <= 10, "Please enter at most 10 prevalence values."),
      need(all(is.finite(x)), "All prevalence values must be numeric."),
      need(all(x > 0 & x < 1), "Each prevalence must be between 0 and 1.")
    )
    
    x
  })
  
  observeEvent(prev_vec(), {
    req(input$comp_type == "order")
    n <- length(prev_vec())
    new_val <- paste(seq_len(n), collapse = ", ")
    
    if (!identical(input$sev_score, new_val)) {
      updateTextInput(session, "sev_score", value = new_val)
    }
  }, ignoreInit = TRUE)
  
  output$n_comp_text <- renderText({
    req(input$comp_type != "covid")
    
    x <- tryCatch({
      vals <- unlist(strsplit(input$prev_input, ","))
      vals <- trimws(vals)
      vals <- vals[nzchar(vals)]
      as.numeric(vals)
    }, error = function(e) numeric(0))
    
    if (length(x) == 0 || any(!is.finite(x))) {
      "Number of components: invalid prevalence input"
    } else {
      paste0("Number of components: ", length(x))
    }
  })
  
  manual_pair_table_init <- reactive({
    req(input$comp_type != "covid")
    n_comp <- length(prev_vec())
    
    pair_df <- generate_vine_pairs(
      n_comp = n_comp,
      vine_type = input$vine_type
    )
    
    if (isTRUE(input$manual_use_decay_init)) {
      tau_base <- tau_values[[input$manual_init_corr]]
      decay <- input$manual_decay_val
      tau_vec <- tau_base * decay^(pair_df$Tree - 1)
      fam_vec <- rep(input$manual_init_fam, nrow(pair_df))
    } else {
      tau_vec <- rep(0.2, nrow(pair_df))
      fam_vec <- rep("Gaussian", nrow(pair_df))
    }
    
    data.frame(
      Tree = pair_df$Tree,
      Pair = pair_df$Pair,
      Copula = fam_vec,
      Tau = round(tau_vec, 3),
      stringsAsFactors = FALSE
    )
  })
  
  manual_tbl <- reactiveVal(NULL)
  
  observeEvent(input$reset_manual_tbl, {
    manual_tbl(manual_pair_table_init())
  }, ignoreInit = FALSE)
  
  observeEvent(list(input$vine_type, prev_vec()), {
    if (input$dep_mode == "manual") {
      manual_tbl(manual_pair_table_init())
    }
  }, ignoreInit = TRUE)
  
  output$paircop_table <- renderRHandsontable({
    req(input$dep_mode == "manual")
    req(manual_tbl())
    
    rhandsontable(manual_tbl(), rowHeaders = FALSE) %>%
      hot_col("Tree", readOnly = TRUE) %>%
      hot_col("Pair", readOnly = TRUE) %>%
      hot_col("Copula", type = "dropdown", 
              source = c("Indep", "Gaussian", "Clayton", "Gumbel", "Frank", "Joe"),
              strict = TRUE,
              allowInvalid = FALSE) %>%
      hot_col("Tau", type = "numeric", format = "0.000")
  })
  
  observeEvent(input$paircop_table, {
    req(input$dep_mode == "manual")
    tbl <- hot_to_r(input$paircop_table)
    req(!is.null(tbl))
    
    old_tbl <- manual_tbl()
    
    if (is.null(old_tbl) || !identical(old_tbl, tbl)) {
      manual_tbl(tbl)
    }
  }, ignoreInit = TRUE)
  
  manual_tbl_clean <- reactive({
    req(input$dep_mode == "manual")
    tbl <- manual_tbl()
    
    shiny::validate(
      need(!is.null(tbl), "Manual pair-copula table is not available."),
      need(all(tbl$Copula %in% copula_choices), "Invalid copula family found in table."),
      need(all(is.finite(tbl$Tau)), "All tau values must be numeric."),
      need(all(tbl$Tau >= 0 & tbl$Tau < 1), "All tau values must be in [0, 1).")
    )
    
    tbl
  })
  
  # --- Reactive to build the RVineMatrix (Sum and Order modes) ---
  rvm_obj <- reactive({
    req(input$comp_type != "covid")
    p_ct <- prev_vec()
    d <- length(p_ct)
    
    # 1. Build the Matrix based on C or D type
    if (input$vine_type == "C") {
      # C-vine: Base variable is the hub (1s in the first column/row)
      mat <- generate_c_vine_matrix(d)
    } else {
      # D-vine: Path-like structure
      mat <- generate_d_vine_matrix(d)
    }
    
    # 2. Get Parameter from Tau
    # Using your name_to_family and tau_values defined earlier
    fam_id <- name_to_family(input$cov_fam) 
    tau_base <- tau_values[[input$corr_level]]
  
    # 3. Build Family and Parameter Matrices
    fam_mat <- matrix(0, d, d)
    par_mat <- matrix(0, d, d)
    
    
    # 3. Fill with Decay Logic
    # In RVineMatrix, rows represent the "Trees"
    if (input$dep_mode == "easy") {
      delta <- if (input$add_decay) input$decay_val else 1
      fam_id <- name_to_family(input$cov_fam)
      tau_base <- tau_values[[input$corr_level]]
      
      for (i in 1:(d - 1)) {
        tau_k <- tau_base * (delta^(i - 1))
        par_k <- if (fam_id == 0 || tau_k <= 0) 0 else BiCopTau2Par(family = fam_id, tau = tau_k)
        
        for (j in (i + 1):d) {
          fam_mat[i, j] <- fam_id
          par_mat[i, j] <- par_k
        }
      }
      
    } else {
      tbl <- manual_tbl_clean()
      
      split_tbl <- split(tbl, tbl$Tree)
      
      for (i in seq_along(split_tbl)) {
        df_tree <- split_tbl[[i]]
        fam_vec <- sapply(df_tree$Copula, name_to_family)
        par_vec <- mapply(function(fam, tau) {
          if (is.na(fam) || fam == 0 || is.na(tau) || tau <= 0) 0
          else BiCopTau2Par(family = fam, tau = tau)
        }, fam = fam_vec, tau = df_tree$Tau)
        
        for (j in seq_along(fam_vec)) {
          col_idx <- i + j
          fam_mat[i, col_idx] <- fam_vec[j]
          par_mat[i, col_idx] <- par_vec[j]
        }
      }
    }
    
    # 4. Create RVineMatrix Object
    RVineMatrix(
      Matrix = mat, 
      family = fam_mat, 
      par = par_mat,
      names = paste0("V", 1:d)
    )
  })
  
  output$vinePlot <- renderPlot({
    if (input$comp_type == "covid") {
      print(
        plot_covid_vine_structure(
          fam1 = input$cov_fam1,
          tau1 = tau_values[[input$cov_corr1]],
          fam2 = input$cov_fam2,
          tau2 = tau_values[[input$cov_corr2]],
          fam3a = input$cov_fam3a,
          tau3a = tau_values[[input$cov_corr3a]],
          fam3b = input$cov_fam3b,
          tau3b = tau_values[[input$cov_corr3b]],
          fam3c = input$cov_fam3c,
          tau3c = tau_values[[input$cov_corr3c]]
        )
      )
    } else {
      rvm <- rvm_obj()
      
      d <- ncol(rvm$Matrix)
      num_trees <- d - 1
      
      cols <- ceiling(sqrt(num_trees))
      rows <- ceiling(num_trees / cols)
      par(mfrow = c(rows, cols), mar = c(2, 2, 4, 2))
      
      for (i in 1:num_trees) {
        plot(rvm, tree = i, type = 0, edge.labels = "family-par")
      }
    }
  },
  height = 900,
  width = 1200,
  res = 96)
  
  rho_grid <- seq(-0.1, 0.1, by = 0.001)
  
  # -----------------------------
  # Parse ordinal distribution
  # -----------------------------
  parse_dist_text <- function(txt, label_prefix = "Distribution") {
    x <- unlist(strsplit(txt, ","))
    x <- trimws(x)
    x <- x[nzchar(x)]
    x <- as.numeric(x)
    
    shiny::validate(
      need(length(x) == 3, paste0(label_prefix, " must contain exactly 3 probabilities.")),
      need(all(is.finite(x)), paste0(label_prefix, " entries must all be numeric.")),
      need(all(x > 0), paste0(label_prefix, " entries must all be > 0.")),
      need(abs(sum(x) - 1) < 1e-8, paste0(label_prefix, " must sum to 1."))
    )
    
    x
  }
  
  dist_vec <- reactive({
    parse_dist_text(
      input$dist_input,
      if (input$dist_source == "control") "Control distribution" else "Overall distribution"
    )
  })
  
  theta_val <- reactive({
    or_val <- as.numeric(trimws(input$OR_text))
    shiny::validate(
      need(length(or_val) == 1 && is.finite(or_val) && or_val > 0,
           "Treatment effect (OR) must be a positive number.")
    )
    log(or_val)
  })
  
  
  
  # -----------------------------
  # Resolved inputs
  # -----------------------------
  resolved_inputs <- reactive({
    shiny::validate(
      need(exists("resolve_3cat_inputs"), "resolve_3cat_inputs() not found. Did you source your R/ files?"),
      need(is.function(resolve_3cat_inputs), "resolve_3cat_inputs exists but is not a function.")
    )
    
    if (input$dist_source == "control") {
      resolve_3cat_inputs(
        p_C_vec = dist_vec(),
        treatment_effect = theta_val()
      )
    } else {
      resolve_3cat_inputs(
        pi_vec = dist_vec(),
        treatment_effect = theta_val()
      )
    }
  })
  
  resolved_display_inputs <- reactive({
    x <- resolved_inputs()
    
    calculate_power_3cat(
      p_C_vec          = x$p_C_vec,
      treatment_effect = x$treatment_effect,
      number_cluster   = 20,   # dummy valid value
      cluster_size     = 50,   # dummy valid value
      rho11            = input$rho11,
      rho12            = input$rho12,
      rho22            = input$rho22,
      add_bin          = input$add_bin,
      bin_rho1         = NULL,
      bin_rho2         = NULL,
      add_WH_DE = input$add_WH_DE,
      WH_DE_ICC = input$WH_DE_ICC,
      add_indep        = input$add_indep,
      alpha            = input$alpha
    )
  })
  
  # -----------------------------
  # WH-DE ICC behavior
  # -----------------------------
  wh_de_icc_touched <- reactiveVal(FALSE)
  
  observeEvent(input$WH_DE_ICC, {
    wh_de_icc_touched(TRUE)
  }, ignoreInit = TRUE)
  
  observeEvent(input$rho11, {
    if (!wh_de_icc_touched()) {
      updateNumericInput(session, "WH_DE_ICC", value = input$rho11)
    }
  }, ignoreInit = FALSE)
  
  # -----------------------------
  # Decide distribution type
  # -----------------------------
  dist_type_info <- reactive({
    x <- resolved_inputs()
    classify_ordinal_distribution(x$pi_vec)
  })
  
  output$resolved_inputs <- renderUI({
    x <- resolved_display_inputs()
    dist_info <- dist_type_info()
    
    fmt_vec <- function(v, digits = 4) {
      paste0("(", paste(sprintf(paste0("%.", digits, "f"), v), collapse = ", "), ")")
    }
    
    fmt_num <- function(v, digits = 4) {
      sprintf(paste0("%.", digits, "f"), v)
    }
    
    tags$div(
      style = "
      font-size: 15px;
      color: #4a4a4a;
      line-height: 1.45;
      margin-top: 6px;
    ",
      
      # line 1: distributions
      tags$div(
        style = "margin-bottom: 10px;",
        tags$span(style = "font-weight:600; color:#34495E;", "Control:"),
        tags$span(paste0(" ", fmt_vec(x$p_C_vec), "   |   ")),
        
        tags$span(style = "font-weight:600; color:#34495E;", "Treatment:"),
        tags$span(paste0(" ", fmt_vec(x$p_T_vec), "   |   ")),
        
        tags$span(style = "font-weight:600; color:#34495E;", "Overall:"),
        tags$span(paste0(" ", fmt_vec(x$pi_vec)))
      ),
      tags$div(
        style = "margin-bottom: 4px;",
        tags$span(style = "font-weight:600; color:#34495E;", "Closest overall ordinal distribution type:"),
        tags$span(
          paste0(
            " ",dist_info$assigned_type,
            "   |   Similarity: ", sprintf("%.1f", dist_info$similarity), "%"
          )
        )
      ),
      # ---- Line 2: OR ----
      tags$div(
        style = "margin-bottom: 10px;",
        tags$span(style = "font-weight:600; color:#34495E;", "Treatment effect (OR):"),
        tags$span(
          HTML(paste0(
            " ", fmt_num(exp(x$treatment_effect)),
            " &nbsp;&nbsp; (log = ", fmt_num(x$treatment_effect), ")"
          ))
        )
      ),
      # rho matrix
      tags$div(
        style = "margin-bottom: 8px;",
        tags$div(
          style = "font-weight:600; color:#34495E; margin-bottom: 4px;",
          "Within-cluster correlation matrix:"
        ),
        tags$div(
          style = "font-family: 'Times New Roman', serif; font-size: 18px; margin-left: 8px;",
          HTML(paste0(
            "<table style='border-collapse:collapse; display:inline-table; vertical-align:middle;'>",
            "<tr>",
            "<td rowspan='2' style='font-size:30px; padding-right:6px;'>(</td>",
            "<td style='padding:2px 12px;'>&rho;<sub>11</sub> = ", fmt_num(input$rho11), "</td>",
            "<td style='padding:2px 12px;'>&rho;<sub>12</sub> = ", fmt_num(input$rho12), "</td>",
            "<td rowspan='2' style='font-size:30px; padding-left:6px;'>)</td>",
            "</tr>",
            "<tr>",
            "<td style='padding:2px 12px;'>&rho;<sub>12</sub> = ", fmt_num(input$rho12), "</td>",
            "<td style='padding:2px 12px;'>&rho;<sub>22</sub> = ", fmt_num(input$rho22), "</td>",
            "</tr>",
            "</table>"
          ))
        )
      ),
      
      # explanations
      tags$div(
        style = "font-size: 14px; margin-top: 6px; margin-bottom: 10px;",
        HTML(paste0(
          "<div>&rho;<sub>11</sub>: correlation between the 1st binary indicator (Category 1) across subunits within the same cluster.</div>",
          
          "<div style='margin-top:4px;'>&rho;<sub>12</sub>: correlation between the 1st and 2nd binary indicators across subunits within the same cluster.</div>",
          
          "<div style='margin-top:4px;'>&rho;<sub>22</sub>: correlation between the 2nd binary indicator (Category 2) across subunits within the same cluster.</div>"
        ))
      ),
      
      if (isTRUE(input$add_bin)) {
        tagList(
          tags$div(
            style = "
              margin-top:14px;
              margin-bottom:8px;
              padding-top:10px;
              border-top:1px solid #D9DEE3;
              font-weight:700;
              color:#34495E;
            ",
            "Dichotomization"
          ),
          
          tags$div(
            style = "margin-bottom:10px;",
            tags$span(style = "font-weight:600; color:#34495E;", "Binary1 log-OR:"),
            tags$span(paste0(" ", fmt_num(x$theta_bin1), " | ")),
            tags$span(style = "font-weight:600; color:#34495E;", "Binary2 log-OR:"),
            tags$span(paste0(" ", fmt_num(x$theta_bin2)))
          ),
          
          tags$div(
            style = "margin-bottom:10px;",
            tags$span(style = "font-weight:600; color:#34495E;", "Binary1 ICC (0 vs 1,2):"),
            tags$span(paste0(" ", fmt_num(x$bin_rho1), " | ")),
            tags$span(style = "font-weight:600; color:#34495E;", "Binary2 ICC (0,1 vs 2):"),
            tags$span(paste0(" ", fmt_num(x$bin_rho2)))
          )#,
      #     tags$div(
      #       style = "margin-top:3px; margin-left:2px; font-size:11px; color:#5D6D7E; line-height:1.45; text-align:left; display:block;",
      #       tags$span(style = "margin-top:2px;","Note: "),
      #       withMathJax(
      #         helpText(
      #           "$$\\text{Binary1 ICC}=\\rho_{11},\\quad\\quad
      # \\text{Binary2 ICC}=
      # \\frac{
      # \\rho_{11}p_1(1-p_1)
      # +2\\rho_{12}\\sqrt{p_1(1-p_1)p_2(1-p_2)}
      # +\\rho_{22}p_2(1-p_2)
      # }{
      # (p_1+p_2)p_3
      # }.$$"
      #         )
      #       )
      #     )
        )
      },
      
      if (isTRUE(input$add_WH_DE)) {
        tags$div(
          style = "margin-top:10px;padding-top:10px;border-top:1px solid #D9DEE3;",
          tags$span(style = "font-weight:600; color:#34495E;", "Whitehead-Design Effect's ICC:"),
          tags$span(paste0(" ", fmt_num(x$WH_DE_ICC)))
        )
      }
      
    )
  })
  
  
  
  
  # -----------------------------
  # Helper for calculate_power_3cat()
  # -----------------------------
  calc_power_once <- function(N, J, x) {
    out <- calculate_power_3cat(
      p_C_vec         = x$p_C_vec,
      treatment_effect = x$treatment_effect,
      number_cluster  = N,
      cluster_size    = J,
      rho11           = input$rho11,
      rho12           = input$rho12,
      rho22           = input$rho22,
      add_bin         = input$add_bin,
      bin_rho1        = NULL,
      bin_rho2        = NULL,
      add_WH_DE = input$add_WH_DE,
      WH_DE_ICC = input$WH_DE_ICC,
      add_indep       = input$add_indep,
      alpha           = input$alpha
    )
    
    to_pct <- function(z) {
      if (is.null(z) || length(z) == 0 || is.na(z)) return(NA_real_)
      round(100 * as.numeric(z), 2)
    }
    
    data.frame(
      N           = as.integer(N),
      J           = as.integer(J),
      GEE_Ordinal = to_pct(out$power_GEE_exch),
      GEE_Binary1 = to_pct(out$power_GEE_bin1),
      GEE_Binary2 = to_pct(out$power_GEE_bin2),
      WH_DE       = to_pct(out$power_WH_DE),
      GEE_Indep   = to_pct(out$power_GEE_indep),
      Whitehead   = to_pct(out$power_WH),
      check.names = FALSE
    )
  }
  
  # -----------------------------
  # X grids
  # -----------------------------
  N_seq <- reactive({
    seq(input$N_range[1], input$N_range[2], by = 1)
  })
  
  J_seq <- reactive({
    seq(input$J_range[1], input$J_range[2], by = 1)
  })
  
  J_seq_req <- reactive({
    seq(input$J_range_req[1], input$J_range_req[2], by = 1)
  })
  
  # N_search_seq <- reactive({
  #   seq(input$N_search_range[1], input$N_search_range[2], by = 1)
  # })
  N_search_seq <- function(J) {
    upper <- max(100, ceiling(2000 / J))  # adaptive upper bound
    seq(2, upper, by = 1)
  }
  
  # -----------------------------
  # Valid rho range
  # -----------------------------
  valid_rho_range_local <- reactive({
    x <- resolved_inputs()
    
    cluster_size_for_check <- switch(
      input$plot_setup,
      "power_vs_N" = input$J_fixed,
      "power_vs_J" = input$J_range[1],
      "N_vs_J"     = input$J_range_req[1] ### need to explain this somewhere!!!!
    )
    
    function(target = c("rho11", "rho12", "rho22")) {
      target <- match.arg(target)
      
      vals <- sapply(rho_grid, function(candidate) {
        r11 <- input$rho11
        r12 <- input$rho12
        r22 <- input$rho22
        
        if (target == "rho11") r11 <- candidate
        if (target == "rho12") r12 <- candidate
        if (target == "rho22") r22 <- candidate
        
        rho_mat <- matrix(c(r11, r12, r12, r22), nrow = 2, byrow = TRUE)
        
        M <- tryCatch(
          gee_info_mat_3cat(
            pi_T1 = x$p_T_vec[1],
            pi_T2 = x$p_T_vec[2],
            pi_C1 = x$p_C_vec[1],
            pi_C2 = x$p_C_vec[2],
            rho_mat = rho_mat,
            cluster_size = cluster_size_for_check
          ),
          error = function(e) NULL
        )
        
        !is.null(M) && is_valid_gee_info_mat(M)
      })
      
      valid_vals <- rho_grid[vals]
      
      if (length(valid_vals) == 0) {
        return(list(min = NA_real_, max = NA_real_, valid_values = numeric(0)))
      }
      
      list(
        min = min(valid_vals),
        max = max(valid_vals),
        valid_values = valid_vals
      )
    }
  })
  
  rho11_range <- reactive(valid_rho_range_local()("rho11"))
  rho12_range <- reactive(valid_rho_range_local()("rho12"))
  rho22_range <- reactive(valid_rho_range_local()("rho22"))
  
  output$rho11_range_text <- renderText({
    rr <- rho11_range()
    if (is.na(rr$min)) {
      "No valid range found on the search grid."
    } else {
      sprintf("Given current \u03c1\u2081\u2082 and \u03c1\u2082\u2082, valid range: [%.3f, %.3f]", rr$min, rr$max)
    }
  })
  
  output$rho12_range_text <- renderText({
    rr <- rho12_range()
    if (is.na(rr$min)) {
      "No valid range found on the search grid."
    } else {
      sprintf("Given current \u03c1\u2081\u2081 and \u03c1\u2082\u2082, valid range: [%.3f, %.3f]", rr$min, rr$max)
    }
  })
  
  output$rho22_range_text <- renderText({
    rr <- rho22_range()
    if (is.na(rr$min)) {
      "No valid range found on the search grid."
    } else {
      sprintf("Given current \u03c1\u2081\u2081 and \u03c1\u2081\u2082, valid range: [%.3f, %.3f]", rr$min, rr$max)
    }
  })
  
  # -----------------------------
  # Current rho combination validity
  # -----------------------------
  current_rho_valid <- reactive({
    x <- resolved_inputs()
    
    cluster_size_for_check <- switch(
      input$plot_setup,
      "power_vs_N" = input$J_fixed,
      "power_vs_J" = max(2, input$J_range[1]),
      "N_vs_J"     = max(2, input$J_range_req[1])
    )
    
    rho_mat <- matrix(c(input$rho11, input$rho12,
                        input$rho12, input$rho22), nrow = 2, byrow = TRUE)
    
    M <- tryCatch(
      gee_info_mat_3cat(
        pi_T1 = x$p_T_vec[1],
        pi_T2 = x$p_T_vec[2],
        pi_C1 = x$p_C_vec[1],
        pi_C2 = x$p_C_vec[2],
        rho_mat = rho_mat,
        cluster_size = cluster_size_for_check
      ),
      error = function(e) NULL
    )
    
    !is.null(M) && is_valid_gee_info_mat(M)
  })
  
  output$rho_validity_text <- renderUI({
    if (current_rho_valid()) {
      HTML("Current &rho; combination is valid and yields positive V<sub>&theta;</sub>.")
    } else {
      HTML("Current &rho; combination is not valid: it does not yield positive V<sub>&theta;</sub>.")
    }
  })
  
 
  
  # -----------------------------
  # Plot data
  # -----------------------------
  title_text <- reactive({
    target_pct <- round(input$target_power, 2)
    OR_val <- round(exp(theta_val()), 2)
    
    if (input$plot_setup == "power_vs_N") {
      paste0(
        "Power vs Total Number of Clusters",
        " (Cluster size = ", input$J_fixed,
        " | Target power = ", target_pct, "%)"
      )
      
    } else if (input$plot_setup == "power_vs_J") {
      paste0(
        "Power vs Cluster Size",
        " (Total clusters = ", input$N_fixed,
        " | Target power = ", target_pct, "%)"
      )
      
    } else {
      paste0(
        "Required Number of Clusters vs Cluster Size",
        " (Target power = ", target_pct, "%)"
      )
    }
  })
  
  method_list <- reactive({
    methods <- c("GEE_Ordinal")
    
    if (input$add_bin) {
      methods <- c(methods, "GEE_Binary1", "GEE_Binary2")
    }
    
    if (input$add_WH_DE) {
      methods <- c(methods, "WH_DE")
    }
    
    if (input$add_indep) {
      methods <- c(methods, "GEE_Indep", "Whitehead")
    }
    
    methods
  })
  
  find_required_N_by_search <- function(J, x, methods, target_pct, 
                                        coarse_step = 5,
                                        N_min = 2,
                                        N_max = max(100, ceiling(2000 / J))) {
    
    found <- setNames(as.list(rep(NA_real_, length(methods))), methods)
    
    coarse_grid <- seq(N_min, N_max, by = coarse_step)
    if (tail(coarse_grid, 1) != N_max) {
      coarse_grid <- c(coarse_grid, N_max)
    }
    
    last_row <- NULL
    last_N <- NA_integer_
    
    for (N in coarse_grid) {
      row <- calc_power_once(N = N, J = J, x = x)
      
      for (m in methods) {
        if (is.na(found[[m]])) {
          val <- row[[m]]
          
          if (!is.null(val) && !is.na(val) && val >= target_pct) {
            # refine within previous coarse interval
            lower_N <- if (is.na(last_N)) N_min else max(N_min, last_N + 1)
            refine_grid <- seq(lower_N, N, by = 1)
            
            for (N_ref in refine_grid) {
              row_ref <- calc_power_once(N = N_ref, J = J, x = x)
              val_ref <- row_ref[[m]]
              
              if (!is.null(val_ref) && !is.na(val_ref) && val_ref >= target_pct) {
                found[[m]] <- N_ref
                break
              }
            }
          }
        }
      }
      
      if (all(!is.na(unlist(found)))) break
      
      last_row <- row
      last_N <- N
    }
    
    cbind(
      data.frame(J = J),
      as.data.frame(found, check.names = FALSE)
    )
  }
  
  power_df <- reactive({
    shiny::validate(
      need(exists("calculate_power_3cat"), "calculate_power_3cat() not found. Did you source your R/ files?"),
      need(is.function(calculate_power_3cat), "calculate_power_3cat exists but is not a function."),
      need(current_rho_valid(), "Current rho combination does not yield positive Vθ.")
    )
    
    x <- resolved_inputs()
    
    if (input$plot_setup == "power_vs_N") {
      res_list <- lapply(N_seq(), function(N) {
        calc_power_once(N = N, J = input$J_fixed, x = x)
      })
      df <- do.call(rbind, res_list)
      df$x_value <- df$N
      df$x_label <- "Total number of clusters"
      
    } else if (input$plot_setup == "power_vs_J") {
      res_list <- lapply(J_seq(), function(J) {
        calc_power_once(N = input$N_fixed, J = J, x = x)
      })
      df <- do.call(rbind, res_list)
      df$x_value <- df$J
      df$x_label <- "Cluster size"
      
    } else {
      target_pct <- round(input$target_power, 2)
      
      prev_solution <- NULL

      res_list <- lapply(J_seq_req(), function(J) {
        
        start_N <- if (is.null(prev_solution)) {
          2
        } else {
          max(2, prev_solution - 10)   # small buffer
        }
        
        result <- find_required_N_by_search(
          J = J,
          x = x,
          methods = method_list(),
          target_pct = target_pct,
          coarse_step = 5,
          N_min = start_N,
          N_max = max(100, ceiling(2000 / J))
        )
        
        # update warm start (use GEE-Ordinal or first method)
        prev_solution <<- result$GEE_Ordinal
        
        result
      })
      
      df <- do.call(rbind, res_list)
    }
    
    df
  })
  
  
 
  
  # -----------------------------
  # Main plot
  # -----------------------------
  output$main_power_plot <- renderPlotly({
    df <- power_df()
    shiny::validate(need(!is.null(df), "Unable to generate plot."))
    
    cols <- list(
      Ordinal = "#78C2AD",
      Bin1    = "#5FA4D0",
      Bin2    = "#E59A9A",
      WH_DE   = "#9A8FB4",
      Indep   = "#A0A0A0",
      WH      = "#E5B64B"
    )

    if (input$plot_setup %in% c("power_vs_N", "power_vs_J")) {
      target_pct <- round(input$target_power, 2)
      idx <- which(df$GEE_Ordinal >= target_pct)
      
      p <- plot_ly(df, x = ~x_value)
      
      p <- p %>% add_lines(
        y = ~GEE_Ordinal,
        name = "GEE-Ordinal",
        line = list(color = cols$Ordinal, width = 3),
        hovertemplate = paste0(df$x_label[1], "=%{x}<br>Power=%{y:.2f}%<extra></extra>")
      )
      
      if (!all(is.na(df$GEE_Binary1))) {
        p <- p %>% add_lines(
          y = ~GEE_Binary1,
          name = "GEE-Binary1 (0 vs 1,2)",
          line = list(color = cols$Bin1, width = 3, dash = "dash"),
          hovertemplate = paste0(df$x_label[1], "=%{x}<br>Power=%{y:.2f}%<extra></extra>")
        )
      }
      
      if (!all(is.na(df$GEE_Binary2))) {
        p <- p %>% add_lines(
          y = ~GEE_Binary2,
          name = "GEE-Binary2 (0,1 vs 2)",
          line = list(color = cols$Bin2, width = 3, dash = "longdash"),
          hovertemplate = paste0(df$x_label[1], "=%{x}<br>Power=%{y:.2f}%<extra></extra>")
        )
      }
      
      if (!all(is.na(df$WH_DE))) {
        p <- p %>% add_lines(
          y = ~WH_DE,
          name = "Whitehead-Design Effect",
          line = list(color = cols$WH_DE, width = 3, dash = "dashdot"),
          hovertemplate = paste0(df$x_label[1], "=%{x}<br>Power=%{y:.2f}%<extra></extra>")
        )
      }
      
      if (!all(is.na(df$Whitehead))) {
        p <- p %>% add_lines(
          y = ~Whitehead,
          name = "Whitehead",
          line = list(color = cols$WH, width = 3, dash = "dashdot"),
          hovertemplate = paste0(df$x_label[1], "=%{x}<br>Power=%{y:.2f}%<extra></extra>")
        )
      }
      
      if (!all(is.na(df$GEE_Indep))) {
        p <- p %>% add_lines(
          y = ~GEE_Indep,
          name = "GEE-Indep",
          line = list(color = cols$Indep, width = 3, dash = "longdash"),
          hovertemplate = paste0(df$x_label[1], "=%{x}<br>Power=%{y:.2f}%<extra></extra>")
        )
      }
      
      p <- p %>% add_lines(
        x = df$x_value,
        y = rep(target_pct, nrow(df)),
        name = sprintf("Target = %.2f%%", target_pct),
        line = list(color = "#717171", dash = "dash", width = 2),
        hoverinfo = "skip"
      )
      
      if (length(idx) > 0) {
        i0 <- min(idx)
        p <- p %>% add_markers(
          x = df$x_value[i0],
          y = df$GEE_Ordinal[i0],
          name = "Threshold point",
          marker = list(color = cols$Ordinal, size = 10),
          hovertemplate = paste0(df$x_label[1], "=%{x}<br>Power=%{y:.2f}%<extra></extra>")
        )
      }
      
      y_max <- 100
      y_min <- 0
      
      # choose break spacing (you can tweak this)
      y_dtick <- 10   # or 5 if you want denser
      
      x_min <- min(df$x_value, na.rm = TRUE)
      x_max <- max(df$x_value, na.rm = TRUE)
      
      p %>% layout(
        hovermode = "x unified",
        title = list(
          text = title_text(),
          x = 0.02,
          xanchor = "left"
        ),
        margin = list(t = 70),   # increase top margin
        xaxis = list(
          title = df$x_label[1],
          showline = TRUE,
          mirror = FALSE,
          linecolor = "black",
          linewidth = 1,
          ticks = "outside",
          tickmode = "linear",
          dtick = if (input$plot_setup == "power_vs_N") 5 else 10
        ),
        
        yaxis = list(
          title = "Power (%)",
          range = c(y_min, y_max),
          showline = TRUE,        # ← this gives you the left axis line
          mirror = FALSE,
          linecolor = "black",
          linewidth = 1,
          ticks = "outside",
          tickmode = "linear",
          tick0 = 0,
          dtick = y_dtick,
          zeroline = FALSE
        ),
        
        legend = list(orientation = "h", x = 0, y = -0.25)
      )
    } else {
      methods <- method_list()
      
      df_long <- pivot_longer(
        df,
        cols = all_of(methods),
        names_to = "method",
        values_to = "N_required"
      )
      
      method_label_map <- c(
        GEE_Ordinal = "GEE-Ordinal",
        GEE_Binary1 = "GEE-Binary1 (0 vs 1,2)",
        GEE_Binary2 = "GEE-Binary2 (0,1 vs 2)",
        WH_DE       = "Whitehead-Design Effect",
        GEE_Indep   = "GEE-Indep",
        Whitehead   = "Whitehead"
      )
      
      df_long$method_label <- unname(method_label_map[df_long$method])
      
      color_map <- c(
        GEE_Ordinal = cols$Ordinal,
        GEE_Binary1 = cols$Bin1,
        GEE_Binary2 = cols$Bin2,
        WH_DE       = cols$WH_DE,
        GEE_Indep   = cols$Indep,
        Whitehead   = cols$WH
      )
      
      dash_map <- c(
        GEE_Ordinal = "solid",
        GEE_Binary1 = "dash",
        GEE_Binary2 = "longdash",
        WH_DE       = "dashdot",
        GEE_Indep   = "longdash",
        Whitehead   = "dashdot"
      )
      
      p <- plot_ly()
      
      for (m in methods) {
        subdf <- df_long[df_long$method == m, , drop = FALSE]
        
        if (!all(is.na(subdf$N_required))) {
          p <- p %>% add_lines(
            data = subdf,
            x = ~J,
            y = ~N_required,
            name = method_label_map[[m]],
            line = list(color = color_map[[m]], width = 3, dash = dash_map[[m]]),
            hovertemplate = "Cluster size=%{x}<br>Required total clusters=%{y}<extra></extra>"
          )
        }
      }
      
      x_min <- min(df$J, na.rm = TRUE)
      x_max <- max(df$J, na.rm = TRUE)
      y_max <- max(df_long$N_required, na.rm = TRUE)
      
      p %>% layout(
        hovermode = "x unified",
        title = list(
          text = title_text(),
          x = 0.02,
          xanchor = "left"
        ),
        margin = list(t = 70),   # increase top margin
        xaxis = list(
          title = "Cluster size",
          showline = TRUE,
          mirror = FALSE,
          linecolor = "black",
          linewidth = 1,
          ticks = "outside",
          tickmode = "linear",
          tick0 = floor(x_min / 10) * 10,
          dtick = 10
        ),
        yaxis = list(
          title = "Required total number of clusters",
          showline = TRUE,
          mirror = FALSE,
          linecolor = "black",
          linewidth = 1,
          ticks = "outside",
          tickmode = "linear",
          tick0 = 0,
          dtick = 5,
          rangemode = "tozero"
        ),
        legend = list(orientation = "h", x = 0, y = -0.25)
      )
    }
    
  })

  # -----------------------------
  # Summary text
  # -----------------------------
  output$main_plot_summary <- renderText({
    df <- power_df()
    target_pct <- round(input$target_power, 2)
    
    if (input$plot_setup == "power_vs_N") {
      idx <- which(df$GEE_Ordinal >= target_pct)
      if (length(idx) == 0) {
        sprintf("GEE-Ordinal: target %.2f%% not reached in the selected range.", target_pct)
      } else {
        sprintf(
          "GEE-Ordinal: minimum total number of clusters to reach target %.2f%% is N = %d.",
          target_pct, df$N[min(idx)]
        )
      }
    } else if (input$plot_setup == "power_vs_J") {
      idx <- which(df$GEE_Ordinal >= target_pct)
      if (length(idx) == 0) {
        sprintf("GEE-Ordinal: target %.2f%% not reached in the selected range.", target_pct)
      } else {
        sprintf(
          "GEE-Ordinal: minimum cluster size to reach target %.2f%% is J = %d.",
          target_pct, df$J[min(idx)]
        )
      }
    } else {
      method_cols <- method_list()
      all_na <- all(sapply(method_cols, function(m) all(is.na(df[[m]]))))
      
      if (all_na) {
        sprintf("Target %.2f%% was not reached within the internal search range of total clusters.", target_pct)
      } else {
        sprintf(
          "For each cluster size, the curves show the minimum total number of clusters needed to reach %.2f%% power for each selected method.",
          target_pct
        )
      }
      
    }
  })
  
  output$whitehead_crt_note <- renderUI({
    req(#input$plot_setup == "N_vs_J", 
        input$add_indep)
    
    tags$div(
      style = "
      margin-top: 8px;
      font-size: 13px;
      color: #5D6D7E;
      line-height: 1.45;
      background: #fbfcfc;
      border-left: 4px solid #E5B64B;
      padding: 8px 10px;
      border-radius: 6px;
    ",
      HTML(
        "Note: Whitehead's formula gives the total individual-level sample size, 
      which corresponds to \\(N \\times J\\) in the CRT setting, where \\(N\\) is the 
     total number of clusters and \\(J\\) is the cluster size."
      )
    )
  })
  # -----------------------------
  # Download
  # -----------------------------
  output$download_power_csv <- downloadHandler(
    filename = function() {
      paste0("power_table_", Sys.Date(), ".csv")
    },
    content = function(file) {
      df <- power_df()
      shiny::validate(need(!is.null(df) && nrow(df) > 0, "No results to download yet."))
      write.csv(df, file, row.names = FALSE)
    }
  )
  
  # ============================================================
  # CRT 4+ Cat. tab: collapse original K-category outcome to 3 categories
  # ============================================================
  
  # -----------------------------
  # Parse original K-category distribution
  # -----------------------------
  crtK_dist_vec <- reactive({
    parse_crtK_dist_text(
      input$crtK_dist_input,
      if (input$crtK_dist_source == "control") {
        "Control distribution"
      } else {
        "Overall distribution"
      }
    )
  })
  
  crtK_theta_val <- reactive({
    or_val <- as.numeric(trimws(input$crtK_OR_text))
    
    shiny::validate(
      need(
        length(or_val) == 1 && is.finite(or_val) && or_val > 0,
        "Treatment effect (OR) must be a positive number."
      )
    )
    
    log(or_val)
  })
  
  crtK_resolved_original <- reactive({
    shiny::validate(
      need(
        exists("resolve_ordinal_inputs"),
        "resolve_ordinal_inputs() not found. Did you source your R/ files?"
      ),
      need(
        is.function(resolve_ordinal_inputs),
        "resolve_ordinal_inputs exists but is not a function."
      )
    )
    
    if (input$crtK_dist_source == "control") {
      resolve_ordinal_inputs(
        p_C_vec = crtK_dist_vec(),
        treatment_effect = crtK_theta_val()
      )
    } else {
      resolve_ordinal_inputs(
        pi_vec = crtK_dist_vec(),
        treatment_effect = crtK_theta_val()
      )
    }
  })
  
  # -----------------------------
  # Dynamic 3-category collapse inputs
  # -----------------------------
  output$crtK_collapse_inputs <- renderUI({
    x <- crtK_resolved_original()
    K <- x$K
    n_way <- input$crtK_n_collapse
    
    tagList(
      lapply(seq_len(n_way), function(i) {
        textInput(
          inputId = paste0("crtK_collapse_", i),
          label = paste0("3-category collapsing scheme ", i),
          value = if (i == 1) default_3cat_collapse_text(K) else "",
          placeholder = "e.g., 1,2; 3,4; 5,6"
        )
      })
    )
  })
  
  # -----------------------------
  # Dynamic binary dichotomization inputs
  # -----------------------------
  output$crtK_binary_inputs <- renderUI({
    req(input$crtK_add_bin)
    
    x <- crtK_resolved_original()
    K <- x$K
    n_way <- input$crtK_n_binary
    
    tagList(
      lapply(seq_len(n_way), function(i) {
        textInput(
          inputId = paste0("crtK_binary_", i),
          label = paste0("Binary dichotomization ", i),
          value = if (i == 1) default_binary_collapse_text(K) else "",
          placeholder = "e.g., 1; 2,3,4,5,6"
        )
      })
    )
  })
  
  # -----------------------------
  # 3-category ordinal collapsing schemes
  # -----------------------------
  crtK_ordinal_schemes <- reactive({
    x <- crtK_resolved_original()
    K <- x$K
    n_way <- input$crtK_n_collapse
    
    out <- lapply(seq_len(n_way), function(i) {
      txt <- input[[paste0("crtK_collapse_", i)]]
      
      groups <- tryCatch(
        parse_collapse_text(txt, K = K, G = 3),
        error = function(e) {
          shiny::validate(
            need(FALSE, paste0("3-category collapsing scheme ", i, ": ", e$message))
          )
        }
      )
      
      p_C_3 <- collapse_probs_by_groups(x$p_C_vec, groups)
      p_T_3 <- collapse_probs_by_groups(x$p_T_vec, groups)
      pi_3  <- collapse_probs_by_groups(x$pi_vec, groups)
      
      theta_3 <- theta_from_probs_K(p_C_3, p_T_3)
      
      list(
        scheme_id = paste0("Ord-", i),
        scheme_no = i,
        type = "Ordinal",
        text = txt,
        groups = groups,
        label = collapse_label_from_groups(groups),
        p_C_vec = p_C_3,
        p_T_vec = p_T_3,
        pi_vec = pi_3,
        treatment_effect = theta_3
      )
    })
    
    out
  })
  
  # -----------------------------
  # Binary dichotomization schemes
  # -----------------------------
  crtK_binary_schemes <- reactive({
    if (!isTRUE(input$crtK_add_bin)) {
      return(list())
    }
    
    x <- crtK_resolved_original()
    K <- x$K
    n_way <- input$crtK_n_binary
    
    out <- lapply(seq_len(n_way), function(i) {
      txt <- input[[paste0("crtK_binary_", i)]]
      
      groups <- tryCatch(
        parse_collapse_text(txt, K = K, G = 2),
        error = function(e) {
          shiny::validate(
            need(FALSE, paste0("Binary dichotomization ", i, ": ", e$message))
          )
        }
      )
      
      p_C_2 <- collapse_probs_by_groups(x$p_C_vec, groups)
      p_T_2 <- collapse_probs_by_groups(x$p_T_vec, groups)
      pi_2  <- collapse_probs_by_groups(x$pi_vec, groups)
      
      list(
        scheme_id = paste0("Bin-", i),
        scheme_no = i,
        type = "Binary",
        text = txt,
        groups = groups,
        label = collapse_label_from_groups(groups),
        p_C_vec = p_C_2,
        p_T_vec = p_T_2,
        pi_vec = pi_2,
        treatment_effect = x$treatment_effect
      )
    })
    
    out
  })
  
  # -----------------------------
  # Power calculation helpers
  # -----------------------------
  crtK_calc_power_ordinal_once <- function(N, J, scheme) {
    
    shiny::validate(
      need(
        exists("calculate_power_3cat"),
        "calculate_power_3cat() not found. Did you source your R/ files?"
      ),
      need(
        is.function(calculate_power_3cat),
        "calculate_power_3cat exists but is not a function."
      )
    )
    
    out <- calculate_power_3cat(
      p_C_vec = scheme$p_C_vec,
      treatment_effect = scheme$treatment_effect,
      number_cluster = N,
      cluster_size = J,
      rho11 = input$crtK_rho11,
      rho12 = input$crtK_rho12,
      rho22 = input$crtK_rho22,
      add_bin = FALSE,
      bin_rho1 = NULL,
      bin_rho2 = NULL,
      add_WH_DE = input$crtK_add_WH_DE,
      WH_DE_ICC = input$crtK_WH_DE_ICC,
      add_indep = FALSE,
      alpha = input$crtK_alpha
    )
    
    to_pct <- function(z) {
      if (is.null(z) || length(z) == 0 || is.na(z)) return(NA_real_)
      round(100 * as.numeric(z), 2)
    }
    
    rows <- list(
      data.frame(
        N = as.integer(N),
        J = as.integer(J),
        scheme_id = scheme$scheme_id,
        scheme_label = scheme$label,
        method = "GEE_Ordinal",
        method_label = paste0(scheme$scheme_id, ": GEE-Ordinal"),
        power = to_pct(out$power_GEE_exch),
        check.names = FALSE
      )
    )
    
    if (isTRUE(input$crtK_add_WH_DE)) {
      rows[[length(rows) + 1]] <- data.frame(
        N = as.integer(N),
        J = as.integer(J),
        scheme_id = scheme$scheme_id,
        scheme_label = scheme$label,
        method = "WH_DE",
        method_label = paste0(scheme$scheme_id, ": Whitehead-DE"),
        power = to_pct(out$power_WH_DE),
        check.names = FALSE
      )
    }
    
    do.call(rbind, rows)
  }
  
  crtK_calc_power_binary_once <- function(N, J, scheme) {
    
    pwr <- power_binary_crt_custom(
      p_C_event = scheme$p_C_vec[1],
      p_T_event = scheme$p_T_vec[1],
      theta_R = scheme$treatment_effect,
      rho = input$crtK_binary_ICC,
      cluster_size = J,
      number_cluster = N,
      alpha = input$crtK_alpha
    )
    
    data.frame(
      N = as.integer(N),
      J = as.integer(J),
      scheme_id = scheme$scheme_id,
      scheme_label = scheme$label,
      method = "GEE_Binary",
      method_label = paste0(scheme$scheme_id, ": GEE-Binary"),
      power = round(100 * as.numeric(pwr), 2),
      check.names = FALSE
    )
  }
  
  crtK_calc_power_once <- function(N, J) {
    
    ord_schemes <- crtK_ordinal_schemes()
    bin_schemes <- crtK_binary_schemes()
    
    ord_rows <- lapply(ord_schemes, function(s) {
      crtK_calc_power_ordinal_once(N = N, J = J, scheme = s)
    })
    
    bin_rows <- if (isTRUE(input$crtK_add_bin) && length(bin_schemes) > 0) {
      lapply(bin_schemes, function(s) {
        crtK_calc_power_binary_once(N = N, J = J, scheme = s)
      })
    } else {
      list()
    }
    
    do.call(rbind, c(ord_rows, bin_rows))
  }
  
  # -----------------------------
  # Required N search
  # -----------------------------
  crtK_find_required_N_by_search <- function(J,
                                             target_pct,
                                             coarse_step = 5,
                                             N_min = 2,
                                             N_max = max(100, ceiling(2000 / J))) {
    
    scheme_rows <- crtK_calc_power_once(N = N_min, J = J)
    
    combos <- unique(
      scheme_rows[, c("scheme_id", "scheme_label", "method", "method_label")]
    )
    
    found <- combos
    found$N_required <- NA_real_
    
    coarse_grid <- seq(N_min, N_max, by = coarse_step)
    if (tail(coarse_grid, 1) != N_max) {
      coarse_grid <- c(coarse_grid, N_max)
    }
    
    last_N <- NA_integer_
    
    for (N in coarse_grid) {
      rows <- crtK_calc_power_once(N = N, J = J)
      
      for (i in seq_len(nrow(found))) {
        if (is.na(found$N_required[i])) {
          hit <- rows[
            rows$scheme_id == found$scheme_id[i] &
              rows$method == found$method[i] &
              !is.na(rows$power) &
              rows$power >= target_pct,
          ]
          
          if (nrow(hit) > 0) {
            lower_N <- if (is.na(last_N)) N_min else max(N_min, last_N + 1)
            
            for (N_ref in seq(lower_N, N, by = 1)) {
              rows_ref <- crtK_calc_power_once(N = N_ref, J = J)
              
              hit_ref <- rows_ref[
                rows_ref$scheme_id == found$scheme_id[i] &
                  rows_ref$method == found$method[i] &
                  !is.na(rows_ref$power) &
                  rows_ref$power >= target_pct,
              ]
              
              if (nrow(hit_ref) > 0) {
                found$N_required[i] <- N_ref
                break
              }
            }
          }
        }
      }
      
      if (all(!is.na(found$N_required))) break
      
      last_N <- N
    }
    
    found$J <- J
    found
  }
  
  
  
  # -----------------------------
  # CRT 4+ plot data
  # -----------------------------
  crtK_power_df <- reactive({
    
    if (input$crtK_plot_setup == "power_vs_N") {
      res_list <- lapply(
        seq(input$crtK_N_range[1], input$crtK_N_range[2], by = 1),
        function(N) {
          df <- crtK_calc_power_once(N = N, J = input$crtK_J_fixed)
          df$x_value <- df$N
          df$x_label <- "Total number of clusters"
          df
        }
      )
      
      do.call(rbind, res_list)
      
    } else if (input$crtK_plot_setup == "power_vs_J") {
      res_list <- lapply(
        seq(input$crtK_J_range[1], input$crtK_J_range[2], by = 1),
        function(J) {
          df <- crtK_calc_power_once(N = input$crtK_N_fixed, J = J)
          df$x_value <- df$J
          df$x_label <- "Cluster size"
          df
        }
      )
      
      do.call(rbind, res_list)
      
    } else {
      target_pct <- round(input$crtK_target_power, 2)
      
      res_list <- lapply(
        seq(input$crtK_J_range_req[1], input$crtK_J_range_req[2], by = 1),
        function(J) {
          crtK_find_required_N_by_search(
            J = J,
            target_pct = target_pct,
            coarse_step = 5,
            N_min = 2,
            N_max = max(100, ceiling(2000 / J))
          )
        }
      )
      
      do.call(rbind, res_list)
    }
  })
  
  # -----------------------------
  # CRT 4+ resolved display
  # -----------------------------
  output$crtK_resolved_inputs <- renderUI({
    x <- crtK_resolved_original()
    ord_schemes <- crtK_ordinal_schemes()
    bin_schemes <- crtK_binary_schemes()
    
    fmt_vec <- function(v, digits = 4) {
      paste0("(", paste(sprintf(paste0("%.", digits, "f"), v), collapse = ", "), ")")
    }
    
    fmt_num <- function(v, digits = 4) {
      sprintf(paste0("%.", digits, "f"), v)
    }
    
    ordinal_blocks <- lapply(ord_schemes, function(s) {
      tags$div(
        style = "margin-top:8px;",
        tags$span(style = "font-weight:600; color:#34495E;", paste0(s$scheme_id, ": ")),
        tags$span(paste0(s$label, "  |  collapsed overall = ", fmt_vec(s$pi_vec))),
        tags$br(),
        tags$span(style = "font-size:13px; color:#7a7a7a;",
                  paste0("Control = ", fmt_vec(s$p_C_vec),
                         " | Treatment = ", fmt_vec(s$p_T_vec),
                         " | log OR = ", fmt_num(s$treatment_effect)))
      )
    })
    
    binary_blocks <- if (isTRUE(input$crtK_add_bin) && length(bin_schemes) > 0) {
      tagList(
        tags$div(
          style = "margin-top:14px; padding-top:10px; border-top:1px solid #D9DEE3; font-weight:700; color:#34495E;",
          "Binary dichotomizations"
        ),
        tags$div(
          style = "font-size:13px; color:#7a7a7a; margin-bottom:6px;",
          "For GEE-Binary, the user-specified odds ratio is applied to each selected binary dichotomization."
        ),
        lapply(bin_schemes, function(s) {
          tags$div(
            style = "margin-top:8px;",
            tags$span(style = "font-weight:600; color:#34495E;", paste0(s$scheme_id, ": ")),
            tags$span(paste0(s$label, "  |  binary overall = ", fmt_vec(s$pi_vec))),
            tags$br(),
            tags$span(style = "font-size:13px; color:#7a7a7a;",
                      paste0("Control = ", fmt_vec(s$p_C_vec),
                             " | Treatment = ", fmt_vec(s$p_T_vec),
                             " | assumed log OR = ", fmt_num(s$treatment_effect)))
          )
        })
      )
    }
    
    tags$div(
      style = "
      font-size: 15px;
      color: #4a4a4a;
      line-height: 1.45;
      margin-top: 6px;
    ",
      
      tags$div(
        style = "margin-bottom: 8px;",
        tags$span(style = "font-weight:600; color:#34495E;", "Original number of categories:"),
        tags$span(paste0(" ", x$K))
      ),
      
      tags$div(
        style = "margin-bottom: 8px;",
        tags$span(style = "font-weight:600; color:#34495E;", "Original control:"),
        tags$span(paste0(" ", fmt_vec(x$p_C_vec)))
      ),
      
      tags$div(
        style = "margin-bottom: 8px;",
        tags$span(style = "font-weight:600; color:#34495E;", "Original treatment:"),
        tags$span(paste0(" ", fmt_vec(x$p_T_vec)))
      ),
      
      tags$div(
        style = "margin-bottom: 8px;",
        tags$span(style = "font-weight:600; color:#34495E;", "Original overall:"),
        tags$span(paste0(" ", fmt_vec(x$pi_vec)))
      ),
      
      tags$div(
        style = "margin-bottom: 8px;",
        tags$span(style = "font-weight:600; color:#34495E;", "Treatment effect (OR):"),
        tags$span(
          HTML(paste0(
            " ", fmt_num(exp(x$treatment_effect)),
            " &nbsp;&nbsp; (log = ", fmt_num(x$treatment_effect), ")"
          ))
        )
      ),
      
      tags$div(
        style = "margin-top:14px; padding-top:10px; border-top:1px solid #D9DEE3; font-weight:700; color:#34495E;",
        "3-category collapsing schemes used for ordinal power calculation"
      ),
      
      ordinal_blocks,
      
      tags$div(
        style = "margin-top:14px; padding-top:10px; border-top:1px solid #D9DEE3;",
        tags$span(style = "font-weight:600; color:#34495E;", "Collapsed 3-category within-cluster correlation matrix: "),
        HTML(
          paste0(
            "&rho;<sub>11</sub> = ", fmt_num(input$crtK_rho11),
            ", &rho;<sub>12</sub> = ", fmt_num(input$crtK_rho12),
            ", &rho;<sub>22</sub> = ", fmt_num(input$crtK_rho22)
          )
        )
      ),
      
      if (isTRUE(input$crtK_add_WH_DE)) {
        tags$div(
          style = "margin-top:8px;",
          tags$span(style = "font-weight:600; color:#34495E;", "Whitehead-DE ICC: "),
          tags$span(fmt_num(input$crtK_WH_DE_ICC))
        )
      },
      
      if (isTRUE(input$crtK_add_bin)) {
        tags$div(
          style = "margin-top:8px;",
          tags$span(style = "font-weight:600; color:#34495E;", "Binary ICC: "),
          tags$span(fmt_num(input$crtK_binary_ICC))
        )
      },
      
      binary_blocks
    )
  })
  
  # -----------------------------
  # CRT 4+ plot
  # -----------------------------
  get_crtK_style <- function(method_label) {
    
    # Extract Ord number if present: Ord-1, Ord1, Ord-2, Ord2, ...
    ord_match <- regmatches(method_label, regexpr("Ord-?[0-9]+", method_label))
    ord_num <- NA_integer_
    if (length(ord_match) > 0 && nzchar(ord_match)) {
      ord_num <- as.integer(gsub("[^0-9]", "", ord_match))
    }
    
    # Extract binary number if present: Bin-1, Bin1, Binary-1, Binary1, ...
    bin_match <- regmatches(method_label, regexpr("Binary-?[0-9]+|Bin-?[0-9]+", method_label))
    bin_num <- NA_integer_
    if (length(bin_match) > 0 && nzchar(bin_match)) {
      bin_num <- as.integer(gsub("[^0-9]", "", bin_match))
    }
    
    # -----------------------------
    # GEE-Ordinal: green/teal family
    # Ord-1 matches CRT (3 Cat.) GEE-Ordinal mint
    # -----------------------------
    if (!is.na(ord_num) && grepl("GEE-Ordinal", method_label)) {
      ord_colors <- c(
        "#78C2AD",  # Ord-1: mint, same as CRT (3 Cat.)
        "#2A9D8F",  # Ord-2: darker teal
        "#4DB6AC",  # Ord-3: teal
        "#A3D9C9",  # Ord-4: light mint
        "#3D8B7D"   # Ord-5: muted dark teal
      )
      
      idx <- ((ord_num - 1) %% length(ord_colors)) + 1
      
      return(list(
        color = ord_colors[idx],
        dash  = "solid"
      ))
    }
    
    # -----------------------------
    # Whitehead-DE: purple family
    # Ord-1 matches CRT (3 Cat.) Whitehead-DE purple
    # Different Ord schemes get slightly different purple shades and dash patterns
    # -----------------------------
    if (!is.na(ord_num) && grepl("Whitehead-DE|Whitehead-Design Effect|WH-DE", method_label)) {
      wh_colors <- c(
        "#9A8FB4",  # Ord-1: original purple
        "#7E6FA8",  # Ord-2: deeper purple
        "#B39DDB",  # Ord-3: lighter purple
        "#6D597A",  # Ord-4: muted plum
        "#C8B6E2"   # Ord-5: pale lavender
      )
      
      wh_dashes <- c(
        "dashdot",
        "dot",
        "dash",
        "longdash",
        "solid"
      )
      
      idx <- ((ord_num - 1) %% length(wh_colors)) + 1
      
      return(list(
        color = wh_colors[idx],
        dash  = wh_dashes[idx]
      ))
    }
    
    # -----------------------------
    # GEE-Binary: keep Bin-1 and Bin-2 consistent with CRT (3 Cat.)
    # -----------------------------
    if (!is.na(bin_num) && bin_num == 1) {
      return(list(
        color = "#5FA4D0",  # blue, same as CRT (3 Cat.) GEE-Binary1
        dash  = "dash"
      ))
    }
    
    if (!is.na(bin_num) && bin_num == 2) {
      return(list(
        color = "#E59A9A",  # pink, same as CRT (3 Cat.) GEE-Binary2
        dash  = "longdash"
      ))
    }
    
    # -----------------------------
    # Additional binary dichotomizations, if any
    # -----------------------------
    if (!is.na(bin_num) && bin_num >= 3) {
      extra_bin_colors <- c(
        "#C77DFF",  # Bin-3
        "#F4A261",  # Bin-4
        "#B56576",  # Bin-5
        "#6D597A"   # extra
      )
      
      extra_bin_dashes <- c(
        "dot",
        "dash",
        "longdash",
        "dashdot"
      )
      
      idx <- ((bin_num - 3) %% length(extra_bin_colors)) + 1
      
      return(list(
        color = extra_bin_colors[idx],
        dash  = extra_bin_dashes[idx]
      ))
    }
    
    # -----------------------------
    # Target line or fallback
    # -----------------------------
    if (grepl("Target", method_label)) {
      return(list(
        color = "#717171",
        dash  = "dash"
      ))
    }
    
    list(
      color = "#666666",
      dash  = "solid"
    )
  }
  output$crtK_power_plot <- renderPlotly({
    df <- crtK_power_df()
    target_pct <- round(input$crtK_target_power, 2)
    
    shiny::validate(
      need(!is.null(df) && nrow(df) > 0, "Unable to generate plot.")
    )
    
    if (input$crtK_plot_setup %in% c("power_vs_N", "power_vs_J")) {
      
      p <- plot_ly()
      
      curve_ids <- unique(df$method_label)
      
      for (cid in curve_ids) {
        subdf <- df[df$method_label == cid, , drop = FALSE]
        
        sty <- get_crtK_style(cid)
        
        p <- p %>% add_lines(
          data = subdf,
          x = ~x_value,
          y = ~power,
          name = cid,
          line = list(color = sty$color, width = 3, dash = sty$dash),
          hovertemplate = paste0(
            subdf$x_label[1],
            "=%{x}<br>Power=%{y:.2f}%<extra></extra>"
          )
        )
      }
      
      p <- p %>% add_lines(
        x = unique(df$x_value),
        y = rep(target_pct, length(unique(df$x_value))),
        name = sprintf("Target = %.2f%%", target_pct),
        line = list(color = "#717171", dash = "dash", width = 2),
        hoverinfo = "skip"
      )
      
      p %>% layout(
        hovermode = "x unified",
        title = list(
          text = if (input$crtK_plot_setup == "power_vs_N") {
            paste0(
              "Power vs Total Number of Clusters",
              " (Collapsed/Dichotomized Outcomes; Cluster size = ",
              input$crtK_J_fixed,
              " | Target power = ",
              target_pct,
              "%)"
            )
          } else {
            paste0(
              "Power vs Cluster Size",
              " (Collapsed/Dichotomized Outcomes; Total clusters = ",
              input$crtK_N_fixed,
              " | Target power = ",
              target_pct,
              "%)"
            )
          },
          x = 0.02,
          xanchor = "left"
        ),
        margin = list(t = 70),
        xaxis = list(
          title = df$x_label[1],
          showline = TRUE,
          mirror = FALSE,
          linecolor = "black",
          linewidth = 1,
          ticks = "outside",
          tickmode = "linear",
          dtick = if (input$crtK_plot_setup == "power_vs_N") 5 else 10
        ),
        yaxis = list(
          title = "Power (%)",
          range = c(0, 100),
          showline = TRUE,
          mirror = FALSE,
          linecolor = "black",
          linewidth = 1,
          ticks = "outside",
          tickmode = "linear",
          tick0 = 0,
          dtick = 10,
          zeroline = FALSE
        ),
        legend = list(orientation = "h", x = 0, y = -0.30)
      )
      
    } else {
      
      p <- plot_ly()
      
      curve_ids <- unique(df$method_label)
      
      for (cid in curve_ids) {
        subdf <- df[df$method_label == cid, , drop = FALSE]
        
        sty <- get_crtK_style(cid)
        
        p <- p %>% add_lines(
          data = subdf,
          x = ~J,
          y = ~N_required,
          name = cid,
          line = list(color = sty$color, width = 3, dash = sty$dash),
          hovertemplate = "Cluster size=%{x}<br>Required total clusters=%{y:.2f}<extra></extra>"
        )
      }
      
      p %>% layout(
        hovermode = "x unified",
        title = list(
          text = paste0(
            "Required Number of Clusters vs Cluster Size",
            " (Collapsed/Dichotomized Outcomes; Target power = ",
            target_pct,
            "%)"
          ),
          x = 0.02,
          xanchor = "left"
        ),
        margin = list(t = 70),
        xaxis = list(
          title = "Cluster size",
          showline = TRUE,
          mirror = FALSE,
          linecolor = "black",
          linewidth = 1,
          ticks = "outside",
          tickmode = "linear",
          dtick = 10
        ),
        yaxis = list(
          title = "Required total number of clusters",
          showline = TRUE,
          mirror = FALSE,
          linecolor = "black",
          linewidth = 1,
          ticks = "outside",
          tickmode = "linear",
          tick0 = 0,
          dtick = 5,
          rangemode = "tozero"
        ),
        legend = list(orientation = "h", x = 0, y = -0.30)
      )
    }
  })
  
  # -----------------------------
  # CRT 4+ summary text
  # -----------------------------
  output$crtK_plot_summary <- renderText({
    df <- crtK_power_df()
    target_pct <- round(input$crtK_target_power, 2)
    
    if (input$crtK_plot_setup %in% c("power_vs_N", "power_vs_J")) {
      ord_df <- df[df$method == "GEE_Ordinal", , drop = FALSE]
      
      if (nrow(ord_df) == 0) {
        return("No GEE-Ordinal results are available.")
      }
      
      summaries <- lapply(split(ord_df, ord_df$scheme_id), function(d) {
        idx <- which(d$power >= target_pct)
        scheme <- unique(d$scheme_id)
        
        if (length(idx) == 0) {
          paste0(scheme, " GEE-Ordinal: target power not reached in the selected range")
        } else {
          if (input$crtK_plot_setup == "power_vs_N") {
            paste0(
              scheme,
              " GEE-Ordinal: minimum total number of clusters to reach the target power is N = ",
              d$N[min(idx)]
            )
          } else {
            paste0(
              scheme,
              " GEE-Ordinal: minimum cluster size to reach the target power is J = ",
              d$J[min(idx)]
            )
          }
        }
      })
      
      paste0(paste(unlist(summaries), collapse = "; "), ".")
      
    } else {
      paste0(
        "For each cluster size, the curves show the minimum total number of clusters needed to reach ",
        target_pct,
        "% power for each selected collapsing/dichotomization scheme and method."
      )
    }
  })
  
  # -----------------------------
  # CRT 4+ collapse summary table
  # -----------------------------
  output$crtK_collapse_table <- renderTable({
    ord_schemes <- crtK_ordinal_schemes()
    bin_schemes <- crtK_binary_schemes()
    
    fmt_vec <- function(v, digits = 4) {
      paste0("(", paste(sprintf(paste0("%.", digits, "f"), v), collapse = ", "), ")")
    }
    
    ord_tbl <- do.call(rbind, lapply(ord_schemes, function(s) {
      data.frame(
        Type = "3-category ordinal collapse",
        ID = s$scheme_id,
        Collapse = s$label,
        `Collapsed overall distribution` = fmt_vec(s$pi_vec),
        `Assumed log OR` = sprintf("%.4f", s$treatment_effect),
        check.names = FALSE
      )
    }))
    
    bin_tbl <- if (isTRUE(input$crtK_add_bin) && length(bin_schemes) > 0) {
      do.call(rbind, lapply(bin_schemes, function(s) {
        data.frame(
          Type = "Binary dichotomization",
          ID = s$scheme_id,
          Collapse = s$label,
          `Collapsed overall distribution` = fmt_vec(s$pi_vec),
          `Assumed log OR` = sprintf("%.4f", s$treatment_effect),
          check.names = FALSE
        )
      }))
    } else {
      NULL
    }
    
    rbind(ord_tbl, bin_tbl)
  })
  
  # -----------------------------
  # CRT 4+ download
  # -----------------------------
  output$download_crtK_power_csv <- downloadHandler(
    filename = function() {
      paste0("crt4plus_power_table_", Sys.Date(), ".csv")
    },
    content = function(file) {
      df <- crtK_power_df()
      shiny::validate(
        need(!is.null(df) && nrow(df) > 0, "No CRT 4+ results to download yet.")
      )
      write.csv(df, file, row.names = FALSE)
    }
  )
  # ============================================================
  # IRT tab: Whitehead formula for general ordinal outcomes
  # ============================================================
  
  # -----------------------------
  # Parse IRT ordinal distribution
  # -----------------------------
  irt_dist_vec <- reactive({
    parse_ordinal_dist_text(
      input$irt_dist_input,
      if (input$irt_dist_source == "control") {
        "Control distribution"
      } else {
        "Overall distribution"
      }
    )
  })
  
  # -----------------------------
  # Parse IRT treatment effect
  # -----------------------------
  irt_theta_val <- reactive({
    or_val <- as.numeric(trimws(input$irt_OR_text))
    
    shiny::validate(
      need(
        length(or_val) == 1 && is.finite(or_val) && or_val > 0,
        "Treatment effect (OR) must be a positive number."
      )
    )
    
    log(or_val)
  })
  
  # -----------------------------
  # Resolve IRT inputs
  # -----------------------------
  irt_resolved_inputs <- reactive({
    shiny::validate(
      need(
        exists("resolve_ordinal_inputs"),
        "resolve_ordinal_inputs() not found. Did you source your R/ files?"
      ),
      need(
        is.function(resolve_ordinal_inputs),
        "resolve_ordinal_inputs exists but is not a function."
      )
    )
    
    if (input$irt_dist_source == "control") {
      resolve_ordinal_inputs(
        p_C_vec = irt_dist_vec(),
        treatment_effect = irt_theta_val()
      )
    } else {
      resolve_ordinal_inputs(
        pi_vec = irt_dist_vec(),
        treatment_effect = irt_theta_val()
      )
    }
  })
  
  # -----------------------------
  # IRT plot data using wh_power()
  # -----------------------------
  irt_power_df <- reactive({
    x <- irt_resolved_inputs()
    
    shiny::validate(
      need(
        exists("wh_power"),
        "wh_power() not found. Did you source your R/ files?"
      ),
      need(
        is.function(wh_power),
        "wh_power exists but is not a function."
      )
    )
    
    n_seq <- seq(
      input$irt_n_range[1],
      input$irt_n_range[2],
      by = 1
    )
    
    power_WH <- wh_power(
      theta_R = x$treatment_effect,
      n = n_seq,
      pr = x$pi_vec,
      A = 1,
      alpha = input$irt_alpha
    )
    
    
    
    data.frame(
      n_total = n_seq,
      Whitehead_raw = 100 * power_WH,
      Whitehead = round(100 * power_WH, 2),
      target_power = round(input$irt_target_power, 2),
      check.names = FALSE
    )
  })
  
  # -----------------------------
  # IRT resolved input display
  # -----------------------------
  output$irt_resolved_inputs <- renderUI({
    x <- irt_resolved_inputs()
    
    fmt_vec <- function(v, digits = 4) {
      paste0("(", paste(sprintf(paste0("%.", digits, "f"), v), collapse = ", "), ")")
    }
    
    fmt_num <- function(v, digits = 4) {
      sprintf(paste0("%.", digits, "f"), v)
    }
    
    tags$div(
      style = "
      font-size: 15px;
      color: #4a4a4a;
      line-height: 1.45;
      margin-top: 6px;
    ",
      
      tags$div(
        style = "margin-bottom: 10px;",
        tags$span(style = "font-weight:600; color:#34495E;", "Number of ordinal categories:"),
        tags$span(paste0(" ", x$K))
      ),
      
      tags$div(
        style = "margin-bottom: 10px;",
        tags$span(style = "font-weight:600; color:#34495E;", "Control:"),
        tags$span(paste0(" ", fmt_vec(x$p_C_vec), "   |   ")),
        
        tags$span(style = "font-weight:600; color:#34495E;", "Treatment:"),
        tags$span(paste0(" ", fmt_vec(x$p_T_vec), "   |   ")),
        
        tags$span(style = "font-weight:600; color:#34495E;", "Overall:"),
        tags$span(paste0(" ", fmt_vec(x$pi_vec)))
      ),
      
      tags$div(
        style = "margin-bottom: 10px;",
        tags$span(style = "font-weight:600; color:#34495E;", "Treatment effect (OR):"),
        tags$span(
          HTML(paste0(
            " ", fmt_num(exp(x$treatment_effect)),
            " &nbsp;&nbsp; (log = ", fmt_num(x$treatment_effect), ")"
          ))
        )
      )
      
    )
  })
  
  output$irt_power_plot <- renderPlotly({
    df <- irt_power_df()
    target_pct <- round(input$irt_target_power, 2)
    
    cols <- list(
      WH = "#E5B64B",
      Target = "#717171"
    )
    
    p <- plot_ly(df, x = ~n_total)
    
    p <- p %>% add_lines(
      y = ~Whitehead,
      name = "Whitehead",
      line = list(color = cols$WH, width = 3),
      hovertemplate = "Total sample size=%{x}<br>Power=%{y:.2f}%<extra></extra>"
    )
    
    p <- p %>% add_lines(
      x = df$n_total,
      y = rep(target_pct, nrow(df)),
      name = sprintf("Target = %.2f%%", target_pct),
      line = list(color = cols$Target, dash = "dash", width = 2),
      hoverinfo = "skip"
    )
    
    idx <- which(df$Whitehead_raw >= target_pct)
    
    if (length(idx) > 0) {
      i0 <- min(idx)
      
      p <- p %>% add_markers(
        x = df$n_total[i0],
        y = df$Whitehead[i0],
        name = "Threshold point",
        marker = list(color = cols$WH, size = 10),
        hovertemplate = "Total sample size=%{x}<br>Power=%{y:.2f}%<extra></extra>"
      )
    }
    
    p %>% layout(
      hovermode = "x unified",
      title = list(
        text = paste0(
          "Power vs Total Sample Size",
          " (Target power = ", target_pct, "%)"
        ),
        x = 0.02,
        xanchor = "left"
      ),
      margin = list(t = 70),
      
      xaxis = list(
        title = "Total sample size",
        showline = TRUE,
        mirror = FALSE,
        linecolor = "black",
        linewidth = 1,
        ticks = "outside",
        tickmode = "linear",
        dtick = 50
      ),
      
      yaxis = list(
        title = "Power (%)",
        range = c(0, 100),
        showline = TRUE,
        mirror = FALSE,
        linecolor = "black",
        linewidth = 1,
        ticks = "outside",
        tickmode = "linear",
        tick0 = 0,
        dtick = 10,
        zeroline = FALSE
      ),
      
      legend = list(
        orientation = "h",
        x = 0,
        y = -0.25
      )
    )
  })
  
  output$irt_plot_summary <- renderText({
    df <- irt_power_df()
    target_pct <- round(input$irt_target_power, 2)
    
    idx <- which(df$Whitehead >= target_pct)
    
    if (length(idx) == 0) {
      sprintf(
        "Whitehead: target %.2f%% was not reached in the selected sample size range.",
        target_pct
      )
    } else {
      sprintf(
        "Whitehead: minimum total sample size to reach target %.2f%% is n = %d.",
        target_pct,
        df$n_total[min(idx)]
      )
    }
  })
  
  output$download_irt_power_csv <- downloadHandler(
    filename = function() {
      paste0("irt_power_table_", Sys.Date(), ".csv")
    },
    content = function(file) {
      df <- irt_power_df()
      shiny::validate(
        need(!is.null(df) && nrow(df) > 0, "No IRT results to download yet.")
      )
      write.csv(df, file, row.names = FALSE)
    }
  )
  
  
  # ============================================================
  # Composite
  # ============================================================
  
  # Define constants for logic
  tau_values <- c(High = .6, Moderate = .4, Low = .2, `Very low` = .05, Indep = 0)
  # Reactive logic for Composite Distribution
  composite_results <- reactive( {
    if (input$comp_type == "covid") {
      # --- COVID SUBGROUP LOGIC ---
      prev_all <- c(hospital = input$p_hosp/100, vent = input$p_vent/100, 
                    oxygen = input$p_oxy/100, symptoms_any = input$p_symp/100, 
                    usual_activity = input$p_act/100)
      
      # Subgroup 1: Hospital-Vent
      p_h <- prev_all['hospital']; p_v <- prev_all['vent']
      par1 <- BiCopTau2Par(family = name_to_family(input$cov_fam1), 
                           tau = tau_values[[input$cov_corr1]])
      pmf1 <- BiCop_PMF(p_h, p_v, parC = par1, specific_method = input$cov_fam1)
      names(pmf1) <- c("X00", "X01", "X10", "X11")
      cat1 <- pmf1['X11']
      
      # Subgroup 2: Hospital-Oxygen (High Corr)
      # Subgroup 2: Hospital-Oxygen
      p_o <- prev_all['oxygen']
      par2 <- BiCopTau2Par(family = name_to_family(input$cov_fam2), 
                           tau = tau_values[[input$cov_corr2]])
      pmf2 <- BiCop_PMF(p_h, p_o, parC = par2, specific_method = input$cov_fam2)
      names(pmf2) <- c("X00", "X01", "X10", "X11")
      cat2 <- pmf2['X11']
      cat3 <- pmf2['X10']
      
      # Subgroup 3: Hospital-Symptoms-Activity
      prev_g3 <- as.numeric(prev_all[c('hospital', 'symptoms_any', 'usual_activity')])
      
      # Tree 1 parameters
      tau_values[[input$cov_corr3a]]
      print(BiCopTau2Par(name_to_family(input$cov_fam3a),tau =as.numeric(tau_values[[input$cov_corr3a]])))
      
      par_t1 <- BiCopTau2Par(family = c(name_to_family(input$cov_fam3a), name_to_family(input$cov_fam3b)), 
                             tau = as.numeric(c(tau_values[[input$cov_corr3a]], tau_values[[input$cov_corr3b]])))
      print(par_t1)
      # Tree 2 parameters
      par_t2 <- BiCopTau2Par(family = name_to_family(input$cov_fam3c), 
                             tau = tau_values[[input$cov_corr3c]])
      
      s_info <- VCO_Initialize(nComp = 3, prev = prev_g3, vine_type = 'D',
                               method_list = list(c(input$cov_fam3a, input$cov_fam3b), input$cov_fam3c),
                               param_list = list(par_t1, par_t2))
      
      res_g3 <- gen_jointpmf(s_info)
      cat4 <- res_g3$jointpmf['X011',]
      cat5 <- res_g3$jointpmf['X010',]
      #cat6 <- res_g3$jointpmf['X000',]
      
      # Use residual probability so the six ordinal categories are exhaustive
      cat6 <- 1 - sum(cat1, cat2, cat3, cat4, cat5)
      
      shiny::validate(
        shiny::need(
          is.finite(cat6) && cat6 >= 0,
          "The COVID-19 category probabilities are not valid: the first five categories sum to more than 1. Please check the prevalence and dependence inputs."
        )
      )
      
      ord_pmf <- matrix(
        as.numeric(c(cat1, cat2, cat3, cat4, cat5, cat6)),
        nrow = 1
      )
      
      colnames(ord_pmf) <- paste0("Cat", 1:6)
      
    } else {
      # --- STANDARD VINE LOGIC (SUM / ORDER) ---
      # Use 20% prevalence as placeholder, or link to a separate UI table
      # p_ct <- rep(0.2, n_comp) 
      p_ct <- prev_vec()
      n_comp <- length(p_ct)
      # n_comp <- input$n_comp
      ntree <- n_comp - 1
      
      
      # --- COMMON PARAMS ---
      if (input$dep_mode == "easy") {
        decay <- if (input$add_decay) input$decay_val else 1
        corr_level <- input$corr_level
        
        tau_seq <- round(tau_values[[corr_level]] * decay^(0:(ntree - 1)), 3)
        
        param <- BiCopTau2Par(
          family = rep(name_to_family(input$cov_fam), ntree),
          tau = tau_seq
        )
        
        s_info <- VCO_Initialize(
          nComp = n_comp,
          prev = p_ct,
          vine_type = input$vine_type,
          method_list = rep(input$cov_fam, ntree),
          param_list = param
        )
      }else{ # input$dep_mode == "manual"
        tbl <- manual_tbl_clean()
        vco_inputs <- manual_table_to_vco_inputs(tbl)
        
        s_info <- VCO_Initialize(
          nComp = n_comp,
          prev = p_ct,
          vine_type = input$vine_type,
          method_list = vco_inputs$method_list,
          param_list = vco_inputs$param_list
        )
      }
      
      j_pmf <- gen_jointpmf(s_info)$jointpmf
      
      if (input$comp_type == "sum") {
        ord_pmf <- Ordinal_PMF(j_pmf, method = 'sum')
      } else {
        # Order logic
        sev <- as.numeric(trimws(unlist(strsplit(input$sev_score, ","))))
        shiny::validate(
          need(length(sev) == n_comp,
               paste0("Severity Order must have exactly ", n_comp, " values.")),
          need(all(is.finite(sev)), "Severity Order must contain numeric values."),
          need(length(unique(sev)) == n_comp, "Severity Order values must be unique.")
        )
        ord_pmf <- Ordinal_PMF(j_pmf, method = 'order', severity_score_increasing = sev)
      }
      colnames(ord_pmf) <- paste0("Cat ", 1:ncol(ord_pmf))
    }
    return(ord_pmf)
  })
  
  observe({
    res <- composite_results()
    probs <- as.numeric(res)
    
    # COVID-19 PMF is displayed from most severe to least severe.
    # Reverse it for copy/paste into IRT/CRT tabs: least severe to most severe.
    probs_copy <- if (isTRUE(input$comp_type == "covid")) {
      rev(probs)
    } else {
      probs
    }
    
    dist_text <- format_prob_vec_for_copy_5(probs_copy, digits = 5)
    
    updateTextAreaInput(
      session,
      "comp_dist_copy",
      value = dist_text
    )
  })
  
  output$comp_dist_copy_note <- renderText({
    res <- composite_results()
    probs <- as.numeric(res)
    
    probs_copy <- if (isTRUE(input$comp_type == "covid")) {
      rev(probs)
    } else {
      probs
    }
    
    dist_text <- format_prob_vec_for_copy_5(probs_copy, digits = 5)
    probs_displayed <- as.numeric(trimws(unlist(strsplit(dist_text, ","))))
    
    base_text <- paste0(
      "Number of categories: ",
      length(probs_displayed),
      " | Sum: ",
      sprintf("%.5f", sum(probs_displayed))
    )
    
    if (isTRUE(input$comp_type == "covid")) {
      paste0(
        base_text,
        " | Note: For the COVID-19 example, the displayed categories are ordered from most severe to least severe, so the copied distribution is reversed to least severe to most severe for use in the IRT/CRT tabs."
      )
    } else {
      base_text
    }
  })
  
  
  # Plotting the Distribution
  output$distPlot <- renderPlot({
    data <- composite_results()
    df <- data.frame(
      Category = colnames(data),
      Probability = as.numeric(data),
      Naive = 1/ncol(data)
    )
    
    ggplot(df, aes(x = Category)) +
      geom_col(aes(y = Probability, fill = "Vine Copula Derived"), alpha = 0.7) +
      geom_hline(aes(yintercept = Naive, linetype = "Naive (Uniform)"), color = "red") +
      scale_fill_manual(values = c("Vine Copula Derived" = "#2c3e50")) +
      labs(title = "Distribution of Ordinal Composite Outcome",
           y = "Probability", fill = "", linetype = "") +
      theme_minimal() +
      theme(
        legend.text = element_text(size = 14),
        legend.title = element_text(size = 14)
      )
  })
  
  output$pmfTable <- renderTable({
    res <- composite_results()
    probs <- as.numeric(res)
    
    if (input$comp_type == "covid") {
      # Specific clinical labels for COVID-19
      cat_names <- c(
        "1) hospitalized on mechanical ventilation or ECMO",
        "2) hospitalized on supplemental oxygen",
        "3) hospitalized not on supplemental oxygen",
        "4) not hospitalized with symptoms and limitation in activity",
        "5) not hospitalized with symptoms but with no limitation in activity",
        "6) not hospitalized without symptoms nor limitation in activity"
      )
    } else if (input$comp_type == "sum") {
      # Dynamic labels for Summation: 0, 1, ..., n-1
      # Note: length(probs) is always (n_components + 1) for the Sum mode
      cat_names <- paste("Number of events = ", 0:(ncol(res) - 1))
    } else {
      # Default for 'Rank by Severity' (Order)
      cat_names <- paste("Category", 1:ncol(res))
    }
    
    # output$bicopPlot <- renderPlot({
    #   if (input$comp_type == "covid") {
    #     print(
    #       plot_covid_bicop_shapes(
    #         fam1 = input$cov_fam1,
    #         tau1 = tau_values[[input$cov_corr1]],
    #         fam2 = input$cov_fam2,
    #         tau2 = tau_values[[input$cov_corr2]],
    #         fam3a = input$cov_fam3a,
    #         tau3a = tau_values[[input$cov_corr3a]],
    #         fam3b = input$cov_fam3b,
    #         tau3b = tau_values[[input$cov_corr3b]],
    #         fam3c = input$cov_fam3c,
    #         tau3c = tau_values[[input$cov_corr3c]]
    #       )
    #     )
    #   } else {
    #     print(plot_blank_message("Bi-copula shapes are currently shown only for COVID mode."))
    #   }
    # })
    
    
    data.frame(
      `Ordinal Category` = cat_names,
      Probability = probs,
      check.names = FALSE
    )
  }, digits = 4)
  
  
}

shinyApp(ui = ui, server = server)