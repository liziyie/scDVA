# Load packages----
getPackage <-
  function(pkg,
           check = TRUE,
           load = TRUE,
           silent = FALSE,
           github = NULL) {
    if (check) {
      if (!suppressMessages(suppressWarnings(require(
        pkg, character.only = TRUE, quietly = TRUE
      )))) {
        if (is.null(github)) {
          try(install.packages(pkg), silent = TRUE)
        }
        else{
          try(remotes::install_github(github))
        }
      }
    }
    if (load)
      suppressPackageStartupMessages(library(pkg, character.only = TRUE, quietly = TRUE))
    if (load & !silent)
      message("Loaded ", pkg)
  }
packages <-
  c("shiny",
    "shinydashboard",
    "shinyWidgets",
    "shinyjs",
    "bsplus",
    "DT",
    "dplyr")
lapply(packages, getPackage)
getPackage(pkg = 'shinydashboardPlus', github = "RinteRface/shinydashboardPlus")
getPackage(pkg = 'shinysky', github = "AnalytixWare/ShinySky")
getPackage(pkg = 'shinyauthr', github = "paulc91/shinyauthr")
getPackage(pkg = 'shinyTree', github = "shinyTree/shinyTree")

# Load functions----
source("./dataprepare_utils.R")
source("./plot_utils.R")

# Close the window and shut down app----
jscode <- "shinyjs.closeWindow = function() { window.close(); }"

# Password for the website----
user_base <- read.csv("./data/accounts.csv")

# Initialize the website with a small dataset----
load("./data/Initialize_expression.rda", envir = .GlobalEnv)
load("./data/Initialize_metadata.rda", envir = .GlobalEnv)
load("./data/Saved_genes_panel.rda", envir = .GlobalEnv)

# Shiny UI----

# Shiny header----
header <- dashboardHeaderPlus(
  title = span(
    span("s", style = "color: #F3746C; font-size: 24px; font-weight: bold"),
    span("ingle "),
    span("c", style = "color: #F3746C; font-size: 24px; font-weight: bold"),
    span("ell RNA-seq "),
    span("D", style = "color: #F3746C; font-size: 24px; font-weight: bold"),
    span("ata "),
    span("V", style = "color: #F3746C; font-size: 24px; font-weight: bold"),
    span("isualization and "),
    span("A", style = "color: #F3746C; font-size: 24px; font-weight: bold"),
    span("nalysis")
  ),
  titleWidth = 550,
  tags$li(
    class = "dropdown",
    extendShinyjs(text = jscode, functions = c("closeWindow")),
    shiny::actionButton(
      inputId = "close",
      label = "Shut down",
      width = "120px",
      icon = icon("window-close"),
      style = "color: #fff; background-color: #337ab7; border-color: #2e6da4"
    )
  )
)

# Shiny sidebar----
sidebar <- dashboardSidebar(
  tags$head(
    tags$style(
      type = 'text/css',
      ".nav-tabs {font-size: 14px} ",
      ".main-sidebar {font-size: 16px}"
    )
  ),
  sidebarMenu(
    menuItem("Embedding", tabName = "Embedding", icon = icon("spinner")),
    menuItem(
      "Distribution",
      tabName = "Distribution",
      icon = icon("bar-chart-o")
    ),
    menuItem(
      "Significance",
      tabName = "Significance",
      icon = icon("braille")
    ),
    menuItem("Heatmap", tabName = "Heatmap", icon = icon("delicious")),
    menuItem(
      "In-silico FACS",
      tabName = "FACS",
      icon = icon("fas fa-fingerprint")
    ),
    menuItem(
      "Metadata",
      tabName = "Metadata",
      icon = icon("stats", lib = "glyphicon")
    ),
    menuItem("DataTable", tabName = "DataTable", icon = icon("table")),
    menuItem("Instructions", tabName = "Instructions", icon = icon("file-alt")),
    menuItem("About", tabName = "About", icon = icon("envelope"))
  )
)

# Shiny body----
body <- dashboardBody(
  useShinyjs(),
  # *Login page----
  shinyauthr::loginUI("login"),
  # *Working page----
  conditionalPanel(condition = "output.Login_logged",
                   fluidRow(
                     column(
                       width = 9,
                       # **Head panel----
                       box(
                         title = "Gene input",
                         width = 4,
                         status = "warning",
                         collapsible = T,
                         collapsed = F,
                         tabBox(
                           id = "GeneInput_panel",
                           side = "left",
                           selected = "Genes",
                           width = 12,
                           tabPanel(
                             "Genes",
                             textAreaInput(
                               inputId = "GeneInput_text",
                               label = "Type a gene or geneset:",
                               value = "CD14"
                             )
                           ),
                           tabPanel(
                             "Saved",
                             textAreaInput(
                               inputId = "GeneInput_saved_text",
                               label = NULL,
                               value = NULL
                             ),
                             selectizeInput(
                               inputId = "GeneInput_saved",
                               label = "Select from the saved geneset:",
                               multiple = TRUE,
                               choice = Saved_genes_panel,
                               selected = NULL
                             )
                           ),
                           tabPanel(
                             "Upload",
                             textAreaInput(
                               inputId = "GeneInput_file_text",
                               label = NULL,
                               value = NULL
                             ),
                             fileInput(
                               inputId = "GeneInput_file",
                               label = "Choose a csv file:",
                               multiple = FALSE
                             ) %>%
                               shinyInput_label_embed(
                                 shiny_iconlink() %>%
                                   bs_embed_popover(
                                     title = "Note:",
                                     placement = "left",
                                     content = "There should be a column named as 'Symbol'."
                                   )
                               )
                           )
                         ),
                         shiny::actionButton(
                           inputId = "GeneInput_submit",
                           label = "Submit",
                           width = "120px",
                           icon = icon("paper-plane"),
                           style = "color: #fff; background-color: #337ab7; border-color: #2e6da4"
                         )
                       ),
                       box(
                         title = "Plot size",
                         width = 4,
                         status = "warning",
                         collapsible = T,
                         collapsed = F,
                         numericInput(
                           inputId = "PlotSize_width",
                           label = "Plot width (px)",
                           value = 960,
                           step = 10
                         ),
                         numericInput(
                           inputId = "PlotSize_height",
                           label = "Plot height (px)",
                           value = 960,
                           step = 10
                         ),
                         p(
                           class = "text-muted",
                           paste(
                             "Note: Please click the submit button in 'Gene Input' box after you change the figure size."
                           )
                         )
                       ),
                       box(
                         title = "Plot parameters",
                         width = 4,
                         status = "warning",
                         collapsible = T,
                         collapsed = F,
                         selectInput(
                           inputId = "PlotPar_facet",
                           label = "Multi gene",
                           choices = c("Seperate" = TRUE, "Geometric mean" = FALSE),
                           selected = "Seperate"
                         ),
                         numericInput(
                           inputId = "PlotPar_nrow",
                           label = "Row number",
                           value = 2
                         ),
                         numericInput(
                           inputId = "PlotPar_fontsize",
                           label = "Font size",
                           value = 16,
                           min = 9,
                           max = 26,
                           step = 1
                         )
                       ),
                       # **Plot region----
                       box(
                         status = "warning",
                         width = 12,
                         tabItems(
                           # ****Embedding Plot----
                           tabItem(
                             tabName = "Embedding",
                             box(
                               status = NULL,
                               width = 12,
                               dropdownButton(
                                 selectInput(
                                   inputId = "Embedding_tag1",
                                   label = "Embedding tag",
                                   choices = c("Global", "Sub"),
                                   selected = "Sub"
                                 ),
                                 selectInput(
                                   inputId = "Embedding_used1",
                                   label = "Embedding used",
                                   choices = c("tSNE", "UMAP"),
                                   selected = "tSNE"
                                 ),
                                 sliderInput(
                                   inputId = "Embedding_dotsize1",
                                   label = "Dot size",
                                   value = 1,
                                   min = 0,
                                   max = 3,
                                   step = 0.2
                                 ),
                                 selectInput(
                                   inputId = "Embedding_colorpanel",
                                   label = "Color profile",
                                   choices = c(
                                     "YlOrRd",
                                     "BrBG",
                                     "RdYlBu",
                                     "RdYlGn",
                                     "Spectral",
                                     "YlOrBr",
                                     "YlGnBu",
                                     "Reds",
                                     "Blues",
                                     "OrRd"
                                   ),
                                   selected = "YlOrRd"
                                 ),
                                 circle = TRUE,
                                 status = "danger",
                                 icon = icon("spinner"),
                                 width = "250px",
                                 tooltip = tooltipOptions(title = "Click to see inputs!")
                               ),
                               uiOutput("Embedding_geneplot.ui")
                             ),
                             box(
                               status = NULL,
                               width = 12,
                               dropdownButton(
                                 selectInput(
                                   inputId = "Embedding_colorby",
                                   label = "Color by",
                                   choices = c(
                                     "Global_Cluster",
                                     "Sub_Cluster",
                                     "Tissue",
                                     "Treatment",
                                     "Day",
                                     "Sample",
                                     "nGene",
                                     "nUMI"
                                   ),
                                   selected = "Sub_Cluster"
                                 ),
                                 selectInput(
                                   inputId = "Embedding_tag2",
                                   label = "Embedding tag",
                                   choices = c("Global", "Sub"),
                                   selected = "Sub"
                                 ),
                                 selectInput(
                                   inputId = "Embedding_used2",
                                   label = "Embedding used",
                                   choices = c("tSNE", "UMAP"),
                                   selected = "tSNE"
                                 ),
                                 sliderInput(
                                   inputId = "Embedding_dotsize2",
                                   label = "Dot size",
                                   value = 1,
                                   min = 0,
                                   max = 3,
                                   step = 0.2
                                 ),
                                 circle = TRUE,
                                 status = "danger",
                                 icon = icon("spinner"),
                                 width = "250px",
                                 tooltip = tooltipOptions(title = "Click to see inputs!")
                               ),
                               uiOutput("Embedding_metaplot.ui")
                             ),
                             downloadButton(
                               outputId = 'Embedding_download',
                               label = 'Download',
                               width = "120px",
                               style = "color: #fff; background-color: #337ab7; border-color: #2e6da4"
                             )
                           ),
                           # ****Distribution Plot----
                           tabItem(
                             tabName = "Distribution",
                             dropdownButton(
                               selectInput(
                                 inputId = "Distribution_plottype",
                                 label = "Plot type",
                                 choices = c("Box Plot", "Violin Plot"),
                                 selected = "Box Plot"
                               ),
                               selectInput(
                                 inputId = "Distribution_groupby",
                                 label = "Group by",
                                 choices = c(
                                   "Global_Cluster",
                                   "Sub_Cluster",
                                   "Tissue",
                                   "Treatment",
                                   "Day",
                                   "Sample"
                                 ),
                                 selected = "Sub_Cluster"
                               ),
                               selectInput(
                                 inputId = "Distribution_colorby",
                                 label = "Color by",
                                 choices = c(
                                   "Exp",
                                   "Global_Cluster",
                                   "Sub_Cluster",
                                   "Tissue",
                                   "Treatment",
                                   "Day",
                                   "Sample"
                                 ),
                                 selected = "Exp"
                               ),
                               selectInput(
                                 inputId = "Distribution_colorpanel",
                                 label = "Color profile",
                                 choices = c(
                                   "YlOrRd",
                                   "BrBG",
                                   "RdYlBu",
                                   "RdYlGn",
                                   "Spectral",
                                   "YlOrBr",
                                   "YlGnBu",
                                   "Reds",
                                   "Blues",
                                   "OrRd",
                                   "RdBu"
                                 ),
                                 selected = "YlOrRd"
                               ),
                               sliderInput(
                                 inputId = "Distribution_dotsize",
                                 label = "Dot size",
                                 value = 0.6,
                                 min = 0,
                                 max = 1,
                                 step = 0.1
                               ),
                               circle = TRUE,
                               status = "danger",
                               icon = icon("bar-chart-o"),
                               width = "250px",
                               tooltip = tooltipOptions(title = "Click to see inputs!")
                             ),
                             uiOutput("Distribution_plot.ui"),
                             downloadButton(
                               outputId = 'Distribution_download',
                               label = 'Download',
                               width = "120px",
                               style = "color: #fff; background-color: #337ab7; border-color: #2e6da4"
                             )
                           ),
                           # ****Significance Plot----
                           tabItem(
                             tabName = "Significance",
                             dropdownButton(
                               selectInput(
                                 inputId = "Significance_groupby",
                                 label = "Group by",
                                 choices = c(
                                   "Global_Cluster",
                                   "Sub_Cluster",
                                   "Tissue",
                                   "Treatment",
                                   "Day",
                                   "Sample"
                                 ),
                                 selected = "Sub_Cluster"
                               ),
                               numericInput(
                                 inputId = "Significance_percutoff",
                                 label = "Expression cutoff",
                                 value = 0
                               ),
                               numericInput(
                                 inputId = "Significance_siglevel",
                                 label = "Significant level",
                                 value = 0.05
                               ),
                               circle = TRUE,
                               status = "danger",
                               icon = icon("braille"),
                               width = "250px",
                               tooltip = tooltipOptions(title = "Click to see inputs!")
                             ),
                             DT::dataTableOutput("Significance_group"),
                             tags$br(),
                             h4("Significance"),
                             verbatimTextOutput("Significance_statistic"),
                             tags$br(),
                             uiOutput("Significance_plot.ui"),
                             downloadButton(
                               outputId = 'Significance_download',
                               label = 'Download',
                               width = "120px",
                               style = "color: #fff; background-color: #337ab7; border-color: #2e6da4"
                             )
                           ),
                           # ****Heatmap Plot----
                           tabItem(
                             tabName = "Heatmap",
                             dropdownButton(
                               selectInput(
                                 inputId = "Heatmap_groupby",
                                 label = "Group by",
                                 choices = c(
                                   "Global_Cluster",
                                   "Sub_Cluster",
                                   "Tissue",
                                   "Treatment",
                                   "Day",
                                   "Sample"
                                 ),
                                 selected = "Sub_Cluster"
                               ),
                               selectInput(
                                 inputId = "Heatmap_colorpanel",
                                 label = "Color profile",
                                 choices = c(
                                   "YlOrRd",
                                   "BrBG",
                                   "RdYlBu",
                                   "RdYlGn",
                                   "Spectral",
                                   "YlOrBr",
                                   "YlGnBu",
                                   "Reds",
                                   "Blues",
                                   "OrRd",
                                   "RdBu"
                                 ),
                                 selected = "RdBu"
                               ),
                               circle = TRUE,
                               status = "danger",
                               icon = icon("delicious"),
                               width = "250px",
                               tooltip = tooltipOptions(title = "Click to see inputs!")
                             ),
                             uiOutput("Heatmap_plot.ui"),
                             downloadButton(
                               outputId = 'Heatmap_download',
                               label = 'Download',
                               width = "120px",
                               style = "color: #fff; background-color: #337ab7; border-color: #2e6da4"
                             )
                           ),
                           # ****In-silico FACS Plot----
                           tabItem(
                             tabName = "FACS",
                             tabBox(
                               id = "FACS",
                               side = "left",
                               selected = "FACS",
                               width = 12,
                               tabPanel("FACS",
                                 dropdownButton(
                                   numericInput(
                                     inputId = "FACS_xcutoff",
                                     label = "x cutoff",
                                     value = 1
                                   ),
                                   numericInput(
                                     inputId = "FACS_ycutoff",
                                     label = "y cutoff",
                                     value = 1
                                   ),
                                   selectInput(
                                     inputId = "FACS_colorpanel",
                                     label = "Color profile",
                                     choices = c(
                                       "YlOrRd",
                                       "BrBG",
                                       "RdYlBu",
                                       "RdYlGn",
                                       "Spectral",
                                       "YlOrBr",
                                       "YlGnBu",
                                       "Reds",
                                       "Blues",
                                       "OrRd",
                                       "RdBu"
                                     ),
                                     selected = "RdBu"
                                   ),
                                   selectInput(
                                     inputId = "FACS_colorby",
                                     label = "Color by",
                                     choices = c(
                                       "None",
                                       "Global_Cluster",
                                       "Sub_Cluster",
                                       "Tissue",
                                       "Treatment",
                                       "Day",
                                       "Sample"
                                     ),
                                     selected = "None"
                                   ),
                                   sliderTextInput(
                                     inputId = "FACS_dotsize",
                                     label = "Dot size",
                                     choices = c(-1, seq(0, 3, 0.5)),
                                     grid = TRUE,
                                     selected = -1
                                   ),
                                   circle = TRUE,
                                   status = "danger",
                                   icon = icon("fas fa-fingerprint"),
                                   width = "250px",
                                   tooltip = tooltipOptions(title = "Click to see inputs!")
                                 ),
                                 uiOutput("FACS_plot.ui"),
                                 downloadButton(
                                   outputId = 'FACS_download',
                                   label = 'Download',
                                   width = "120px",
                                   style = "color: #fff; background-color: #337ab7; border-color: #2e6da4")
                               ),
                               tabPanel("DEG",
                                 box(width = 12,
                                   column(width = 4,
                                     selectInput(
                                       inputId = "FACS_DEG_group1",
                                       label = "The 1st group",
                                       choices = c("Group1", "Group2", "Group3", "Group4"),
                                       selected = "Group1"
                                     ),
                                     selectInput(
                                       inputId = "FACS_DEG_group2",
                                       label = "The 2nd group",
                                       choices = c("Group1", "Group2", "Group3", "Group4"),
                                       selected = "Group2"
                                     ),
                                     p(
                                       class = "text-muted",
                                       paste(
                                         "Note: The number of cells in each group will be downsampled to 1000 randomly."
                                       )
                                     )
                                   ),
                                   column(width = 4,
                                     numericInput(
                                       inputId = "FACS_Volcano_Pvalue_cutoff",
                                       label = "adj.P.Val cutoff",
                                       value = 0.05
                                     ),
                                     numericInput(
                                       inputId = "FACS_Volcano_logFC_cutoff",
                                       label = "logFC cutoff",
                                       value = 0.5
                                     ),
                                     numericInput(
                                       inputId = "FACS_Volcano_ngenes_labeled",
                                       label = "Labeled genes",
                                       value = 25
                                     )
                                   ),
                                   column(width = 4,
                                     numericInput(
                                       inputId = "FACS_Volcano_fontsize",
                                       label = "Font size",
                                       value = 15
                                     ),
                                     sliderTextInput(
                                       inputId = "FACS_Volcano_dotsize",
                                       label = "Dot size",
                                       choices = seq(0, 10, 0.5),
                                       grid = TRUE,
                                       selected = 2
                                     ),
                                     shiny::actionButton(
                                       inputId = "FACS_DEG_submit",
                                       label = "Calculate",
                                       width = "120px",
                                       icon = icon("calculator"),
                                       style = "color: #fff; background-color: #337ab7; border-color: #2e6da4"
                                     )
                                   )
                                 ),
                                 box(width = 12,
                                   uiOutput("FACS_volcano_plot.ui"),
                                   downloadButton(
                                     outputId = 'FACS_volcano_download',
                                     label = 'Download',
                                     width = "120px",
                                     style = "color: #fff; background-color: #337ab7; border-color: #2e6da4"
                                   )
                                 ),
                                 box(width = 12,
                                   h3("Differentially expressed genes"),
                                   p(
                                     class = "text-muted",
                                     paste(
                                       "Note: Here only shows 2000 genes with the smallest adj.P.Val. If you want the full gene list, you can download it."
                                     )
                                   ),
                                   DT::dataTableOutput("FACS_DEGDataTable_output"),
                                   downloadButton(
                                     'FACS_DEGDataTable_download',
                                     label = 'Download',
                                     width = "120px",
                                     style = "color: #fff; background-color: #337ab7; border-color: #2e6da4"
                                   )
                                 )
                               )
                             )
                           ),
                           # ****Metadata Plot----
                           tabItem(
                             tabName = "Metadata",
                             dropdownButton(
                               selectInput(
                                 inputId = "Metadata_plottype",
                                 label = "Plot type",
                                 choices = c("Bar Plot", "Pie Plot"),
                                 selected = "Bar Plot"
                               ),
                               selectInput(
                                 inputId = "Metadata_groupby1",
                                 label = "Group by1",
                                 choices = c(
                                   "Global_Cluster",
                                   "Sub_Cluster",
                                   "Tissue",
                                   "Treatment",
                                   "Day",
                                   "Sample"
                                 ),
                                 selected = "Sub_Cluster"
                               ),
                               selectInput(
                                 inputId = "Metadata_colorby",
                                 label = "Color by",
                                 choices = c(
                                   "Global_Cluster",
                                   "Sub_Cluster",
                                   "Tissue",
                                   "Treatment",
                                   "Day",
                                   "Sample"
                                 ),
                                 selected = "Sample"
                               ),
                               selectInput(
                                 inputId = "Metadata_groupby2",
                                 label = "Group by2",
                                 choices = c(
                                   "None",
                                   "Global_Cluster",
                                   "Sub_Cluster",
                                   "Tissue",
                                   "Treatment",
                                   "Day",
                                   "Sample"
                                 ),
                                 selected = "None"
                               ),
                               selectInput(
                                 inputId = "Metadata_quantify",
                                 label = "Quantified by",
                                 choices = c("Proportion", "Count"),
                                 selected = "Proportion"
                               ),
                               materialSwitch(
                                 inputId = "Metadata_flip",
                                 label = "Coordinates flipped",
                                 value = TRUE,
                                 status = "primary",
                                 right  = TRUE
                               ),
                               numericInput(
                                 inputId = "Metadata_cutoff",
                                 label = "Proportion shown in pie plot >",
                                 value = 10
                               ),
                               circle = TRUE,
                               status = "danger",
                               icon = icon("stats", lib = "glyphicon"),
                               width = "250px",
                               tooltip = tooltipOptions(title = "Click to see inputs!")
                             ),
                             uiOutput("Metadata_plot.ui"),
                             downloadButton(
                               outputId = 'Metadata_download',
                               label = 'Download',
                               width = "120px",
                               style = "color: #fff; background-color: #337ab7; border-color: #2e6da4"
                             )
                           ),
                           # ****DataTable Output----
                           tabItem(
                             tabName = "DataTable",
                             DT::dataTableOutput("DataTable_output"),
                             downloadButton(
                               'DataTable_download',
                               label = 'Download',
                               width = "120px",
                               style = "color: #fff; background-color: #337ab7; border-color: #2e6da4"
                             )
                           ),
                           # ****Instructions Output----
                           tabItem(
                             tabName = "Instructions",
                             tags$iframe(style="height:1200px; width:100%", 
                                         src="User Manual of scDVA.pdf")
                           ),
                           # ****About Output----
                           tabItem(
                             tabName = "About",
                             h2("scDVA"),
                             HTML(
                               "<p style = 'margin-left:3%'><font size='4'>We developed this interactive web server scDVA (short for single cell RNA-seq data visualization and analyzation) for users to explore and analyze the single cell RNA-seq data produced in our lab.</font></p>"
                             ),
                             HTML(
                               "<p style = 'margin-left:3%'><font size='4'>This tool is developed by Ziyi Li, Tao Wang of <a href = 'http://cancer-pku.cn/'>Zemin Zhang's lab</a>, Peking University.</font></p>"
                             ),
                             h2("GEPIA2"),
                             HTML(
                               "<p style = 'margin-left:3%'><font size='4'>Please try <a href = 'http://gepia2.cancer-pku.cn/#index')>GEPIA2</a> to analyze the bulk RNA-seq data from the TCGA and GTEx projects.</font></p>"
                             ),
                             h2("Contact us"),
                             HTML(
                               "<p style = 'margin-left:3%'><font size='4'>If any question about the use of scDVA or find any bug, please contact us: <a href = 'mailto:liziyie@pku.edu.cn'>liziyie@pku.edu.cn</a></font></p>"
                             )
                           )
                         )
                       ),
                       busyIndicator(text = "Please wait...", wait = 1000)
                     ),
                     # **Right dataset selected panel----
                     column(
                       width = 3,
                       box(
                         title = "Load data",
                         width = NULL,
                         status = "danger",
                         collapsible = T,
                         collapsed = F,
                         shinyTree(
                           "LoadDataTree", checkbox = TRUE, search = TRUE,
                           searchtime = 1000, theme = "proton"
                         ),
                         p(
                           class = "text-muted",
                           paste("Note: please select one or more dataset(s) to load.")
                         ),
                         shiny::actionButton(
                           inputId = "LoadButton",
                           label = "Load",
                           width = 120,
                           icon = icon("upload"),
                           style = "color: #fff; background-color: #337ab7; border-color: #2e6da4"
                         ),
                         shinyalert("LoadDataset_alert", FALSE, auto.close.after = 5)
                       ),
                       box(
                         title = "Select a dataset",
                         width = NULL,
                         status = "warning",
                         collapsible = T,
                         collapsed = F,
                         selectInput(
                           inputId = "SelectData_dataset",
                           label = "Dataset",
                           choices = "Initialize"
                         ),
                         selectInput(
                           inputId = "SelectData_normalization",
                           label = "Normalization",
                           choices = c("tpm", "counts")
                         )
                       ),
                       box(
                         title = "Subset the dataset",
                         width = NULL,
                         status = "warning",
                         pickerInput(
                           inputId = "SubsetData_globalcluster",
                           label = "Global cluster",
                           choices = NULL,
                           multiple = TRUE,
                           options = list(
                             `actions-box` = TRUE,
                             `live-search` = TRUE,
                             size = 8
                           )
                         ),
                         pickerInput(
                           inputId = "SubsetData_subcluster",
                           label = "Sub cluster",
                           choices = NULL,
                           multiple = TRUE,
                           options = list(
                             `actions-box` = TRUE,
                             `live-search` = TRUE,
                             size = 8
                           )
                         ),
                         pickerInput(
                           inputId = "SubsetData_tissue",
                           label = "Tissue",
                           choices = NULL,
                           multiple = TRUE,
                           options = list(
                             `actions-box` = TRUE,
                             `live-search` = TRUE,
                             size = 8
                           )
                         ),
                         pickerInput(
                           inputId = "SubsetData_treatment",
                           label = "Treatment",
                           choices = NULL,
                           multiple = TRUE,
                           options = list(
                             `actions-box` = TRUE,
                             `live-search` = TRUE,
                             size = 8
                           )
                         ),
                         pickerInput(
                           inputId = "SubsetData_day",
                           label = "Day",
                           choices = NULL,
                           multiple = TRUE,
                           options = list(
                             `actions-box` = TRUE,
                             `live-search` = TRUE,
                             size = 8
                           )
                         ),
                         pickerInput(
                           inputId = "SubsetData_sample",
                           label = "Sample",
                           choices = NULL,
                           multiple = TRUE,
                           options = list(
                             `actions-box` = TRUE,
                             `live-search` = TRUE,
                             size = 8
                           )
                         ),
                         p(
                           class = "text-muted",
                           paste(
                             "Note: Please click the submit button in 'Gene Input' box after you subset the dataset."
                           )
                         )
                       )
                     )
                   ))
)


# Define UI----
ui <- dashboardPagePlus(header, sidebar, body, use_bs_popover())

# Shiny Server----
server <- function(input, output, session) {
  # Close the server----
  observeEvent(input$close, {
    js$closeWindow()
    stopApp()
  })
  
  # # Loggin Page----
  credentials <- callModule(
    shinyauthr::login,
    "login",
    data = user_base,
    user_col = user,
    pwd_col = password,
    log_out = reactive(logout_init())
  )
  
  logout_init <-
    callModule(shinyauthr::logout,
               "logout",
               reactive(credentials()$user_auth))
  
  output$Login_logged <- reactive({
    req(credentials()$user_auth)
  })
  outputOptions(output, 'Login_logged', suspendWhenHidden = FALSE)
  
  # Load data tree-----
  output$LoadDataTree <- renderTree({
    dataset_tree
  })
  
  # Interactive load button----
  observeEvent(input$LoadButton, {
    req(input$LoadDataTree)
    LoadData_selected <- get_selected(input$LoadDataTree, format = "classid") # Record the dataset to load
    LoadData_selected <- LoadData_selected[LoadData_selected %in% dataset_map$DatasetName]
    if (is.null(LoadData_selected)) {
      LoadData_selected <- character(0)
    }
    LoadData.status <-
      loadDataSet(LoadData_selected) # Load the dataset
    showshinyalert(session, "LoadDataset_alert", "Done.",
                   styleclass = "info")
    updateSelectInput(
      session,
      inputId = "SelectData_dataset",
      choices = LoadData_selected,
      selected = LoadData_selected[1]
    )
  })
  
  # Update the panel----
  # Update the SubsetData panel
  observe({
    SubsetData_panel_list <- getPanelList(input$SelectData_dataset)
    updatePickerInput(
      session = session,
      inputId = "SubsetData_globalcluster",
      choices = SubsetData_panel_list$GlobalCluster,
      selected = SubsetData_panel_list$GlobalCluster[1]
    )
    updatePickerInput(
      session = session,
      inputId = "SubsetData_subcluster",
      choices = SubsetData_panel_list$SubCluster_list,
      selected = unlist(SubsetData_panel_list$SubCluster_list)
    )
    updatePickerInput(
      session = session,
      inputId = "SubsetData_tissue",
      choices = SubsetData_panel_list$Tissue,
      selected = SubsetData_panel_list$Tissue
    )
    updatePickerInput(
      session = session,
      inputId = "SubsetData_sample",
      choices = SubsetData_panel_list$Sample,
      selected = SubsetData_panel_list$Sample
    )
    updatePickerInput(
      session = session,
      inputId = "SubsetData_treatment",
      choices = SubsetData_panel_list$Treatment,
      selected = SubsetData_panel_list$Treatment
    )
    updatePickerInput(
      session = session,
      inputId = "SubsetData_day",
      choices = SubsetData_panel_list$Day,
      selected = SubsetData_panel_list$Day
    )
  })
  # Update the GeneInput panel
  observe({
    updateTextAreaInput(
      session = session,
      inputId = "GeneInput_saved_text",
      value = input$GeneInput_saved
    )
    if (length(input$GeneInput_file$datapath) != 0) {
      df <-
        read.csv(
          input$GeneInput_file$datapath,
          header = T,
          stringsAsFactors = F
        )
      genes_upload <- paste(as.character(df$Symbol), collapse = ",")
      updateTextAreaInput(session = session,
                          inputId = "GeneInput_file_text",
                          value = genes_upload)
    }
  })
  
  # Initialize InputValue----
  InputValue <- reactiveValues(
    geneinput = "CD14",
    geneinfo = checkGeneList(isolate("CD14"), "Initialize", "tpm"),
    normalization = "tpm",
    plot_data = getPlotData(
      Genes = "CD14",
      SelectData_dataset = "Initialize",
      SelectData_normalization = "tpm",
      GlobalCluster.selected = "Myeloid cell",
      SubCluster.selected = getPanelList("Initialize")$SubCluster_list[["Myeloid cell"]],
      Sample.selected = getPanelList("Initialize")$Sample,
      Tissue.selected = getPanelList("Initialize")$Tissue,
      Treatment.selected = getPanelList("Initialize")$Treatment,
      Day.selected = getPanelList("Initialize")$Day
    ),
    plot_width = 960,
    plot_height = 960
  )
  
  # Initialize in-silico FACS DEGenes table----
  FACS_DEGenes <- reactiveValues(
    limma = data.frame(logFC = 0, Warning = "Please click the 'Calculate' button!")
  )
  
  # Update InputValue----
  observeEvent(input$GeneInput_submit, {
    if (input$GeneInput_panel == "Genes") {
      InputValue$geneinput <- toupper(input$GeneInput_text)
    } else if (input$GeneInput_panel == "Saved") {
      InputValue$geneinput <- toupper(input$GeneInput_saved_text)
    } else if (input$GeneInput_panel == "Upload") {
      InputValue$geneinput <- toupper(input$GeneInput_file_text)
    }
    InputValue$dataset <- input$SelectData_dataset
    InputValue$normalization <- input$SelectData_normalization
    InputValue$geneinfo <-
      checkGeneList(isolate(InputValue$geneinput),
                    InputValue$dataset,
                    InputValue$normalization)
    if (InputValue$geneinfo$gene_number >= 1) {
      if (!InputValue$geneinfo$all_genes_avaliable_flag) {
        message = paste(
          paste(InputValue$geneinfo$wrong_gene, collapse = ","),
          "not found in the dataset!"
        )
        sendSweetAlert(session = session,
                       title = message,
                       type = "warning")
      }
      InputValue$plot_data <-
        getPlotData(
          Genes = InputValue$geneinfo$right_gene,
          SelectData_dataset = InputValue$dataset,
          SelectData_normalization = InputValue$normalization,
          GlobalCluster.selected = input$SubsetData_globalcluster,
          SubCluster.selected = input$SubsetData_subcluster,
          Sample.selected = input$SubsetData_sample,
          Tissue.selected = input$SubsetData_tissue,
          Treatment.selected = input$SubsetData_treatment,
          Day.selected = input$SubsetData_day
        )
    } else{
      InputValue$geneinfo <-
        checkGeneList(isolate("CD14"), "Initialize", "tpm")
      message = ("No input genes could be used!")
      sendSweetAlert(session = session,
                     title = message,
                     type = "error")
    }
    InputValue$plot_width <- input$PlotSize_width
    InputValue$plot_height <- input$PlotSize_height
  })
  
  # Update in-silico FACS DEG calculation----
  observeEvent(input$FACS_DEG_submit, {
    FACS_DEGenes$limma <- getLIMMAData(
        Plot.data = InputValue$plot_data, 
        SelectData_dataset = InputValue$dataset, 
        SelectData_normalization = InputValue$normalization, 
        genes = InputValue$geneinfo$right_gene, 
        x.cutoff = input$FACS_xcutoff, 
        y.cutoff = input$FACS_xcutoff, 
        group1 = input$FACS_DEG_group1, 
        group2 = input$FACS_DEG_group2
      )
  })
  
  # Embedding Plot----
  output$Embedding_geneplot <- renderPlot({
    DoEmbeddingPlot(
      plot.data = InputValue$plot_data,
      genes = InputValue$geneinfo$right_gene,
      color.by = "Genes",
      color.panel = input$Embedding_colorpanel,
      embedding.tag = input$Embedding_tag1,
      embedding.used = input$Embedding_used1,
      seperate.plot = input$PlotPar_facet,
      nrow = input$PlotPar_nrow,
      dot.size = input$Embedding_dotsize1,
      font.size = input$PlotPar_fontsize
    )
  })
  
  output$Embedding_geneplot.ui <- renderUI({
    plotOutput(
      "Embedding_geneplot",
      width = InputValue$plot_width,
      height = InputValue$plot_height
    )
  })
  
  output$Embedding_metaplot <- renderPlot({
    DoEmbeddingPlot(
      plot.data = InputValue$plot_data,
      color.by = input$Embedding_colorby,
      embedding.tag = input$Embedding_tag2,
      embedding.used = input$Embedding_used2,
      seperate.plot = input$PlotPar_facet,
      nrow = input$PlotPar_nrow,
      dot.size = input$Embedding_dotsize2,
      font.size = input$PlotPar_fontsize
    )
  })
  
  output$Embedding_metaplot.ui <- renderUI({
    plotOutput(
      "Embedding_metaplot",
      width = InputValue$plot_width,
      height = InputValue$plot_height
    )
  })
  
  output$Embedding_download <- downloadHandler(
    filename = function() {
      paste0(
        input$SelectData_dataset,
        "_",
        input$SelectData_normalization,
        "_EmbeddingPlot_",
        stringi::stri_rand_strings(1, 10),
        ".pdf"
      )
    },
    content = function(file) {
      pdf(file,
          width = InputValue$plot_width / 80,
          height = InputValue$plot_height / 80)
      plot(
        DoEmbeddingPlot(
          plot.data = InputValue$plot_data,
          genes = InputValue$geneinfo$right_gene,
          color.by = "Genes",
          color.panel = input$Embedding_colorpanel,
          embedding.tag = input$Embedding_tag1,
          embedding.used = input$Embedding_used1,
          seperate.plot = input$PlotPar_facet,
          nrow = input$PlotPar_nrow,
          dot.size = input$Embedding_dotsize1,
          font.size = input$PlotPar_fontsize
        )
      )
      plot(
        DoEmbeddingPlot(
          plot.data = InputValue$plot_data,
          color.by = input$Embedding_colorby,
          embedding.tag = input$Embedding_tag2,
          embedding.used = input$Embedding_used2,
          seperate.plot = input$PlotPar_facet,
          nrow = input$PlotPar_nrow,
          dot.size = input$Embedding_dotsize2,
          font.size = input$PlotPar_fontsize
        )
      )
      dev.off()
    }
  )
  
  # Distribution Plot----
  output$Distribution_plot <- renderPlot({
    DoDistributionPlot(
      plot.data = InputValue$plot_data,
      genes = InputValue$geneinfo$right_gene,
      plot.type = input$Distribution_plottype,
      group.by = input$Distribution_groupby,
      color.by = input$Distribution_colorby,
      color.panel = input$Distribution_colorpanel,
      dot.size = input$Distribution_dotsize,
      seperate.plot = input$PlotPar_facet,
      nrow = input$PlotPar_nrow,
      font.size = input$PlotPar_fontsize
    )
  })
  
  output$Distribution_plot.ui <- renderUI({
    plotOutput(
      "Distribution_plot",
      width = 1.2 * InputValue$plot_width,
      height = 0.5 * InputValue$plot_height
    )
  })
  
  output$Distribution_download <- downloadHandler(
    filename = function() {
      paste0(
        input$SelectData_dataset,
        "_",
        input$SelectData_normalization,
        "_DistributionPlot_",
        stringi::stri_rand_strings(1, 10),
        ".pdf"
      )
    },
    content = function(file) {
      pdf(
        file,
        width = 1.2 * InputValue$plot_width / 80,
        height = 0.5 * InputValue$plot_height / 160
      )
      plot(
        DoDistributionPlot(
          plot.data = InputValue$plot_data,
          genes = InputValue$geneinfo$right_gene,
          plot.type = input$Distribution_plottype,
          group.by = input$Distribution_groupby,
          color.by = input$Distribution_colorby,
          dot.size = input$Distribution_dotsize,
          color.panel = input$Distribution_colorpanel,
          seperate.plot = input$PlotPar_facet,
          nrow = input$PlotPar_nrow,
          font.size = input$PlotPar_fontsize
        )
      )
      dev.off()
    }
  )
  
  # Significance plot----
  output$Significance_group <- DT::renderDataTable({
    CorrList <-
      getSigData(
        InputValue$plot_data,
        input$Significance_groupby,
        input$Significance_percutoff
      )
    if (!is.null(CorrList)) {
      DT::datatable(
        t(CorrList$group_def[, c("Group.cha", "Exp.percent", "Exp.mean", "Exp.sd")]),
        rownames = c(
          "Cluster",
          "Cell Percentage",
          "Expression Mean",
          "Expression Sd"
        ),
        class = "stripe table-condensed",
        options = list(
          scrollX = TRUE,
          processing = FALSE,
          pageLength = 10
        ),
        autoHideNavigation = TRUE
      )
    } else{
      DT::datatable(data.frame(Warnings = "There needs at least two groups!"))
    }
  })
  
  output$Significance_statistic <- renderText({
    CorrList <-
      getSigData(
        InputValue$plot_data,
        input$Significance_groupby,
        input$Significance_percutoff
      )
    sprintf("F-value is %s\np-value is %s",
            CorrList$Fvalue[1],
            CorrList$Fvalue[2])
  })
  
  output$Significance_plot <- renderPlot({
    DoSignificancePlot(
      plot.data = InputValue$plot_data,
      group.by = input$Significance_groupby,
      per.cutoff = input$Significance_percutoff,
      font.size = input$PlotPar_fontsize / 10,
      sig.level = input$Significance_siglevel
    )
  })
  
  output$Significance_plot.ui <- renderUI({
    plotOutput(
      "Significance_plot",
      width = 1.2 * InputValue$plot_width,
      height = InputValue$plot_width
    )
  })
  
  output$Significance_download <- downloadHandler(
    filename = function() {
      paste0(
        input$SelectData_dataset,
        "_",
        input$SelectData_normalization,
        "_SignificancePlot_",
        stringi::stri_rand_strings(1, 10),
        ".pdf"
      )
    },
    content = function(file) {
      pdf(
        file,
        width = 1.2 * InputValue$plot_width / 60,
        height = InputValue$plot_height / 60
      )
      DoSignificancePlot(
        plot.data = InputValue$plot_data,
        group.by = input$Significance_groupby,
        per.cutoff = input$Significance_percutoff,
        font.size = input$PlotPar_fontsize / 10,
        sig.level = input$Significance_siglevel
      )
      dev.off()
    }
  )
  
  # Heatmap plot----
  output$Heatmap_plot <- renderPlot({
    DoHeatmapPlot(
      plot.data = InputValue$plot_data,
      genes = InputValue$geneinfo$right_gene,
      group.by = input$Heatmap_groupby,
      color.panel = input$Heatmap_colorpanel,
      font.size = input$PlotPar_fontsize
    )
  })
  
  output$Heatmap_plot.ui <- renderUI({
    plotOutput(
      "Heatmap_plot",
      width = 1.2 * InputValue$plot_width,
      height = 0.6 * InputValue$plot_height
    )
  })
  
  output$Heatmap_download <- downloadHandler(
    filename = function() {
      paste0(
        input$SelectData_dataset,
        "_",
        input$SelectData_normalization,
        "_HeatmapPlot_",
        stringi::stri_rand_strings(1, 10),
        ".pdf"
      )
    },
    content = function(file) {
      pdf(
        file,
        width = 1.2 * InputValue$plot_width / 80,
        height = 0.6 * InputValue$plot_height / 160
      )
      DoHeatmapPlot(
        plot.data = InputValue$plot_data,
        genes = InputValue$geneinfo$right_gene,
        group.by = input$Heatmap_groupby,
        color.panel = input$Heatmap_colorpanel,
        font.size = input$PlotPar_fontsize
      )
      dev.off()
    }
  )
  
  # In-silico FACS Plot----
  output$FACS_plot <- renderPlot({
    DoISFACSPlot(
      plot.data = InputValue$plot_data,
      genes = InputValue$geneinfo$right_gene,
      color.by = input$FACS_colorby,
      color.panel = input$FACS_colorpanel,
      dot.size = input$FACS_dotsize,
      x.cutoff = input$FACS_xcutoff,
      y.cutoff = input$FACS_ycutoff,
      font.size = input$PlotPar_fontsize
    )
  })
  
  output$FACS_plot.ui <- renderUI({
    plotOutput("FACS_plot",
               width = InputValue$plot_width,
               height = InputValue$plot_height)
  })
  
  output$FACS_download <- downloadHandler(
    filename = function() {
      paste0(
        input$SelectData_dataset,
        "_",
        input$SelectData_normalization,
        "_FACSPlot_",
        stringi::stri_rand_strings(1, 10),
        ".pdf"
      )
    },
    content = function(file) {
      pdf(file,
          width = InputValue$plot_width / 80,
          height = InputValue$plot_height / 80)
      DoISFACSPlot(
        plot.data = InputValue$plot_data,
        genes = InputValue$geneinfo$right_gene,
        color.by = input$FACS_colorby,
        color.panel = input$FACS_colorpanel,
        dot.size = input$FACS_dotsize,
        x.cutoff = input$FACS_xcutoff,
        y.cutoff = input$FACS_ycutoff,
        font.size = input$PlotPar_fontsize
      )
      dev.off()
    }
  )
  
  output$FACS_volcano_plot <- renderPlot({
    DoISFACSVolcanoPlot(
      limma.data = FACS_DEGenes$limma,
      dot.size = input$FACS_Volcano_dotsize,
      adj.P.Val.cutoff = input$FACS_Volcano_Pvalue_cutoff,
      logFC.cutoff = input$FACS_Volcano_logFC_cutoff,
      ngenes.labeled = input$FACS_Volcano_ngenes_labeled,
      font.size = input$FACS_Volcano_fontsize
    )
  })
  
  output$FACS_volcano_plot.ui <- renderUI({
    plotOutput("FACS_volcano_plot",
               width = InputValue$plot_width,
               height = InputValue$plot_height)
  })
  
  output$FACS_volcano_download <- downloadHandler(
    filename = function() {
      paste0(
        input$SelectData_dataset,
        "_",
        input$SelectData_normalization,
        "_Volcano_Plot_",
        stringi::stri_rand_strings(1, 10),
        ".pdf"
      )
    },
    content = function(file) {
      pdf(file,
          width = InputValue$plot_width / 80,
          height = InputValue$plot_height / 80)
      plot(
        DoISFACSVolcanoPlot(
          limma.data = FACS_DEGenes$limma,
          dot.size = input$FACS_Volcano_dotsize,
          adj.P.Val.cutoff = input$FACS_Volcano_Pvalue_cutoff,
          logFC.cutoff = input$FACS_Volcano_logFC_cutoff,
          ngenes.labeled = input$FACS_Volcano_ngenes_labeled,
          font.size = input$FACS_Volcano_fontsize
        )
      )
      dev.off()
    }
  )
  
  output$FACS_DEGDataTable_output <- DT::renderDataTable({
    DT::datatable(
      FACS_DEGenes$limma %>% arrange(adj.P.Val) %>% head(2000) %>% arrange(desc(logFC)) %>% mutate_if(is.character, as.factor),
      class = "stripe table-condensed",
      options = list(pageLength = 10, scrollX = TRUE)
    )
  })
  
  output$FACS_DEGDataTable_download <- downloadHandler(
    filename = function() {
      paste0(
        input$SelectData_dataset,
        "_",
        input$SelectData_normalization,
        "_",
        InputValue$geneinfo$right_gene,
        "_X",
        input$FACS_xcutoff,
        "_Y",
        input$FACS_ycutoff,
        "_",
        input$FACS_DEG_group1,
        "_",
        input$FACS_DEG_group2,
        "_",
        stringi::stri_rand_strings(1, 10),
        ".csv"
      )
    },
    content = function(file) {
      write.csv(FACS_DEGenes$limma, file, row.names = TRUE)
    }
  )
  
  # Metadata plot----
  output$Metadata_plot <- renderPlot({
    DoMetadataPlot(
      plot.data = InputValue$plot_data,
      plot.type = input$Metadata_plottype,
      group.by1 = input$Metadata_groupby1,
      color.by = input$Metadata_colorby,
      group.by2 = input$Metadata_groupby2,
      quantify = input$Metadata_quantify,
      coord.flip = input$Metadata_flip,
      show.cutoff = input$Metadata_cutoff,
      nrow = input$PlotPar_nrow,
      font.size = input$PlotPar_fontsize
    )
  })
  
  output$Metadata_plot.ui <- renderUI({
    plotOutput("Metadata_plot",
               width = InputValue$plot_width,
               height = InputValue$plot_height)
  })
  
  output$Metadata_download <- downloadHandler(
    filename = function() {
      paste0(
        input$SelectData_dataset,
        "_",
        input$SelectData_normalization,
        "_MetadataPlot_",
        stringi::stri_rand_strings(1, 10),
        ".pdf"
      )
    },
    content = function(file) {
      pdf(file,
          width = InputValue$plot_width / 80,
          height = InputValue$plot_height)
      plot(
        DoMetadataPlot(
          plot.data = InputValue$plot_data,
          plot.type = input$Metadata_plottype,
          group.by1 = input$Metadata_groupby1,
          color.by = input$Metadata_colorby,
          group.by2 = input$Metadata_groupby2,
          quantify = input$Metadata_quantify,
          coord.flip = input$Metadata_flip,
          show.cutoff = input$Metadata_cutoff,
          nrow = input$PlotPar_nrow,
          font.size = input$PlotPar_fontsize
        )
      )
      dev.off()
    }
  )
  
  # DataTable output----
  output$DataTable_output <- DT::renderDataTable({
    DT::datatable(
      InputValue$plot_data %>% mutate_if(is.character, as.factor),
      class = "stripe table-condensed",
      options = list(pageLength = 10, scrollX = TRUE)
    )
  })
  
  output$DataTable_download <- downloadHandler(
    filename = function() {
      paste0(
        input$SelectData_dataset,
        "_",
        input$SelectData_normalization,
        "_DataTablePlot_",
        stringi::stri_rand_strings(1, 10),
        ".csv"
      )
    },
    content = function(file) {
      write.csv(InputValue$plot_data, file, row.names = TRUE)
    }
  )
  
}

# Run the application----
shinyApp(ui = ui, server = server)
