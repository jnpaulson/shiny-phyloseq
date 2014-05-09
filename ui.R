source("animation-parameters.R", local=FALSE)
# Type for distance/network/etc. Samples or Taxa
uitype = function(id="type", selected="taxa"){
  selectInput(inputId=id, label="Calculation: Samples or Taxa?",
              selected=selected,
              choices=list("Taxa"="taxa", "Samples"="samples"))
}
# ui for point size slider
uiptsz = function(id="size"){
  numericInput(inputId=id, label="Point Size:", min=1, value=7, step=2)
}
# ui for point opacity slider
uialpha = function(id="alpha"){
  sliderInput(inputId=id, label="Opacity:", min=0, max=1, value=1, step=0.1)
}
#   Function to reate ui for distance method selection
#   NOTE: not all distance methods are supported if "taxa" selected for type. 
#   For example, the UniFrac distance and DPCoA cannot be calculated for taxa-wise 
#   distances, because they use a taxa-wise tree as part of their calculation 
#   between samples, and there is no transpose-equivalent for this tree
uidist = function(id, selected="bray"){
  distlist = as.list(unlist(phyloseq::distance("list")))
  names(distlist) <- distlist
  return(selectInput(id, "Distance Method:", distlist, selected=selected))
}
# Whether to use proportions or counts
uicttype = function(id="uicttype"){
  radioButtons(inputId=id, label="Count Type",
               choices=c("Counts", "Proportions"),
               selected="Counts")
} 
################################################################################
# sbp of plot_network
################################################################################
sbp_net = sidebarPanel(uitype("type_net", "samples"),
                       uidist("dist_net"),
                       uiOutput("network_uix_color"),
                       uiOutput("network_uix_shape"),
                       sliderInput("dispdist", "Edge Distance Threshold:",
                                   animate=animationOptions(interval=interval, loop=loop),
                                   min=0.0,
                                   max=netdist,
                                   value=0.5*netdist,
                                   step=step),
                       uiptsz("size_net"), uialpha("alpha_net")
)
################################################################################
# Define each fluid page
################################################################################
# Define in a single function, a standard definition
make_fluidpage = function(fptitle="", sbp, outplotid){
  fluidPage(
    titlePanel(fptitle),
    sidebarLayout(
      sidebarPanel=sbp,
      mainPanel=mainPanel(
        textOutput("textdist"),
        tags$hr(),
        plotOutput(outplotid)
      )
    )
  )
}
netpage = make_fluidpage("", sbp_net, "network")
# Data I/O page
datapage = fluidPage(
  titlePanel(""),
  sidebarLayout(
    sidebarPanel(
      h4('Select Dataset'),
      uiOutput("phyloseqDataset"),
      tags$hr(),
      h4('Upload Local Data'),
      fileInput('file1', ""),
      tags$hr(),
      h4('microbio.me/qiime public data'),
      p('(requires phyloseq 1.9.5+)'),
      textInput("qiime_server_ID", "Identifier for QIIME server", value="NULL"),
      radioButtons("qiime_server_ext", "File Extension",
                   choices=list(".zip", ".tgz", ".tar.gz")),
      tags$hr(),
      p('Plot parameters:'),
      numericInput("dataset_count_threshold", "Count Threshold", value=3, min=0, step=1)
    ),
    mainPanel(
      plotOutput("library_sizes"),
      plotOutput("OTU_count_thresh_hist"),
      tags$hr(),
      p("Data Summary:"),
      htmlOutput('contents')
    )
  )
)
# Data Filter page
filterpage = fluidPage(
  titlePanel(""),
  sidebarLayout(
    sidebarPanel(
      p("  "),
      p('Filtering Parameters:'),
      tags$hr(),
      p("Subset Expressions:"),
      textInput("filter_subset_taxa_expr", "`subset_taxa()` Expression: (e.g. Phylum=='Firmicutes')", value="NULL"),
      textInput("filter_subset_samp_expr", "`subset_samples()` Expression: (e.g. SampleType %in% 'Feces')", value="NULL"),
      tags$hr(),
      p('Total Sums Filtering:'),
      numericInput("filter_sample_sums_threshold", "Sample Sums Count Threshold", value=0, min=0, step=1),
      numericInput("filter_taxa_sums_threshold", "OTU Sums Count Threshold", value=0, min=0, step=1),
      tags$hr(),
      p('kOverA OTU Filtering:'),
      numericInput("filter_kOverA_count_threshold", "`A` - The Count Value Threshold", value=0, min=0, step=1),
      uiOutput("filter_ui_kOverA_k")
    ),
    mainPanel(
      plotOutput("filter_summary_plot"),
      tags$hr(),
      p("Original Data:"),
      htmlOutput('filtered_contents0'),
      tags$hr(),
      p("Filtered Data:"),
      htmlOutput('filtered_contents'),
      tags$hr(),
      p("Sample Variables:"),
      textOutput('sample_variables'),
      tags$hr(),
      p("Taxonomic Ranks:"),
      textOutput('rank_names')      
    )
  )
)
################################################################################
# Define the full user-interface, `ui`
################################################################################
ui = navbarPage("Shiny + phyloseq",
                tabPanel("Select Dataset", datapage),
                tabPanel("Filter", filterpage),
                tabPanel("Network", netpage)
)
shinyUI(ui)
