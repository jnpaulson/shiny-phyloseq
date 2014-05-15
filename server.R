################################################################################
# Options, default settings, and load packages
################################################################################
# load packages
library("shiny"); packageVersion("shiny")
library("phyloseq"); packageVersion("phyloseq")
library("ggplot2"); packageVersion("ggplot2")
library("data.table"); packageVersion("data.table")
library("d3Network"); packageVersion("d3Network")
# Animation options
source("animation-parameters.R", local=FALSE)
# By default, the file size limit is 5MB. It can be changed by
# setting this option. Here we'll raise limit to 9MB.
options(shiny.maxRequestSize = 100*1024^2)
theme_set(theme_bw())
pal = "Set1"
scale_colour_discrete <- function(palname = pal, ...) {
  scale_colour_brewer(palette = palname, ...)
}
scale_fill_discrete <- function(palname = pal, ...) {
  scale_fill_brewer(palette = palname, ...)
}
################################################################################
# Begin Shiny Server definition.
################################################################################
shinyServer(function(input, output){
  ################################################################################
  # Data-Reactive UI Definitions.
  ################################################################################
  # Placeholder for future interactive filtering:
  physeq = reactive({
    return(prune_taxa(taxa_sums(phytcf) > OTUSumDefault, phytcf))
  })
  # Define data-reactive variable lists
  rankNames = reactive({
    rankNames = as.list(rank_names(physeq(), errorIfNULL=FALSE))
    names(rankNames) <- rankNames
    return(rankNames)
  })
  variNames = reactive({
    variNames = as.list(sample_variables(physeq(), errorIfNULL=FALSE))
    names(variNames) <- variNames
    return(variNames)
  })
  sampvarlist = reactive({c(variNames(), list("NULL"="NULL"))})
  observe({print(paste0("Variables (`sampvarlist()`): ", sampvarlist(), collapse=" "))})
  specvarlist = reactive({c(rankNames(), list("NULL"="NULL"))})
  observe({print(paste0("Variables (`specvarlist()`): ", specvarlist(), collapse=" "))})
  vars = reactive({c(rankNames(), variNames(), list("NULL"="NULL"))})
  observe({print(paste0("Variables (`vars()`): ", vars(), collapse=" "))})
  # A generic selectInput UI. Plan is to pass a reactive argument to `choices`.
  uivar = function(id, label="Variable:", choices, selected="NULL"){
    selectInput(inputId=id, label=label, choices=choices, selected=selected)
  }
  ################################################################################
  # d3network uix
  ################################################################################
  output$d3_uix_color <- renderUI({
    if(input$type_d3=="samples"){
      return(uivar("color_d3", "Color Variable:", sampvarlist(), d3NetworkColorVar))
    } else if(input$type_d3=="taxa"){
      return(uivar("color_d3", "Color Variable:", specvarlist(), d3NetworkColorVar))
    } else {
      # Some kind of fail, throw all variables up.
      return(uivar("color_d3", "Color Variable:", vars(), d3NetworkColorVar))
    }
  })
  output$d3_uix_node_label <- renderUI({
    # "d3ShowRanks"
    if(input$type_d3=="samples"){
      return(uivar("d3_node_label", "Node Label:", sampvarlist(), d3NetworkColorVar))
    } else if(input$type_d3=="taxa"){
      return(uivar("d3_node_label", "Node Label:", specvarlist(), d3NetworkColorVar))
    } else {
      # Some kind of fail, throw all variables up.
      return(uivar("d3_node_label", "Node Label:", vars(), d3NetworkColorVar))
    }
  })  
  ################################################################################
  # Plot Rendering Stuff.
  ################################################################################
  # Define a proportions-only version of input phyloseq object
  physeqProp = reactive({transform_sample_counts(physeq(), function(x){x / sum(x)})})
  # Define a dummy "failed plot" to return if render section cannot build valid plot.
  fail_gen = function(main="Graphic Fail.", 
                      subtext="Please change parameters and try again."){
    qplot(x=0, y=0, main=main) + 
      annotate("text", 0, 0, label=":-(",
               size=75, angle=270, hjust=0.5, vjust=0.35) +
      annotate("text", 0, 0, label=subtext, size=10, hjust=0.5, vjust=-7) +
      theme_bw() + 
      theme(panel.border=element_blank(), axis.line=element_blank(),
            axis.text=element_blank(), axis.ticks=element_blank())
  }
  failp = fail_gen()
  # Define a default controlled ggplot printing check for all print rendering
  shiny_phyloseq_print = function(p, f=failp){
    if(inherits(p, "ggplot")){
      # Check that rendering will work
      printout = NULL
      try(printout <- print(p), silent=TRUE)
      if(is.null(printout)){
        # If still NULL, the print-render failed,
        # otherwise print() would have returned a 'list'
        # Nothing was printed. Print the fail graphic in its place.
        print(f)
      }
    } else {
      print(f)
    }
  }    
  # Define generic function to access/clean variables
  # This especially converts "NULL" to NULL
  av = function(x){
    if( isTRUE(all.equal(x, "")) | isTRUE(all.equal(x, "NULL")) ){
      return(NULL)
    }
    return(x)
  }
  ################################################################################
  # Generate a d3 network graphic 
  ################################################################################
  # Define global reactive distance matrix. 
  # Re-calc only if method or plot-type change.
  gdist <- reactive({
    try({idist <- distance(physeq(), method=input$dist_d3, type=input$type_d3)}, silent=TRUE)
    if(is.null(idist)){warning("gdist: Could not calculate distance matrix with these settings.")}
    return(idist)
  })
  calculate_links_data = reactive({
    d3dist <- as.matrix(gdist())
    # Set duplicate entries and self-links to Inf
    d3dist[upper.tri(d3dist, diag = TRUE)] <- Inf
    # Create data.table.
    d3LinkNames = c("OTU", "OTU_target")
    LinksData = data.table(reshape2::melt(d3dist, varnames=d3LinkNames, as.is = TRUE))
    # Remove entries above the threshold
    # (This will also remove self-links and duplicate links)
    LinksData <- LinksData[value < input$dist_d3_threshold, ]
    # Don't sort yet, instead create mapping variable from OTU ID to link node ID
    # d3link nodes are numbered from 0.
    nodeUnion = union(LinksData$OTU, LinksData$OTU_target)
    d3lookup = (0:(length(nodeUnion)-1))
    names(d3lookup) <- nodeUnion
    # In-place replacement.
    LinksData[, OTU:=d3lookup[OTU]]
    LinksData[, OTU_target:=d3lookup[OTU_target]]
    # Order by the `d3lookup` node ID, in this case, the OTU label
    setkey(LinksData, OTU)
    # Create covariates table (taxa in this case)
    NodeData = data.frame(OTU=nodeUnion, tax_table(physeq())[nodeUnion, ], stringsAsFactors = FALSE)
    #NodeData$FullRank <- apply(NodeData[, d3ShowRanks], 1, paste0, collapse=" ")
    return(list(link=data.frame(LinksData), node=NodeData))
  })
  default_OTU = function(x){
    if(is.null(av(x))){
      return("OTU")
    } else {
      return(x)
    }
  }
  # The d3Network output definition.
  output$networkPlot <- renderPrint({
    d3Network::d3ForceNetwork(
      Links = calculate_links_data()$link, 
      Nodes = calculate_links_data()$node,
      Source = "OTU", Target = "OTU_target",
      Value = "value",
      NodeID = default_OTU(input$d3_node_label),
      Group = default_OTU(input$color_d3),
      opacity = input$d3_opacity,
      zoom = FALSE,
      standAlone = FALSE, 
      width = input$width_d3, height = input$height_d3,
      parentElement = "#networkPlot")
  })
})
