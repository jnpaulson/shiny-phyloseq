################################################################################
# Options, default settings, and load packages
################################################################################
# Animation options
source("animation-parameters.R", local=FALSE)
# By default, the file size limit is 5MB. It can be changed by
# setting this option. Here we'll raise limit to 9MB.
options(shiny.maxRequestSize = 100*1024^2)
# load packages
library("shiny"); packageVersion("shiny")
library("phyloseq"); packageVersion("phyloseq")
library("ggplot2"); packageVersion("ggplot2")
library("data.table"); packageVersion("data.table")
library("d3Network"); packageVersion("d3Network")
theme_set(theme_bw())
pal = "Set1"
scale_colour_discrete <- function(palname = pal, ...) {
  scale_colour_brewer(palette = palname, ...)
}
scale_fill_discrete <- function(palname = pal, ...) {
  scale_fill_brewer(palette = palname, ...)
}
################################################################################
# Included Data
# Define the named list of datasets to choose from
################################################################################
includedDatasets = c("GlobalPatterns", "enterotype", "esophagus", "soilrep")
data(list=includedDatasets)
datalist = list(GlobalPatterns=GlobalPatterns, 
                enterotype=enterotype,
                esophagus=esophagus,
                soilrep=soilrep)
# The mouse study data:
load("data/phytcf.RData")
# One-Off fix infected variable
sample_data(phytcf)$mouse.type <- factor(get_variable(phytcf, "mouse.type"))
sample_data(phytcf)$infected <- as.logical(get_variable(phytcf, "infected"))
sample_data(phytcf)$mouse.number <- factor(get_variable(phytcf, "mouse.number"))
# Add phytcf as 'mouse' in the available data list.
if(inherits(phytcf, "phyloseq")){
  datalist <- c(list(mouse=phytcf), datalist)
}
################################################################################
# Begin Shiny Server definition.
################################################################################
shinyServer(function(input, output){
  ################################################################################
  # Define the available phyloseq datasets for plotting.
  ################################################################################
  get_qiime_data = reactive({
    qiime_data = NULL
    if(!is.null(av(input$qiime_server_ID))){
      if( !is.na(as.integer(input$qiime_server_ID)) ){
        observe({print(paste0("Attempting integer ID import: ", input$qiime_server_ID))})
        zipftp = as(isolate({input$qiime_server_ID}), "integer")
        studyname = input$qiime_server_ID
      } else {
        observe({print(paste0("Attempting character ID import: ", input$qiime_server_ID))})
        zipftp = as(isolate({input$qiime_server_ID}), "character")
        studyname = gsub("\\_split\\_.+$", "", basename(zipftp))
      }
      observe({print(paste0("Extension Chosen: ", input$qiime_server_ext))})
      trash = try({qiime_data <- microbio_me_qiime(zipftp, ext=input$qiime_server_ext)}, silent=TRUE)
    }
    if(inherits(qiime_data, "phyloseq")){
      qiime_data <- list(qiime_data)
      names(qiime_data) <- studyname
      datalist <<- c(qiime_data, datalist)
    } else {
      observe({print("Attempt made to access qiime server data, but didn't work this pass...")})
    }
    return(NULL)
  })
  get_loaded_data = reactive({
    if(!is.null(input$file1$name)){
      # Added uploaded data, if provided, and it is phyloseq-class.
      objectNames = load(input$file1$datapath)
      loadedObjects = mget(objectNames)
      arePhyloseq = sapply(loadedObjects, inherits, "phyloseq")
      if(any(arePhyloseq)){
        loadedObjects <- loadedObjects[which(arePhyloseq)]
      } else {
        loadedObjects <- NULL
      }
      datalist <<- c(loadedObjects, datalist)
      observe({print(paste("Available objects in datalist:", names(datalist), collapse=", "))})
    }
    return(NULL)
  })
  output$phyloseqDataset <- renderUI({
    # Expect the side-effect of these two functions to be to add
    # elements to the datalist, if appropriate
    get_loaded_data()
    get_qiime_data()
    return(radioButtons("physeqSelect", "Available Datasets:", names(datalist)))
  })
  get_phyloseq_data = reactive({
    ps0 = NULL
    if(!is.null(input$physeqSelect)){
      if(input$physeqSelect %in% names(datalist)){
        ps0 <- datalist[[input$physeqSelect]]
      }
    }
    observe({print(ps0)})
    if(inherits(ps0, "phyloseq")){
      return(ps0)
    } else {
      observe({print("ps0 is NULL in get_phyloseq_data()")})
      return(NULL)
    }
  })
  ################################################################################
  # Filtering
  ################################################################################
  physeq = reactive({
    ps0 = get_phyloseq_data()
    if(inherits(ps0, "phyloseq")){
      observe({print(paste("filter_kOverA_count_threshold:", input$filter_kOverA_count_threshold))})
      observe({print(paste("filter_kOverA_sample_threshold:", input$filter_kOverA_sample_threshold))})
      observe({print(paste("filter_sample_sums_threshold:", input$filter_sample_sums_threshold))})
      observe({print(paste("filter_taxa_sums_threshold:", input$filter_taxa_sums_threshold))})
      observe({print(paste("filter_subset_taxa_expr:", input$filter_subset_taxa_expr))})
      observe({print(paste("filter_subset_samp_expr:", input$filter_subset_samp_expr))})
      # Expression filters
      if( !is.null(av(input$filter_subset_taxa_expr)) ){
        ps0 = eval(parse(text=paste0("subset_taxa(ps0, ", input$filter_subset_taxa_expr, ")")))
        observe({print("subset_taxa...")})
        observe({print(ps0)})
      }
      if( !is.null(av(input$filter_subset_samp_expr)) ){
        ps0 = eval(parse(text=paste0("subset_samples(ps0, ", input$filter_subset_samp_expr, ")")))
        observe({print("subset_samples...")})
        observe({print(ps0)})
      }
      if( input$filter_taxa_sums_threshold > 0 ){
        # OTU sums filter
        ps0 <- prune_taxa({taxa_sums(ps0) > input$filter_taxa_sums_threshold}, ps0)
        observe({print("prune OTUs...")})
        observe({print(ps0)})
      }
      if( input$filter_sample_sums_threshold > 0 ){
        # Sample sums filtering
        ps0 <- prune_samples({sample_sums(ps0) > input$filter_sample_sums_threshold}, ps0)
        observe({print("prune samples...")})
        observe({print(ps0)})
      }
      if(inherits(input$filter_kOverA_sample_threshold, "numeric")){
        if(input$filter_kOverA_sample_threshold > 1){
          # kOverA OTU Filtering
          flist = genefilter::filterfun(
            genefilter::kOverA(input$filter_kOverA_sample_threshold,
                               input$filter_kOverA_count_threshold, na.rm=TRUE)
          )
          ps0 <- filter_taxa(ps0, flist, prune=TRUE)
        }
      }
      return(ps0)
    } else {
      return(NULL)
    }
  })
  # kOverA `k` Filter UI
  maxSamples = reactive({
    # Create logical indicated the samples to keep, or dummy logical if nonsense input
    if(inherits(get_phyloseq_data(), "phyloseq")){
      return(nsamples(get_phyloseq_data()))
    } else {
      # Dummy response.
      return(NULL)
    }
  })
  output$filter_ui_kOverA_k <- renderUI({
    numericInput("filter_kOverA_sample_threshold",
                "`k` - Number of Samples that Must Exceed `A`",
                min=0, max=maxSamples(), value=0, step=1)    
  })  
  output_phyloseq_print_html <- reactive({
    HTML(paste0(capture.output(print(get_phyloseq_data())), collapse=" <br/> "))
  })
  output$contents <- renderUI({
    output_phyloseq_print_html()
  })
  output$filtered_contents0 <- renderUI({
    output_phyloseq_print_html()
  })
  output$filtered_contents <- renderUI({
    HTML(paste0(capture.output(print(physeq())), collapse=" <br/> "))
  })
  # Generic Function for plotting marginal histograms
  sums_hist = function(thesums=NULL, xlab="", ylab=""){
    if(is.null(thesums)){
      p = qplot(0)
    } else {
      p = ggplot(data.frame(sums=thesums), aes(x=sums))
      p = p + geom_histogram()
      p = p + xlab(xlab) + ylab(ylab) 
      p = p + scale_x_log10(labels = scales::comma)
    }
    return(p)
  }
  lib_size_hist = reactive({
    xlab = "Number of Reads (Counts)"
    ylab = "Number of Libraries"
    return(sums_hist(sample_sums(get_phyloseq_data()), xlab, ylab))
  })
  otu_sum_hist = reactive({
    xlab = "Number of Reads (Counts)"
    ylab = "Number of OTUs"
    return(sums_hist(taxa_sums(get_phyloseq_data()), xlab, ylab))    
  })
  output$library_sizes <- renderPlot({
    if(inherits(get_phyloseq_data(), "phyloseq")){
      libtitle = "Histogram of Library Sizes in Selected Data"
      p1 = lib_size_hist() + ggtitle(libtitle)
      otusumtitle = "Histogram of OTU total counts in Selected Data"
      p2 = otu_sum_hist() + ggtitle(otusumtitle)
      gridExtra::grid.arrange(p1, p2, ncol=2)
    } else {
      libfailtext = "Press the `Load Selection` button \n to load/refresh data."
      print(qplot(x=0, y=0, main="") + xlim(-1, 1) + ylim(-1, 1) + 
              geom_segment(aes(x=0, y=0, xend=-0.5, yend=0.15), size=3,
                           arrow=grid::arrow(length=grid::unit(0.5, "cm"))) +
              annotate("text", 0, 0, label=libfailtext, size=12, hjust=0.5, vjust=-1) +
              theme(panel.border=element_blank(), axis.line=element_blank(),
                    axis.text=element_blank(), axis.ticks=element_blank())
      )
    }
  })
  output$OTU_count_thresh_hist <- renderPlot({
    ps0 = get_phyloseq_data()
    if(inherits(get_phyloseq_data(), "phyloseq")){
      mx = as(otu_table(ps0), "matrix")
      if(!taxa_are_rows(ps0)){mx <- t(mx)}
      thresh = input$dataset_count_threshold
      df = data.frame(x=apply(mx, 1, function(x, thresh){sum(x>thresh)}, thresh))
      p = ggplot(df, aes(x=x)) + geom_histogram()
      p = p + xlab("Number of Samples with Count Above Threshold") + ylab("Number of OTUs")
      p = p + ggtitle(paste("Histogram of OTUs Observed More Than", thresh, "Times"))
      return(shiny_phyloseq_print(p))
    } else {
      return(qplot(0))
    }
  })
  output$sample_variables <- renderText({return(
    paste0(sample_variables(get_phyloseq_data(), errorIfNULL=FALSE), collapse=", ")
  )})
  output$rank_names <- renderText({return(
    paste0(rank_names(get_phyloseq_data(), errorIfNULL=FALSE), collapse=", ")
  )})
  output$filter_summary_plot <- renderPlot({
    plib0 = lib_size_hist() + ggtitle("Original Data")
    potu0 = otu_sum_hist() + ggtitle("Original Data")
    if(inherits(physeq(), "phyloseq")){
      potu1 = sums_hist(taxa_sums(physeq()), xlab = "Number of Reads (Counts)",
                       ylab = "Number of OTUs"
                       ) + 
        ggtitle("Filtered Data")
      plib1 = sums_hist(sample_sums(physeq()), xlab = "Number of Reads (Counts)",
                        ylab = "Number of Libraries"
                        ) + 
        ggtitle("Filtered Data")
    } else {
      potu1 = plib1 = failp
    }
    gridExtra::grid.arrange(plib0, potu0, plib1, potu1, ncol=2, 
                            main="Histograms: Before and After Filtering")
  })
  ################################################################################
  # Data-Reactive UI Definitions.
  ################################################################################
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
  # plot_network() ui
  ################################################################################
  output$network_uix_color <- renderUI({
    if(input$type_net=="samples"){
      return(uivar("color_net", "Color Variable:", sampvarlist(), netThreshColorVariableDefault))
    } else if(input$type_net=="taxa"){
      return(uivar("color_net", "Color Variable:", specvarlist(), netThreshColorVariableDefault))
    } else {
      # Some kind of fail, throw all variables up.
      return(uivar("color_net", "Color Variable:", vars()))
    }
  })
  output$network_uix_shape <- renderUI({
    if(input$type_net=="samples"){
      return(uivar("shape_net", "Shape Variable:", sampvarlist(), netThreshShapeVariableDefault))
    } else if(input$type_net=="taxa"){
      return(uivar("shape_net", "Shape Variable:", specvarlist(), netThreshShapeVariableDefault))
    } else {
      # Some kind of fail, throw all variables up.
      return(uivar("shape_net", "Shape Variable:", vars()))
    }
  })
  ################################################################################
  # ordination uix
  ################################################################################
  output$ord_uix_color <- renderUI({
    if(input$type_ord=="samples"){
      return(uivar("color_ord", "Color Variable:", sampvarlist(), ordinationColorVariableDefault))
    } else if(input$type_ord=="taxa"){
      return(uivar("color_ord", "Color Variable:", specvarlist(), ordinationColorVariableDefault))
    } else {
      # Some kind of fail, throw all variables up.
      return(uivar("color_ord", "Color Variable:", vars()))
    }
  })
  output$ord_uix_shape <- renderUI({
    if(input$type_ord=="samples"){
      return(uivar("shape_ord", "Shape Variable:", sampvarlist(), ordinationShapeVariableDefault))
    } else if(input$type_ord=="taxa"){
      return(uivar("shape_ord", "Shape Variable:", specvarlist(), ordinationShapeVariableDefault))
    } else {
      # Some kind of fail, throw all variables up.
      return(uivar("shape_ord", "Shape Variable:", vars()))
    }
  })  
  # Time variable for ordination animation
  output$ord_uix_timevar <- renderUI({
    if(input$type_ord=="samples"){
      return(uivar("timevar_ord", "Time Variable:", sampvarlist(), ordinationTimeVariableDefault))
    } else if(input$type_ord=="taxa"){
      return(uivar("timevar_ord", "Time Variable:", specvarlist(), ordinationTimeVariableDefault))
    } else {
      # Some kind of fail... dummy ui
      return(selectInput("timevar_ord", "Time Variable", choices = "NULL"))
    }
  })
  get_timevar_ord = reactive({
    if(input$type_ord=="taxa" | is.null(av(input$timevar_ord))){
      return(NULL)
    }
    return(as.numeric(get_variable(physeq_ord(), input$timevar_ord)))
  })
  output$ord_uix_timeslider <- renderUI({
    observe({print(paste0(get_timevar_ord(), collapse=", "))})
    if(is.null(av(input$timevar_ord)) | is.null(get_timevar_ord())){
      return(selectInput("timeslider_ord", "Time Variable", choices = "NULL"))
    }
    timerange = range(get_timevar_ord())
    return(sliderInput("timeslider_ord", label = paste0("Animate: ", input$timevar_ord),
                       min=timerange[1], 
                       max=timerange[2],
                       value=timerange[1],
                       animate=animationOptions(interval=interval, loop=loop))
    )
  })
  # Dynamic filtering for ordination
  output$ord_uix_subsetvar <- renderUI({
    if(input$type_ord=="samples"){
      return(uivar("subsetvar_ord", "Subset Variable:", sampvarlist()))
    } else if(input$type_ord=="taxa"){
      return(uivar("subsetvar_ord", "Subset Variable:", specvarlist()))
    } else {
      # Some kind of fail, throw all variables up.
      return(uivar("subsetvar_ord", "Subset Variable:", vars()))
    }
  })
  output$ord_uix_selectelem <- renderUI({
    observe({print(paste0("ord_uix_selectelem input$subsetvar_ord: ", input$subsetvar_ord))})
    if(is.null(av(input$subsetvar_ord))){
      return(selectInput("elem_select", "Subset by Elements", 
                  choices="NULL", selected="NULL", multiple=TRUE))
    }
    if(input$type_ord=="samples"){
      choices = unique(get_variable(physeq(), input$subsetvar_ord))
      return(
        selectInput("elem_select", "Subset by Elements", 
                    choices=choices, selected=choices, multiple=TRUE)
      )
    } else if(input$type_ord=="taxa"){
      choices = get_taxa_unique(physeq(), input$subsetvar_ord)
      return(
        selectInput("elem_select", "Subset by Elements",
                    choices=choices, selected=choices, multiple=TRUE)
      )
    }
  })
  physeq_ord = reactive({
    if(is.null(av(input$subsetvar_ord))){
      return(physeq())
    }
    if(input$type_ord=="samples"){
      keepSamples = get_variable(physeq(), input$subsetvar_ord) %in% input$elem_select
      return(prune_samples(keepSamples, physeq()))
    } else if(input$type_ord=="taxa"){
      keepTaxa = as(tax_table(physeq()), "matrix")[, input$subsetvar_ord] %in% input$elem_select
      return(prune_taxa(keepTaxa, physeq()))
    }
  })
  ################################################################################
  # d3network uix
  ################################################################################
  output$d3_uix_color <- renderUI({
    if(input$type_d3=="samples"){
      return(uivar("color_d3", "Color Variable:", choices = c("Sample", variNames()), "mouse.number"))
    } else if(input$type_d3=="taxa"){
      return(uivar("color_d3", "Color Variable:", 
                   c(rankNames(), list(OTU="OTU")), 
                   d3NetworkColorVar))
    } else {
      # Some kind of fail, throw all variables up.
      return(uivar("color_d3", "Color Variable:", vars(), d3NetworkColorVar))
    }
  })
  output$d3_uix_node_label <- renderUI({
    if(input$type_d3=="samples"){
      return(selectInput("d3_node_label", "Node Label (hover):", 
                         choices = c("Sample", variNames()),
                         selected = c("Sample", "infected"), 
                         multiple = TRUE))
    } else if(input$type_d3=="taxa"){
      return(selectInput("d3_node_label", "Node Label (hover):",
                         choices = c(rankNames(), list(OTU="OTU")), 
                         multiple = TRUE,
                         selected = d3NodeLabelVar))
    } else {
      # Some kind of fail, throw all variables up.
      return(uivar("d3_node_label", "Node Label:", vars(), d3NetworkColorVar))
    }
  }) 
  ################################################################################
  # Plot Rendering Stuff.
  ################################################################################
  default_Source = function(x){
    if(is.null(av(x))){
      if(input$type_d3=="taxa"){
        return("OTU")
      } else {
        return("Sample")
      }
    } else {
      return(x)
    }
  }
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
  # Generate a network plot 
  ################################################################################
  # The reactive value version of the network.
  # Only returns distance matrix, regardless of distance-method argument
  distonly <- reactive({
    idist = NULL
    try({idist <- distance(physeq(), method=input$dist_net, type=input$type_net)}, silent=TRUE)
    if(is.null(idist)){warning("distonly: Could not calculate distance matrix with these settings.")}
    # rescale the distance matrix to be [0, 1]
    idist <- idist / max(idist, na.rm=TRUE)
    idist <- idist - min(idist, na.rm=TRUE)
    return(idist)
  })
  ig <- reactive({
    make_network(physeq(),
                 type=input$type_net,
                 distance=distonly(),
                 max.dist=netdist,
                 keep.isolates=FALSE)
  })
  initial_plot_network = reactive({
    plot_network(ig(), physeq(),
                 type=input$type_net,
                 color=isolate(input$color_net),
                 shape=isolate(input$shape_net),
                 point_size=isolate(input$size_net),
                 alpha=isolate(input$alpha_net)
    )
  })
  get_edge_df = reactive({
    p = initial_plot_network()
    whichEdge = which(sapply(p$layers, function(x){x$geom$objname=="line"}))
    edgeDF0 = p$layers[[whichEdge]]$data
    # Add the distance associated with each edge entry
    edgeDF0 = plyr::ddply(edgeDF0, "id", function(df, dmat){
      df$dist <- dmat[df$value[1], df$value[2]]
      return(df)
    }, dmat = as.matrix(distonly()))
    return(edgeDF0)
  })
  update_plot_network = reactive({
    p = initial_plot_network()
    # New edge layer
    # Define the layer that contains the edges and vertices
    whichEdge = which(sapply(p$layers, function(x){x$geom$objname=="line"}))
    # Define the edge data frame
    edgeDF0 = get_edge_df()
    # Subset newEdgeDF according to max allowed distance
    newEdgeDF = edgeDF0[edgeDF0$dist <= input$dispdist, ]
    newEdgeMap = aes_string(x="x", y="y", group="id", colour=input$color_net)
    newEdgeLayer = geom_line(mapping=newEdgeMap, data=newEdgeDF, alpha=input$alpha_net)
    p$layers[[whichEdge]] <- newEdgeLayer
    # New vertex layer
    whichVert = which(sapply(p$layers, function(x){x$geom$objname=="point"}))
    # Subset
    newVertDF = p$data[as.character(p$data$value) %in% as.character(newEdgeDF$value), ]
    newVertMap = aes_string(x="x", y="y", colour=input$color_net, shape=input$shape_net)
    newVertLayer = geom_point(mapping=newVertMap, data=newVertDF,
                              size=input$size_net, alpha=input$alpha_net)
    p$layers[[whichVert]] <- newVertLayer
    # New label layer
    # Re-define the subset 
    whichLabel = which(sapply(p$layers, function(x){x$geom$objname=="text"}))
    p$layers[[whichLabel]]$data <- newVertDF
    # Updates legend labels.
    p = update_labels(p, list(colour=input$color_net))
    p = update_labels(p, list(shape=input$shape_net))
    # Fix the coordinate ranges based on the original.
    p = p + xlim(I(range(edgeDF0$x, na.rm=TRUE, finite=TRUE)))
    p = p + ylim(I(range(edgeDF0$y, na.rm=TRUE, finite=TRUE)))
    p = p + ggtitle(paste("Edge Distance Threshold: ", round(input$dispdist, digits=3)))
    return(p)
  })
  output$network <- renderPlot({
    shiny_phyloseq_print(update_plot_network())
  }, width=700, height=500)
  ################################################################################
  # Ordination section
  ################################################################################
  get_formula <- reactive({
    if(is.null(av(input$formula)) | input$formula=="NULL"){
      return(NULL)
    } else {
      return(as.formula(input$formula))
    }
  })
  # Define global reactive distance matrix. Will re-calc if method or plot-type change.
  gdist <- reactive({
    if(input$dist_ord %in% distance("list")$vegdist){
      return(input$dist_ord)
    } else {
      idist = NULL
      try({idist <- distance(physeq_ord(), method=input$dist_ord, type=input$type_ord)}, silent=TRUE)
      if(is.null(idist)){warning("gdist: Could not calculate distance matrix with these settings.")}
      return(idist)
    }
  })
  # Define reactive ordination access
  get_ord = reactive({
    ordinate(physeq_ord(), method=input$ord_method, distance=gdist(), formula=get_formula())
  })
  make_ord_plot = reactive({
    p1 = NULL
    try(p1 <- plot_ordination(physeq_ord(), get_ord(), type=input$type_ord), silent=TRUE)
    return(p1)
  })
  # ordination plot definition
  output$ordination <- renderPlot({
    p1 = make_ord_plot()
    if(inherits(p1, "ggplot")){
      # define the full xlim, ylim
      xlim0 = range(p1$data[[as.character(p1$mapping$x)]])
      ylim0 = range(p1$data[[as.character(p1$mapping$y)]])
      # Reactive size and opacity.
      p1$layers[[1]]$geom_params$size <- av(input$size_ord)
      p1$layers[[1]]$geom_params$alpha <- av(input$alpha_ord)
      # Reactive color and shape
      if(!is.null(av(input$color_ord))){
        p1$mapping$colour <- as.symbol(av(input$color_ord))
      }
      if(!is.null(av(input$shape_ord))){
        p1$mapping$shape  <- as.symbol(av(input$shape_ord))
      }
      # facetting 
      if(!is.null(c(av(input$color_ord), av(input$shape_ord)))){
        facetForm = as.formula(paste0("~", paste0(c(av(input$color_ord), av(input$shape_ord)), collapse = " + ")))
        p1 = p1 + facet_wrap(facets = facetForm)
      }
      # optional time-variable animation with path.
      if(!is.null(get_timevar_ord()) & !is.null(av(input$timeslider_ord))){
        # Subset to just the data up to the current time variable
        p1$data <- p1$data[as.numeric(p1$data[[input$timevar_ord]]) <= input$timeslider_ord, ]
        if(!is.null(c(av(input$color_ord), av(input$shape_ord)))){
          currentDF = plyr::ddply(p1$data, .variables=c(av(input$color_ord), av(input$shape_ord)), function(x){
            x[which.max(as.numeric(x[[input$timevar_ord]])), , drop=FALSE]
          })
          p1 = p1 + geom_path()
          p1 = p1 + geom_point(data=currentDF, 
                               color="black",
                               shape="+",
                               size=0.9*av(input$size_ord))
        }
        # Add a marker to indicate the "current" point(s)
        p1$data[as.numeric(p1$data[[input$timevar_ord]]) <= input$timeslider_ord, ]
        # Update range to keep full range at every time point.
        p1 = p1 + xlim(xlim0)
        p1 = p1 + ylim(ylim0)
        # Display current time in title
        title_ord = paste0(input$timevar_ord, ": ", round(as.numeric(input$timeslider_ord), 2))
        p1 = p1 + ggtitle(title_ord)
      }
      # Updates legend labels.
      p1 = update_labels(p1, list(colour=input$color_ord))
      p1 = update_labels(p1, list(shape=input$shape_ord))
    }
    shiny_phyloseq_print(p1)
  }, width=800, height=650)
  ################################################################################
  # d3 interactive network graphic 
  ################################################################################
  # Define global reactive distance matrix. 
  # Re-calc only if method or plot-type change.
  d3distReact <- reactive({
    try({idist <- distance(physeq(), method=input$dist_d3, type=input$type_d3)}, silent=TRUE)
    if(is.null(idist)){warning("d3dist: Could not calculate distance matrix with these settings.")}
    return(idist)
  })  
  calculate_links_data = reactive({
    d3dist <- as.matrix(d3distReact())
    # Set duplicate entries and self-links to Inf
    d3dist[upper.tri(d3dist, diag = TRUE)] <- Inf
    # Create data.table.
    d3LinkNames = c("Source", "target")
    LinksData = data.table(reshape2::melt(d3dist, varnames=d3LinkNames, as.is = TRUE))
    # Remove entries above the threshold
    # (This will also remove self-links and duplicate links)
    LinksData <- LinksData[value < input$dist_d3_threshold, ]
    # Rescale remaining links
    LinksData[, value:=(0.1+input$d3_link_scale*(value-min(value))/max(value))]
    # Don't sort yet, instead create mapping variable from Source ID to link node ID
    # d3link nodes are numbered from 0.
    nodeUnion = union(LinksData$Source, LinksData$target)
    d3lookup = (0:(length(nodeUnion)-1))
    names(d3lookup) <- nodeUnion
    # In-place replacement.
    LinksData[, Source:=d3lookup[Source]]
    LinksData[, target:=d3lookup[target]]
    # Order by the `d3lookup` node ID, in this case, the Source label
    setkey(LinksData, Source)
    # Create covariates table (taxa in this case)
    if(input$type_d3 == "taxa"){
      NodeData = data.frame(OTU=nodeUnion, tax_table(physeq())[nodeUnion, ], stringsAsFactors = FALSE)
    } else {
      NodeData = data.frame(Sample=nodeUnion, sample_data(physeq())[nodeUnion, ], stringsAsFactors = FALSE)      
    }
    NodeData$ShowLabels <- apply(NodeData[, input$d3_node_label, drop=FALSE], 1, paste0, collapse="; ")
    return(list(link=data.frame(LinksData), node=NodeData))
  })
  # The d3Network output definition.
  output$networkPlot <- renderPrint({
    d3Network::d3ForceNetwork(
      Links = calculate_links_data()$link, 
      Nodes = calculate_links_data()$node,
      Source = "Source", Target = "target",
      Value = "value",
      NodeID = "ShowLabels",
      Group = default_Source(input$color_d3),
      linkColour = input$d3_link_color,
      opacity = input$d3_opacity,
      zoom = TRUE, #zoom = as.logical(input$d3_zoom),
      standAlone = FALSE, 
      width = input$width_d3, height = input$height_d3,
      parentElement = "#networkPlot")
  })
})
