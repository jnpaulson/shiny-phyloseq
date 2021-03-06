# d3 code section
################################################################################
# d3network uix
################################################################################
output$d3_uix_color <- renderUI({
  selectInput("color_d3", "Color Variable:",
              choices = vars(input$type_d3, TRUE, TRUE),
              selected = d3NetworkColorVar)
})
output$d3_uix_node_label <- renderUI({
  selectInput("d3_node_label", "Node Label (hover):",
              choices = vars(input$type_d3, TRUE, TRUE),
              selected = d3NodeLabelVar,
              multiple = TRUE)
})
################################################################################
# d3 interactive network graphic 
################################################################################
# Define global reactive distance matrix. 
# Re-calc only if method or plot-type change.
d3distReact <- reactive({
  idist = NULL
  try({idist <- distance(physeq(), method=input$dist_d3, type=input$type_d3)}, silent=TRUE)
  if(is.null(idist)){warning("d3dist: Could not calculate distance matrix with these settings.")}
  return(idist)
})  
calculate_links_data = reactive({
  LinksData = dist_to_edge_table(d3distReact(),
                                 input$dist_d3_threshold,
                                 c("Source", "target"))
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
# The d3Network output definition.
output$networkPlot <- renderPrint({
  d3Network::d3ForceNetwork(
    Links = calculate_links_data()$link, 
    Nodes = calculate_links_data()$node,
    Source = "Source", Target = "target",
    Value = "Distance",
    NodeID = "ShowLabels",
    Group = default_Source(input$color_d3),
    linkColour = input$d3_link_color,
    opacity = input$d3_opacity,
    zoom = FALSE, 
    standAlone = FALSE, 
    width = input$width_d3, height = input$height_d3,
    parentElement = "#networkPlot")
})
# Downloadable standalone HTML file.
content_d3 = function(file){
  d3Network::d3ForceNetwork(
    Links = calculate_links_data()$link, 
    Nodes = calculate_links_data()$node,
    Source = "Source",
    Target = "target",
    Value = "value",
    NodeID = "ShowLabels",
    Group = default_Source(input$color_d3),
    linkColour = input$d3_link_color,
    opacity = input$d3_opacity,
    zoom = FALSE, 
    standAlone = TRUE, 
    width = input$width_d3,
    height = input$height_d3,
    file = file
  )
}
output$downloadd3 <- downloadHandler(filename = function(){paste0("d3_", simpletime(), ".html")},
                                     content = content_d3)
#zoom = as.logical(input$d3_zoom),