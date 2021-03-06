################################################################################
# sbp of plot_ordination 
################################################################################
ordtypelist = as.list(phyloseq::plot_ordination("list"))
names(ordtypelist) <- c("Samples", "Species", "Biplot", "Split Plot", "Scree Plot")
sbp_ord = sidebarPanel(
  uitype("type_ord", "samples"),
  uidist("dist_ord"),
  uiOutput("ord_uix_color"),
  uiOutput("ord_uix_shape"),
  selectInput("ord_method", "Ordination Method:", ordlist, selected="DCA"),
  selectInput("ord_plot_type", "Ordination Plot Type:", ordtypelist), 
  textInput("formula", "Ordination Constraint Formula", value="NULL"),
  h4('Figure Details'),
  uiptsz("size_ord"),
  uialpha("alpha_ord"),
  uipal("pal_ord"),
  uitheme("theme_ord"),
  tags$hr(),
  h4('Figure Dimensions'),
  numericInput("width_ord", "Figure Width (inches)", 8, 1, 100, 1),
  numericInput("height_ord", "Figure Height (inches)", 8, 1, 100, 1),
  graphicTypeUI("downtype_ord"),
  downloadButton('downloadOrdination', 'Download Graphic')
)
################################################################################
# ordpage = make_fluidpage("", sbp_ord, "ordination")
dummyText = paste0(rep("  \n  ", 7), collapse = "")
ordpage = fluidPage(
  titlePanel(""),
  sidebarLayout(
    sidebarPanel=sbp_ord,
    mainPanel=mainPanel(
      plotOutput("ordination"),
      tags$hr(dummyText),
      tags$br(dummyText),
      tags$hr(dummyText),
      tags$br(dummyText),
      tags$hr(dummyText),
      plotOutput("scree_ord")
    )
  )
)
