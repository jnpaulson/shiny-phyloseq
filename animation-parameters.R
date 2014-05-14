################################################################################
# Animation Parameters
################################################################################
# The frame delay interval in milliseconds
interval=1000
# Whether or not to loop the value range, or stop after one pass
loop=TRUE
################################################################################
# Filter Default Parameters
################################################################################
# sample_sums_threshold
# input$filter_sample_sums_threshold
SampleSumDefault = 1000
# taxa_sums threshold 
# input$filter_taxa_sums_threshold
OTUSumDefault = 250
################################################################################
# Network Threshold Animation Default Parameters
################################################################################
# The maximum distance with which the network structure is based. 
netdist = 0.20
# The distance step size to use in the animation.
step = netdist/15
netThreshColorVariableDefault = "mouse.number"
netThreshShapeVariableDefault = "infected"
netThreshDistanceMethod = "jaccard"
################################################################################
# Ordination Default Parameters
################################################################################
# Default Selections (have no impact if not in choices list)
ordinationTimeVariableDefault = "timeOrderMouse"
ordinationColorVariableDefault = "mouse.number"
ordinationShapeVariableDefault = "infected"
################################################################################
# d3 default Default Parameters
################################################################################
LinkDistThreshold = 0.6
d3DefaultDistance = "jaccard"
d3NetworkColorVar = "Phyllum"
d3ShowRanks = c("Phyllum", "Genus", "Species")
################################################################################