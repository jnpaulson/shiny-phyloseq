################################################################################
# Animation Parameters
################################################################################
# The frame delay interval in milliseconds
interval=1000
# Whether or not to loop the value range, or stop after one pass
loop=TRUE
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
