library(igraph)
library(ggplot2)
n <- 83
m <- 734
transitivity_distribution_macro <- vector()
for(i in 1:10000) {
  gn <- sample_gnm(n, m, directed = FALSE, loops = FALSE)
  transitivity_distribution_macro[i] = transitivity(gn)
}
hist(transitivity_distribution_macro)

macronet <- read.delim('MICTools_Network_Macro_PresenceAbsence/MIC_tools_output/strength_MACRO_PA.txt', header = TRUE, sep = "\t")
macronet <- macronet[,2:7]
g=graph.data.frame(macronet,directed = FALSE)
print(transitivity(g,type = "global"))

n <- 619
m <- 40126
transitivity_distribution_micro <- vector()
for(i in 1:10000) {
  gn <- sample_gnm(n, m, directed = FALSE, loops = FALSE)
  transitivity_distribution_micro[i] = transitivity(gn)
}
hist(transitivity_distribution_micro)

micronet <- read.delim('MICTools_Network_Micro_PresenceAbsence/MIC_tools_output/strength_MICRO_PresenceAbsence_Network.txt', header = TRUE, sep = "\t")
micronet <- micronet[,2:7]
g=graph.data.frame(micronet,directed = FALSE)
print(transitivity(g,type = "global"))

tmicro <- data.frame(transitivity = transitivity_distribution_micro)
tmacro <- data.frame(transitivity = transitivity_distribution_macro)
tmicro$type <- 'micro'
tmacro$type <- 'macro'
combo <- rbind(tmicro, tmacro)
ggplot(combo, aes(transitivity, fill = type)) + geom_density(alpha = 0.2)