library(igraph)
library(reshape2)
library(ggplot2)
micronet <- read.delim('MICTools_Network_Micro_PresenceAbsence/MIC_tools_output/strength_MICRO_PresenceAbsence_Network.txt', header = TRUE, sep = "\t")
micronet <- micronet[,2:7]
micronet$weight <- micronet$PearsonR
g=graph.data.frame(micronet,directed = FALSE)
dg <- decompose.graph(g)
g = dg[[3]]
microlist <- read.csv('2.Tabla_Indval_Final_Filter100_AllASVs_MICRO.csv', header = TRUE)
rownames(microlist) <- microlist[,1]
n <- V(g)$name
microlist = microlist[intersect(n,rownames(microlist)),] # remove species not in network
membership <- rep(NA, 619)
membership[(microlist$IntertidalZone == "LowZone") & (microlist$p.value <= 0.05)] = 1
membership[(microlist$IntertidalZone == "MiddleZone") & (microlist$p.value <= 0.05)] = 2
membership[(microlist$IntertidalZone == "HighZone") & (microlist$p.value <= 0.05)] = 3
membership[microlist$p.value > 0.05] = 4
modularity(g,membership)
cle <- cluster_leading_eigen(g,weights = abs(E(g)$weight))
#cle <- cluster_fast_greedy(g,weights = abs(E(g)$weight))
T <- data.frame(row.names = c(0,1,2,3))
for (i in 1:max(cle$membership)){
  for (j in c(0,1,2,3)){
    T[j+1,i] = sum(microlist[cle$names[cle$membership == i],]$GroupColor == j)/sum(microlist$GroupColor == j)
  }
}
T$id <- c('1. N','2. H','3. M','4. L')
df2 <- melt(T, id.vars='id')
ggplot(df2, aes(x=id, y=value, fill=variable)) + geom_bar(stat='identity', position='dodge')

macronet <- read.delim('MICTools_Network_Macro_PresenceAbsence/MIC_tools_output/strength_MACRO_PA.txt', header = TRUE, sep = "\t")
macronet <- macronet[,2:7]
macronet$weight <- macronet$PearsonR
g=graph.data.frame(macronet,directed = FALSE)
macrolist <- read.csv('2.Table_Indval_Sessile&Mobile_MACRO.csv', header = TRUE)
rownames(macrolist) <- macrolist[,1]
n <- V(g)$name
macrolist = macrolist[intersect(n,rownames(macrolist)),] # remove species not in network
membership <- rep(NA, 619)
membership[(macrolist$IntertidalZone == "LowZone") & (macrolist$p.value <= 0.05)] = 1
membership[(macrolist$IntertidalZone == "MiddleZone") & (macrolist$p.value <= 0.05)] = 2
membership[(macrolist$IntertidalZone == "HighZone") & (macrolist$p.value <= 0.05)] = 3
membership[macrolist$p.value > 0.05] = 4
modularity(g,membership)

cle <- cluster_leading_eigen(g,weights = abs(E(g)$weight))
#cle <- cluster_fast_greedy(g)
T <- data.frame(row.names = c(0,1,2,3))
for (i in 1:max(cle$membership)){
  for (j in c(0,1,2,3)){
    T[j+1,i] = sum(macrolist[cle$names[cle$membership == i],]$GroupColor == j)/sum(macrolist$GroupColor == j)
  }
}
T$id <- c('1. N','2. H','3. M','4. L')
df2 <- melt(T, id.vars='id')
ggplot(df2, aes(x=id, y=value, fill=variable)) + geom_bar(stat='identity', position='dodge')