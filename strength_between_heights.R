library(igraph)
library(ggplot2)
library(reshape2)
micronet <- read.delim('MICTools_Network_Micro_PresenceAbsence/MIC_tools_output/strength_MICRO_PresenceAbsence_Network.txt', header = TRUE, sep = "\t")
micronet <- micronet[,2:7]
g=graph.data.frame(micronet)
g1=get.adjacency(g,attr = "PearsonR",sparse=FALSE) # adjacency matrix of microbes
g1=g1+t(g1) # symmetric matrix

microlist <- read.csv('2.Tabla_Indval_Final_Filter100_AllASVs_MICRO.csv', header = TRUE)
rownames(microlist) <- microlist[,1]
microlist = microlist[intersect(rownames(g1),rownames(microlist)),] # remove species not in network
L <- microlist[(microlist[,"IntertidalZone"] == "LowZone") & (microlist[,"Significant"] == 0),]
M <- microlist[(microlist[,"IntertidalZone"] == "MiddleZone") & (microlist[,"Significant"] == 0),]
H <- microlist[(microlist[,"IntertidalZone"] == "HighZone") & (microlist[,"Significant"] == 0),]
N <- microlist[microlist[,"Significant"] == 1,]

strength_micro <- matrix(data = 0,nrow = 4,ncol = 4)
# compute connectivity within each subgraph, matrices are symmetric
gt = g1[rownames(L),rownames(L)]
strength_micro[1,1] <- mean(gt)
gt = g1[rownames(L),rownames(M)]
strength_micro[1,2] <- mean(gt)
strength_micro[2,1] <- mean(gt)
gt = g1[rownames(L),rownames(H)]
strength_micro[1,3] <- mean(gt)
strength_micro[3,1] <- mean(gt)
gt = g1[rownames(L),rownames(N)]
strength_micro[1,4] <- mean(gt)
strength_micro[4,1] <- mean(gt)
gt = g1[rownames(M),rownames(M)]
strength_micro[2,2] <- mean(gt)
gt = g1[rownames(M),rownames(H)]
strength_micro[2,3] <- mean(gt)
strength_micro[3,2] <- mean(gt)
gt = g1[rownames(M),rownames(N)]
strength_micro[2,4] <- mean(gt)
strength_micro[4,2] <- mean(gt)
gt = g1[rownames(H),rownames(H)]
strength_micro[3,3] <- mean(gt)
gt = g1[rownames(H),rownames(N)]
strength_micro[3,4] <- mean(gt)
strength_micro[4,3] <- mean(gt)
gt = g1[rownames(N),rownames(N)]
strength_micro[4,4] <- mean(gt)
smicro1 <- strength_micro

macronet <- read.delim('MICTools_Network_Macro_PresenceAbsence/MIC_tools_output/strength_MACRO_PA.txt', header = TRUE, sep = "\t")
macronet <- macronet[,2:7]
g=graph.data.frame(macronet)
g2=get.adjacency(g,attr = "PearsonR",sparse=FALSE) # adjacency matrix of microbes
g2=g2+t(g2) # symmetric matrix

macrolist <- read.csv('2.Table_Indval_Sessile&Mobile_MACRO.csv', header = TRUE)
rownames(macrolist) <- macrolist[,1]
macrolist = macrolist[intersect(rownames(g2),rownames(macrolist)),] # remove species not in network
L <- macrolist[(macrolist[,"IntertidalZone"] == "Low") & (macrolist[,"pvalue"] < 0.05),]
M <- macrolist[(macrolist[,"IntertidalZone"] == "Middle") & (macrolist[,"pvalue"] < 0.05),]
H <- macrolist[(macrolist[,"IntertidalZone"] == "High") & (macrolist[,"pvalue"] < 0.05),]
N <- macrolist[macrolist[,"pvalue"] < 0.05,]

strength_macro <- matrix(data = 0,nrow = 4,ncol = 4)
# compute connectivity within each subgraph, matrices are symmetric
gt = g2[rownames(L),rownames(L)]
strength_macro[1,1] <- mean(gt)
gt = g2[rownames(L),rownames(M)]
strength_macro[1,2] <- mean(gt)
strength_macro[2,1] <- mean(gt)
gt = g2[rownames(L),rownames(H)]
strength_macro[1,3] <- mean(gt)
strength_macro[3,1] <- mean(gt)
gt = g2[rownames(L),rownames(N)]
strength_macro[1,4] <- mean(gt)
strength_macro[4,1] <- mean(gt)
gt = g2[rownames(M),rownames(M)]
strength_macro[2,2] <- mean(gt)
gt = g2[rownames(M),rownames(H)]
strength_macro[2,3] <- mean(gt)
strength_macro[3,2] <- mean(gt)
gt = g2[rownames(M),rownames(N)]
strength_macro[2,4] <- mean(gt)
strength_macro[4,2] <- mean(gt)
gt = g2[rownames(H),rownames(H)]
strength_macro[3,3] <- mean(gt)
gt = g2[rownames(H),rownames(N)]
strength_macro[3,4] <- mean(gt)
strength_macro[4,3] <- mean(gt)
gt = g2[rownames(N),rownames(N)]
strength_macro[4,4] <- mean(gt)
smacro1 <- strength_macro

strength_micro = data.frame(round(strength_micro,digits = 2))
colnames(strength_micro) <- c('1. L','2. M','3. H','4. N')
strength_micro$row <- c('1. L','2. M','3. H','4. N')
micro_melt <- melt(data = strength_micro, id.vars = "row")
d <- ggplot(data=micro_melt,
       aes(x=variable, y=row, fill=value)) + geom_tile() + 
  geom_text(aes(label=value), color='white') + theme_bw() 
d + scale_fill_distiller(palette = "RdBu",limits = c(-0.4,0.4))

strength_macro = data.frame(round(strength_macro,digits = 2))
colnames(strength_macro) <- c('1. L','2. M','3. H','4. N')
strength_macro$row <- c('1. L','2. M','3. H','4. N')
macro_melt <- melt(data = strength_macro, id.vars = "row")
d <- ggplot(data=macro_melt,
       aes(x=variable, y=row, fill=value)) + geom_tile() + 
  geom_text(aes(label=value), color='white') + theme_bw()
d + scale_fill_distiller(palette = "RdBu",limits = c(-0.27,0.27))