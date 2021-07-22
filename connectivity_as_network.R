library(igraph)
library(RColorBrewer)
source('/Volumes/Macintosh HD/Users/mara_freilich/Dropbox (Personal)/10.MacroMicro_Networks_Chile/2.MACRO&MICRO/network_analysis_codes/links_between_heights.R')
source('/Volumes/Macintosh HD/Users/mara_freilich/Dropbox (Personal)/10.MacroMicro_Networks_Chile/2.MACRO&MICRO/network_analysis_codes/strength_between_heights.R')
rdbucols <- rev(brewer.pal(11, 'RdBu'))
newcol <- colorRampPalette(rdbucols)
rdbucols2 <- newcol(100)
coords <- matrix(data = 0,nrow = 4,ncol = 2)
coords[,1] <- c(1,1,2,2)
coords[,2] <- c(1,2,2,1)

net1<-graph_from_adjacency_matrix(smacro1,weighted = TRUE)
net2<-graph_from_adjacency_matrix(cmacro1,weighted = TRUE)
valcol <- (E(net1)$weight + 0.4)/(0.4*2)
plot(net1,vertex.color = c(brewer.pal(3,"YlOrBr"),"#636363"),edge.arrow.mode = 0,vertex.size = 40,
     vertex.label = c("low","mid","high","none"),edge.width = E(net2)$weight*15,
     edge.color = rdbucols2[valcol*100],layout = coords,edge.label = round(E(net2)$weight,2))

net1<-graph_from_adjacency_matrix(smicro1,weighted = TRUE)
net2<-graph_from_adjacency_matrix(cmicro1,weighted = TRUE)
valcol <- (E(net1)$weight + 0.4)/(0.4*2)
plot(net1,vertex.color = c(brewer.pal(3,"Purples"),"#636363"),edge.arrow.mode = 0,vertex.size = 40,
     vertex.label = c("low","mid","high","none"),edge.width = E(net2)$weight*15,
     edge.color = rdbucols2[valcol*100],layout = coords,edge.label = round(E(net2)$weight,2))