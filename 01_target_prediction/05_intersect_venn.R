library(ggvenn)
compoundName="DEP"
diseaseName="Brain injury"
setwd(".")
geneList=list()
rt=read.table("Compound.txt", header=F, sep="\t", check.names=F)
geneNames=as.vector(rt[,1])
geneList[[compoundName]]=geneNames
rt=read.table("Disease.txt", header=F, sep="\t", check.names=F)
geneNames=as.vector(rt[,1])
geneList[[diseaseName]]=geneNames
pdf(file="venn.pdf", width=6, height=6)
ggvenn(geneList,show_percentage = T,
	stroke_color = "white", stroke_size = 0.5,
	fill_color = c("#E41A1C","#1E90FF"),
	set_name_color =c("#E41A1C","#1E90FF"),
	set_name_size=6, text_size=4.5)
dev.off()
interGenes=Reduce(intersect, geneList)
write.table(file="interGenes.txt", interGenes, sep="\t", quote=F, col.names=F, row.names=F)
networkTab=rbind(cbind(compoundName,interGenes, "Compound"), cbind(diseaseName,interGenes, "Disease"))
colnames(networkTab)=c("Node1", "Node2", "Type")
write.table(file="net.network.txt", networkTab, sep="\t", quote=F, row.names=F)
nodeTab=rbind(cbind(compoundName,"Compound"), cbind(diseaseName,"Disease"), cbind(interGenes,"Gene"))
colnames(nodeTab)=c("Node", "Type")
write.table(file="net.node.txt", nodeTab, sep="\t", quote=F, row.names=F)
