#### Result 2 code for validation ####

library(ggplot2)
library(ggsci)
library(ggpubr)
library(data.table)
library(edgeR)
library(DESeq2)
library(ggrepel)
library(dplyr)

#### Figure 2A ####

#### Mapping Rate Figure
mapping_data <- data.frame(fread("./Data/mapping_rate_tissue.txt", sep = "\t",header=T))


de_novo_data <- mapping_data[which(mapping_data$Mapping_way=="De novo"),]
reference_data <- mapping_data[which(mapping_data$Mapping_way=="Reference"),]


colnames(mapping_data)[1] <- "Approach"
mapping_data$Approach <- factor(mapping_data$Approach,levels = c("Reference","De novo"))


figure2 <- ggplot(mapping_data, aes(x = Approach, y = Rate,color = Approach, fill=Approach)) +
  geom_jitter(size = 7, alpha = 0.8, position = position_jitter(seed = 1)) +
  geom_boxplot(outlier.shape = NA, alpha = 0, size = 0.5) +
  theme_classic() +
  scale_color_manual(values = c("#f4b41a","#194a6b")) +
  scale_fill_manual(values = c("#f4b41a","#194a6b")) +
  ylab("Alignment rate (%)") +
  theme(axis.title = element_text(size =14), axis.text = element_text(size = 14), 
        axis.title.x = element_blank(), legend.position = "none") 

figure2

ggsave("Figure 2A.pdf",
       plot = figure2, 
       height = 4 , width = 4)




#### Figure 2B ####
meta_data <- data.frame(fread("./Data/meta_data.txt", sep = "\t",header=T))


fileList <- paste0(meta_data$Run, '.txt')

temp <- data.frame(fread(paste0("./Data/Quantification_quail/",fileList[1])))



geneSymbol <- temp$Geneid
count <- data.frame(temp[, 7])


for (i in 2:length(fileList)) {
  temp <- data.frame(fread(paste0("./Data/Quantification_quail/",fileList[i])))
  count <- data.frame(count, temp[, 7])
}



# Loading gene annotation
Annotation <- data.frame(fread("./Data/mart_export.txt", sep = "\t", head = T))
temp <- data.frame(geneSymbol)


temp <-  merge(
  x = temp,
  y = Annotation,
  by.x = "geneSymbol",
  by.y = "Gene.stable.ID",
  sort = F
)

Annotation <- temp

all.equal(as.character(Annotation$geneSymbol), geneSymbol)

data.frame(as.character(Annotation$geneSymbol), geneSymbol)



## Removing unnecessary chromosome

indexRemoval <- c(grep("KQ",Annotation$Chromosome.scaffold.name),grep("LSZ",Annotation$Chromosome.scaffold.name),grep("LGE",Annotation$Chromosome.scaffold.name),grep("MT",Annotation$Chromosome.scaffold.name))



## Removing all zero counted genes

indexRemoval <- c(indexRemoval,(which(rowSums(count) == 0)))

count <- count[-indexRemoval, ]
Annotation <- Annotation[-indexRemoval, ]
geneSymbol <- geneSymbol[-indexRemoval]


colnames(count) <- meta_data$Run

row.names(count) <- Annotation$geneSymbol


row.names(count)



## Defining experimental variables
Tissue <- factor(meta_data$Tissue, levels = c("Gonad", "Brain"))
Sex <- factor(meta_data$Sexual, levels = c("Female", "Male"))



## TMM normalization using edgeR package
Y <- DGEList(counts = count, genes = geneSymbol)
Y <- calcNormFactors(Y, method = "TMM") ## TMM Normalization



## MDS plotting
MDS_data <- plotMDS(Y)

plot_data <- data.frame(Sex:Tissue, Tissue, Sex, X = MDS_data$x, Y = MDS_data$y)
length(plot_data$Tissue)

plot_data$Tissue <- c("Female brain","Male brain","Ovary","Testis","Female brain","Male brain","Ovary","Testis","Female brain","Male brain","Ovary","Testis")

## Whole MDS plot

plot_data$Tissue <- factor(plot_data$Tissue,levels = c("Ovary","Testis","Female brain","Male brain"))


MDSplot <- ggplot(plot_data, aes(x = X, y = Y ,color = Tissue)) +
  geom_jitter(size = 5, alpha = 0.8) +
  theme_classic() +
  scale_color_manual(values = c("#fecb3e","#fc8370","#c2549d","#7e549e")) +
  xlab("MDS1") +
  ylab("MDS2") +
  coord_cartesian(xlim = c(min(plot_data$X)-0.3,max(plot_data$X)+0.3),
                  ylim = c(min(plot_data$Y)-0.3,max(plot_data$Y)+0.3)) +
  theme(axis.title = element_text(size =14), axis.text = element_text(size = 14), legend.title = element_blank(),legend.position = "None") 

MDSplot

ggsave("Figure 1B.pdf",
       plot = MDSplot,
       width=4,height=4)







#### Figure 2C ####
data <- fread("./Data/RefvsDenovo.txt")


col <- c()
gglabel <- c()
legend <- c()



for(i in 1:nrow(data)){
  if(data$adjP_brain[i]<=0.05 && data$FB_MB[i]<=0.05){
    col <- c(col,"#FF0000")
    gglabel <- c(gglabel,as.character(data$GeneName[i]))
    legend <- c(legend,"Common")
    next
  }
  
  if(data$adjP_brain[i]<=0.05){
    col <- c(col,"#E4F7BA")
    gglabel <- c(gglabel,"")
    legend <- c(legend,"Reference")
    next
  }
  
  if(data$FB_MB[i]<=0.05){
    col <- c(col,"#D9E5FF")
    gglabel <- c(gglabel,"")
    legend <- c(legend,"Denovo")
    next
    
  }
  
  col <- c(col,"#EAEAEA")
  gglabel <- c(gglabel,"")
  legend <- c(legend,"Non")
  
}


data$color <- col
data$legend <- legend



cor(data$FB_MB,data$adjP_brain)

data <- data %>% dplyr::mutate(label = ifelse(FB_MB<=0.05|adjP_brain<=0.05, Gene.name,NA)) 

figure1C <- ggplot(data, aes(x=log10(adjP_brain)*-1,y=log10(FB_MB)*-1)) +
  geom_jitter(size = 5, alpha = 0.7, aes(color = legend)) +
  geom_smooth(method = "lm", color = "red",size=0.7) +
  geom_text_repel(aes(label=label),fontface="italic") +
  theme_classic() +
  scale_color_manual(values = c("#62bd70","#194a6b","grey60","#f4b41a"),
                     labels = c("Commonly significant genes at FDR adjusted P<=0.05","Significant genes only in de novo approach",
                                "Genes with no significance","Significant genes only in reference genome based approach")) +
  xlab("-log10 (FDR adjusted P-value) of the test for difference between sexes\n in reference genome-based approach.") +
  ylab("-log10 (FDR adjusted P-value) of the test for difference between sexes\n in de novo approach.") +
  theme(axis.title = element_text(size =12), axis.text = element_text(size = 14), legend.title = element_blank(),
        legend.text = element_text(size = 12),legend.position = "None"       ) 

figure1C

ggsave("Figure1C.pdf",
       plot = figure1C,
       width=6,height=3)

#### Figure 2D ####

data <- fread("Data/RefvsDenovo.txt")
col <- c()
gglabel <- c()
legend <- c()


for(i in 1:nrow(data)){
  if(data$adjP_gonad[i]<=0.05 && data$O_T[i]<=0.05){
    col <- c(col,"#FF0000")
    gglabel <- c(gglabel,as.character(data$GeneName[i]))
    legend <- c(legend,"Common")
    next
  }
  
  if(data$adjP_gonad[i]<=0.05){
    col <- c(col,"#E4F7BA")
    gglabel <- c(gglabel,"")
    legend <- c(legend,"Reference")
    next
  }
  
  if(data$O_T[i]<=0.05){
    col <- c(col,"#D9E5FF")
    gglabel <- c(gglabel,"")
    legend <- c(legend,"Denovo")
    next
    
  }
  
  col <- c(col,"#EAEAEA")
  gglabel <- c(gglabel,"")
  legend <- c(legend,"Non")
  
}

data$color <- col
data$legend <- legend



data$legend <- factor(data$legend,levels = c("Common","Denovo","Reference","Non"))


data <- data %>% 
  arrange(match(legend, levels(legend)))



data <- data %>% dplyr::mutate(label = ifelse(O_T<=0.05|adjP_gonad<=0.05, Gene.name,NA)) 



figure1D <- ggplot(data, aes(x=log10(adjP_gonad)*-1,y=log10(O_T)*-1)) +
  geom_jitter(size = 5, alpha = 0.8, aes(color = legend)) +
  geom_smooth(method = "lm", color = "red",size=0.7) +
  geom_text_repel(aes(label=label),fontface="italic") +
  theme_classic() +
  scale_color_manual(values = c("#62bd70","#194a6b","#f4b41a","grey60"),
                     labels = c("Common DEG in both method","De novo transcript assembly specific DEG",
                                "Reference guided specific DEG","Non significant gene in both method")) +
  xlab("-log10(FDR adjusted p-value)\nin Reference guided method") +
  ylab("-log10(FDR adjusted p-value)\nin De novo transcript assembly method") +
  theme(axis.title = element_text(size =12), axis.text = element_text(size = 14), legend.title = element_blank(),
        legend.position = "none") 

figure1D

ggsave("Figure1D.pdf",
       plot = figure1D,
       width=6,height=3)



#### Figure 2E ####
min_max_normalize <- function(x) {
  return((x - min(x)) / (max(x) - min(x)))
}



gene_expression <- fread("Data/GSE64961_Processed_data.csv")






## Loading featureCount data
meta_data <- data.frame(fread("Data/meta_data.txt", sep = "\t",header=T))


fileList <- paste0(meta_data$Run, '.txt')

temp <- data.frame(fread(paste0("./Data/Quantification_quail/",fileList[1])))


geneSymbol <- temp$Geneid
count <- data.frame(temp[, 7])



for (i in 2:length(fileList)) {
  temp <- data.frame(fread(paste0("./Data/Quantification_quail/",fileList[i])))
  count <- data.frame(count, temp[, 7])
}



## Loading gene annotation
Annotation <- data.frame(fread("Data/mart_export.txt", sep = "\t", head = T))
temp <- data.frame(geneSymbol)

temp <-  merge(
  x = temp,
  y = Annotation,
  by.x = "geneSymbol",
  by.y = "Gene.stable.ID",
  sort = F
)

Annotation <- temp

gene_list <- "VIP"

gene_expression <- gene_expression[gene_expression$Gene_Symbol %in% gene_list,]

## Get raw count data of VIP gene ##
idx <- which(Annotation$Gene.name %in% gene_list)

count <- count[idx, ]
Annotation <- Annotation[idx, ]
geneSymbol <- geneSymbol[idx]



method <- c(rep("Denovo",12),rep("Reference",12))


Tissue <- c("Female brain","Male brain","Ovary","Testis","Female brain","Male brain","Ovary","Testis","Female brain","Male brain","Ovary","Testis")


Tissue <- rep(Tissue,2)


count_matrix <- data.frame(gene_expression[,-1],count)



normalized_df <- as.data.frame(t(apply(count_matrix, 1, min_max_normalize)))


temp <- t(normalized_df)


df <- data.frame(Method=method,expression=temp[,1],Tissue=Tissue)

df$Method <- factor(df$Method,levels = c("Reference","Denovo"))


df$SampleID <- as.factor(c(meta_data$Run,meta_data$Run))

df$Tissue <- as.factor(df$Tissue)
color <- c("#A4193D","#FFDFB9")
# 바 플롯 그리기
figure1E <- ggplot(df, aes(x=Tissue, y=expression, fill=Method)) +
  geom_bar(stat="identity", position=position_dodge()) +
  scale_fill_manual(values=color) +
  scale_y_continuous(expand=c(0,0))+
  labs(x="Tissue", y="Min-max normalzied\ngene expression level") +
  theme_classic() +
  theme(legend.title = element_blank(),
        panel.grid.major.x = element_line(color = "#9E9B9A", linetype = "dashed"),
        panel.grid.minor.x = element_line(color = "#9E9B9A", linetype = "dashed",size=0.5),
        axis.title.x = element_blank(), axis.title.y = element_text(size = 12), axis.text.y = element_text(size = 12),
        legend.text = element_text(size = 12), 
        axis.text.x = element_text(size = 12, angle = 20, vjust = 0.7),
        legend.position = "None"
  )

figure1E
ggsave("Figure1E.pdf",
       plot = figure1E,
       width=4,height=3)




#### Figure 2F ####

min_max_normalize <- function(x) {
  return((x - min(x)) / (max(x) - min(x)))
}

gene_expression <- fread("Data/GSE64961_Processed_data.csv")





## Loading featureCount data
meta_data <- data.frame(fread("Data/meta_data.txt", sep = "\t",header=T))


fileList <- paste0(meta_data$Run, '.txt')

temp <- data.frame(fread(paste0("./Data/Quantification_quail/",fileList[1])))


geneSymbol <- temp$Geneid
count <- data.frame(temp[, 7])



for (i in 2:length(fileList)) {
  temp <- data.frame(fread(paste0("./Data/Quantification_quail/",fileList[i])))
  count <- data.frame(count, temp[, 7])
}



## Loading gene annotation
Annotation <- data.frame(fread("Data/mart_export.txt", sep = "\t", head = T))
temp <- data.frame(geneSymbol)

temp <-  merge(
  x = temp,
  y = Annotation,
  by.x = "geneSymbol",
  by.y = "Gene.stable.ID",
  sort = F
)

Annotation <- temp

gene_list <- "CCNH"

gene_expression <- gene_expression[gene_expression$Gene_Symbol %in% gene_list,]

## Get raw count data of CCNH gene ##
idx <- which(Annotation$Gene.name %in% gene_list)

count <- count[idx, ]
Annotation <- Annotation[idx, ]
geneSymbol <- geneSymbol[idx]



method <- c(rep("Denovo",12),rep("Reference",12))


Tissue <- c("Female brain","Male brain","Ovary","Testis","Female brain","Male brain","Ovary","Testis","Female brain","Male brain","Ovary","Testis")


Tissue <- rep(Tissue,2)


count_matrix <- data.frame(gene_expression[,-1],count)



normalized_df <- as.data.frame(t(apply(count_matrix, 1, min_max_normalize)))


temp <- t(normalized_df)


df <- data.frame(Method=method,expression=temp[,1],Tissue=Tissue)

df$Method <- factor(df$Method,levels = c("Reference","Denovo"))


df$SampleID <- as.factor(c(meta_data$Run,meta_data$Run))

df$Tissue <- as.factor(df$Tissue)
color <- c("#A4193D","#FFDFB9")
# 바 플롯 그리기
figure1F <- ggplot(df, aes(x=Tissue, y=expression, fill=Method)) +
  geom_bar(stat="identity", position=position_dodge()) +
  scale_fill_manual(values=color) +
  scale_y_continuous(expand=c(0,0))+
  labs(x="Tissue", y="Min-max normalzied\ngene expression level") +
  theme_classic() +
  theme(legend.title = element_blank(),
        panel.grid.major.x = element_line(color = "#9E9B9A", linetype = "dashed"),
        panel.grid.minor.x = element_line(color = "#9E9B9A", linetype = "dashed",size=0.5),
        axis.title.x = element_blank(), axis.title.y = element_text(size = 12), axis.text.y = element_text(size = 12),
        legend.text = element_text(size = 12), 
        axis.text.x = element_text(size = 12, angle = 20, vjust = 0.7),
        legend.position = "None"
  )

figure1F
ggsave("Figure1F.pdf",
       plot = figure1F,
       width=4,height=3)


