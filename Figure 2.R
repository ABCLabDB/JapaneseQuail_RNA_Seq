#### Result 3 code for validation ####

library(ggplot2)
library(ggsci)
library(ggpubr)
library(data.table)
library(edgeR)
library(DESeq2)
library(ggrepel)
library(dplyr)
library(RColorBrewer)
library(GenomicRanges)
library(data.table)
library(Gviz)
library(tidyverse)

#### Test for finding sex-biased genes in brain tissue ####

meta_data <- data.frame(fread("./Data/Figure2_data/meta_data_FB_MB.txt", sep = "\t",header=T))

fileList <- paste0(meta_data$Run, '.txt')

temp <- data.frame(fread(paste0("./Data/Quantification_quail/",fileList[1])))



geneSymbol <- temp$Geneid
count <- data.frame(temp[, 7])


for (i in 2:length(fileList)) {
  temp <- data.frame(fread(paste0("./Data/Quantification_quail/",fileList[i])))
  count <- data.frame(count, temp[, 7])
}




## Loading gene annotation
temp <- data.frame(geneSymbol)
Annotation <- data.frame(fread("./Data/mart_export.txt", sep = "\t", head = T))
temp <-  merge(
  x = temp,
  y = Annotation,
  by.x = "geneSymbol",
  by.y = "Gene.stable.ID",
  sort = F
)
temp
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



## Defining experimental variables


Tissue <- factor(meta_data$Tissue, levels = c("Female brain", "Male brain"))



design <- model.matrix( ~ Tissue)



## TMM normalization using edgeR package
Y <- DGEList(counts = count, genes = geneSymbol)


Y <- calcNormFactors(Y, method = "TMM") ## TMM Normalizatio

Y <- estimateDisp(Y, design)



fit <- glmFit(Y, design)## GLM


logCPM <- cpm(Y, normalized.lib.sizes = TRUE, log = T)



Result_Brain <- glmLRT(fit, coef = 2) 

Result_Brain_table <- topTags(Result_Brain,n = dim(logCPM)[1], sort.by = "none")$table



#### Figure 3A ####
data <- fread("./Data/Figure2_data/Result_FB_MB.txt")



# Calculating -log10 of P-value for better visualization
data$negLogPval <- -log10(data$PValue)


data$color <- with(data, ifelse(logFC> 0 & FDR < 0.05, "#43A1F7",
                                ifelse(logFC < 0 & FDR < 0.05, "#F74C43", "grey")))



# Creating the volcano plot
ggplot(data, aes(x = logFC, y = negLogPval)) +
  geom_jitter(aes(color=color),alpha = 0.7,size=5) +  # Highlight significant points
  scale_color_identity() +  
  theme_classic()+
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "black",size=1) +  # Highlight fold change cutoff
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black",size=1) +
  ggrepel::geom_text_repel(aes(label = ifelse(color %in% c("#43A1F7", "#F74C43"), as.character(Gene.name), "")),color="black",fontface="italic")+
  theme(axis.text.x = element_text(size = 10,color="black"), 
        axis.text.y = element_text(size = 10,color="black"),
        axis.title.x=element_blank(),
        axis.title.y=element_blank()
  )  



ggsave(filename="Volcano.pdf",width = 5,height = 4)



#### Figure 3B ####

data <- fread("./Data/Figure2_data/Result_FB_MB.txt")
data <- data[data$FDR<=0.05,]

temp <- data

temp <- temp %>%
  mutate(Pattern=ifelse(logFC>=0,"Male biased gene","Female biased gene"))

temp_table <- table(temp$Pattern)

gender_df <- as.data.frame(temp_table)

gender_df$percentage <- with(gender_df, round(Freq / sum(Freq) * 100, 1))


temp <- temp %>%
  mutate(chromosome_type=ifelse(Chromosome.scaffold.name=="Z","Z chromosome","Autosome"))

temp_table <- table(temp$Chromosome.scaffold.name)

chromosome_df <- as.data.frame(temp_table)


chromosome_df$Var1 <- factor(chromosome_df$Var1,levels=c("1","2","3","5","9","26","Z"))

chromosome_df$percentage <- with(chromosome_df, round(Freq / sum(Freq) * 100, 1))


## Color mapping ##
gender_colors <- rev(pal_lancet()(2))
names(gender_colors) <- gender_df$Var1

chromosome_colors <- pal_simpsons()(7) 
names(chromosome_colors) <- chromosome_df$Var1

# Drawing piechart
ggplot() +
  geom_bar(data=chromosome_df, aes(x = 2, y = percentage, fill = Var1), stat = "identity", width = 1,color="black",size=1) +
  geom_bar(data=gender_df, aes(x = 1.5, y = percentage, fill = Var1), stat = "identity", width = 1,color="black",size=1) +
  coord_polar(theta = "y") +
  scale_fill_manual(values = c(gender_colors, chromosome_colors)) +
  theme_void() +
  theme(legend.position = "None")


ggsave("Piechart.pdf",width=3,height=3)

#### Figure 3C-3D, and 3F ####
load(file="./Data/Figure2_data/FB_MB.RData")

temp <- Result_Brain_table[which(Result_Brain_table$FDR<0.05),]


idx <- c()



idx <- c(idx,which(Annotation$geneSymbol=="ENSCJPG00005008895" | Annotation$geneSymbol=="ENSCJPG00005012300" | Annotation$geneSymbo=="ENSCJPG00005003198")
)


temp_name <- paste0(Annotation$Gene.name[idx])


for (i in 1:length(temp_name)) {
  if (temp_name[i] == "") {
    temp_name[i] <- geneSymbol[idx][i]
  }
}




# Chromosome Information Add

plot_data <- data.frame(
  Expression = as.numeric(logCPM[idx[1], ]),
  GeneSymbol = rep(temp_name[1], ncol(logCPM)),
  Tissue
)


for (i in 2:length(temp_name)) {
  temp <- data.frame(
    Expression = as.numeric(logCPM[idx[i], ]),
    GeneSymbol = rep(temp_name[i], ncol(logCPM)),
    Tissue
  )
  
  
  plot_data <- rbind(plot_data, temp)
}

plot_data$Tissue <- factor(plot_data$Tissue)


genes_of_interest <- c("ENSCJPG00005008895","ENSCJPG00005012300","ENSCJPG00005003198")



for(i in 1:length(genes_of_interest)){
  
  filtered_data <- plot_data[plot_data$GeneSymbol %in% genes_of_interest[i], ]
  

  p2 <- ggplot(filtered_data, aes(x = Tissue, y = Expression, fill = Tissue)) +
    geom_boxplot(outlier.shape = NA, alpha = 0.8) +  # 박스플롯 설정
    geom_jitter(aes(color = Tissue),alpha = 0.6, size = 5) +  # 지터 플롯 추가
    scale_fill_manual(values = c("Female brain" = "#ed6a66", "Male brain" = "#00239C")) +  
    scale_color_manual(values = c("Female brain" = "#ed6a66", "Male brain" = "#00239C")) +  
    theme_classic2()+
    theme(
      legend.position = "none",
      axis.title.x = element_blank(),
      axis.title.y = element_blank(),
      axis.text = element_text(size=12, color="black"),
      strip.text = element_text(size=12, color="black")
    )
  
  p2
  
  ggsave(paste0(genes_of_interest[i],".pdf"),width=3,height=3)
  
}

#### Figure 3E ####

options(scipen=30)
color_label <- function(data){
  label <- c()
  
  color <- c()
  #ABF200
  for(i in 1:nrow(data)){
    
    if(abs(data$logFC[i]) == 0){
      label <- c(label,"FC=0")
    }else if(abs(data$logFC[i]) > 0  && abs(data$logFC[i]) <= 0.5){
      label <- c(label,"0<FC<0.5")
    } else if(abs(data$logFC[i]) > 0.5  && abs(data$logFC[i]) <= 1){
      label <- c(label,"0.5<FC<1")
    } else if(abs(data$logFC[i]) > 1 && abs(data$logFC[i]) <= 2 ){
      label <- c(label,"1<FC<2")
    }else if(abs(data$logFC[i]) > 2 && abs(data$logFC[i]) <= 3){
      label <- c(label,"2<FC<3")
    }else if(abs(data$logFC[i]) > 3 && abs(data$logFC[i]) <= 4){
      label <- c(label,"3<FC<4")
    }else if(abs(data$logFC[i]) > 4 && abs(data$logFC[i]) <= 5){
      label <- c(label,"4<FC<5")
    }else if(abs(data$logFC[i]) > 5 && abs(data$logFC[i]) <= 6){
      label <- c(label,"5<FC<6")
    }else if(abs(data$logFC[i]) > 6){
      label <- c(label,"FC>6")
    }
    
  }
  
  return(label)
}

data <- fread("./Data/Figure2_data/Result_FB_MB.txt")

Chromosome_index <- which(data$Chromosome.scaffold.name=="Z")


data <- data[Chromosome_index]

label_list <- color_label(data)

data$label <- factor(label_list,levels = c("0<FC<0.5","0.5<FC<1","1<FC<2","2<FC<3","3<FC<4","4<FC<5","FC>6"))
# Reordering #
data <- data %>% arrange(desc(label))

# Replace gene name with geneSymbol #
data$Gene.name[which(data$Gene.name == "")] <- data$geneSymbol[which(data$Gene.name == "")]


num_colors <- length(unique(data$label))
colors <- rev(brewer.pal(num_colors, "RdYlGn"))


Dosage_brain <- ggplot(data=data,aes(x=Gene.start..bp.,y=logFC,label=Gene.name,color=label))+
  geom_jitter(alpha=0.8,size=5)+
  scale_color_manual(values=colors)+
  coord_cartesian(ylim = c(-7, 5))+
  geom_hline(yintercept = 0.512,size=1,linetype='longdash',color="red")+
  geom_hline(yintercept = 0.839,size=1,linetype='longdash',color="blue")+
  geom_text_repel(fontface="italic",color="black")+
  theme_bw()+
  theme_classic()+
  theme(
    axis.title = element_blank(),
    axis.text = element_text(size=12,color="black"),
    strip.text=element_text(size=12,color="black"),
    legend.position = "None"
  )

Dosage_brain

ggsave("ScatterPlot.pdf",width=12,height=4)






#### Figure 3G ####

# Load annotation data of chicken #
data <- fread("./Data/Figure2_data/Galgal6_annotation.txt")
data <- data[data$`Chromosome/scaffold name` == "Z", ]

# Gene list #
gene_list <- c("SLC1A1", "CDC37L1", "AK3", "RCL1")

# Extraction of genes #
idx <- which(data$`Gene stable ID` == "ENSGALG00000051419" | data$`Gene name` %in% gene_list)

# get gene start base pair #
temp <- data[idx, ]
temp <- temp[order(temp$`Gene start (bp)`), ]


gene_data <- data.frame(
  gene_name = temp$`Gene name`,
  start = temp$`Gene start (bp)`,
  end = temp$`Gene end (bp)`,
  strand = temp$Strand
)

gene_data$gene_name[2] <- "ENSGALG00000051419"


genes_on_z <- GRanges(
  seqnames = rep("chrZ", nrow(gene_data)),
  ranges = IRanges(start = gene_data$start, end = gene_data$end),
  strand = gene_data$strand
)

mcols(genes_on_z)$gene_name <- gene_data$gene_name

genomeAxis <- GenomeAxisTrack()
geneTrack <- GeneRegionTrack(
  range = genes_on_z,
  name = "Genes on Z",
  symbol = as.factor(gene_data$gene_name),
  showId = TRUE,
  geneSymbol=TRUE,
  shape="arrow",
  height=0.2,
  lwd=0.1
)


plotTracks(list(genomeAxis, geneTrack), from = 27000000, to = 27520000)


# Load the gene annotation of Japanse quail #


data <- fread("./Data/Figure2_data/quail_annotation.txt")
data <- data[data$Chromosome.scaffold.name == "Z", ]


gene_list <- c("SLC1A1", "CDC37L1", "AK3", "RCL1")


idx <- which(data$geneSymbol == "ENSCJPG00005003198" | data$Gene.name %in% gene_list)

temp <- data[idx, ]
temp <- temp[order(temp$Gene.start..bp.), ]


gene_data <- data.frame(
  gene_name = temp$Gene.name,
  start = temp$Gene.start..bp.,
  end = temp$Gene.end..bp.,
  strand = temp$Strand
)

gene_data$gene_name[2] <- "ENSCJPG00005003198"


genes_on_z <- GRanges(
  seqnames = rep("chrZ", nrow(gene_data)),
  ranges = IRanges(start = gene_data$start, end = gene_data$end),
  strand = gene_data$strand
)


mcols(genes_on_z)$gene_name <- gene_data$gene_name


genomeAxis <- GenomeAxisTrack()
geneTrack <- GeneRegionTrack(
  range = genes_on_z,
  name = "Genes on Z",
  symbol = as.factor(gene_data$gene_name),
  showId = TRUE,
  geneSymbol=TRUE,
  shape="arrow",
  height=0.2,
  lwd=0.1
)

plotTracks(list(genomeAxis, geneTrack), from = 24500000, to = 24760000)








#### Figure 3H ####

data <- fread("./Data/Figure2_data/David_result_brain.txt")

palette <- c("#537c78","#7ba591","#CC222B","#F15B4C","#FAA41b","#FFD45B")

palette <- c("#ff595e","#ff924c", "#ffca3a",  "#c5ca30",  "#8ac926",  "#52a675",  "#1982c4",  "#4267ac",  "#6a4c93","#b5a6c9")

setorder(data, PValue)


data$Term <- factor(data$Term, levels = unique(data$Term))

Zplot <- ggplot(data, aes(x=Term, y=-log10(data$PValue), fill=Term)) +
  geom_bar(stat="identity", colour="black", size=1,width=0.7) + 
  scale_fill_manual(values = palette) + 
  theme_classic() +
  scale_y_continuous(expand = c(0,0)) +
  geom_hline(yintercept=-log10(0.01),linetype=2,size=1.3) +
  theme(legend.position = "none",
        axis.title.y=element_blank(),
        axis.text = element_text(size = 12),
        panel.grid.major.y = element_line(color = "#d4d0cf", linetype = "dashed"),
        panel.grid.minor.y = element_line(color = "#d4d0cf", linetype = "dashed"),
        axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 15, hjust=1, color="black"),
        axis.text.y = element_text(size=15, color="black")
  )

Zplot

ggsave("GOterm.pdf",width=10,height=5)





#### Figure 3I ####

load("Data/Figure2_data/FB_MB.RData")


idx <- c()


data <- fread("Data/Figure2_data/WD.txt")

idx <- which(Annotation$geneSymbol %in% data$geneSymbol)


temp <- Annotation[idx,]


temp_name <- Annotation$Gene.name[idx]
for (i in 1:length(temp_name)) {
  if (temp_name[i] == "") {
    temp_name[i] <- geneSymbol[idx][i]
  }
}


# Chromosome Information Add

plot_data <- data.frame(
  Expression = as.numeric(logCPM[idx[1], ]),
  GeneSymbol = rep(temp_name[1], ncol(logCPM)),
  Tissue
)


for (i in 2:length(temp_name)) {
  temp <- data.frame(
    Expression = as.numeric(logCPM[idx[i], ]),
    GeneSymbol = rep(temp_name[i], ncol(logCPM)),
    Tissue
  )
  
  
  plot_data <- rbind(plot_data, temp)
}



long_data <- plot_data %>%
  group_by(GeneSymbol, Tissue) %>%
  mutate(Observation = paste(Tissue, row_number(), sep="_")) %>%
  ungroup()


wide_data <- long_data %>%
  select(-Tissue) %>%
  pivot_wider(names_from = Observation, values_from = Expression)



matrix_data <- wide_data %>% select(-GeneSymbol)



rownames(matrix_data) <- wide_data$GeneSymbol

name <- wide_data$GeneSymbol

colname <- c(rep("Female brain",3),rep("Male brain",3))


colMeans(matrix_data)


selected_columns <- matrix_data[, c(1, 3, 5)]

averages <- colMeans(selected_columns, na.rm = TRUE)
mean(averages)

selected_columns <- matrix_data[, c(2, 4, 6)]
averages <- colMeans(selected_columns, na.rm = TRUE)
mean(averages)


pheatmap::pheatmap(
  matrix_data,
  border_color = "black",
  angle_col = 0,

  color=colorRampPalette(c("#4575b4", "#ffffbf", "#d73027"))(10),
  cellwidth = 20,  # 셀의 폭 조절
  cellheight = 10  # 셀의 높이 조절
)



