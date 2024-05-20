#### Result 4 code for validation ####

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
library(viridis)

#### Test for finding sex-biased genes in reproductive organ ####

meta_data <- data.frame(fread("./Data/meta_data_O_T.txt", sep = "\t",header=T))

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

dim(count)
dim(Annotation)
length(geneSymbol)


colSums(count)



## Defining experimental variables


Tissue <- factor(meta_data$Tissue, levels = c("Ovary", "Testis"))
design <- model.matrix( ~ Tissue)



## TMM normalization using edgeR package
Y <- DGEList(counts = count, genes = geneSymbol)
Y <- calcNormFactors(Y, method = "TMM") ## TMM Normalizatio
Y <- estimateDisp(Y, design)
fit <- glmFit(Y, design)## GLM
logCPM <- cpm(Y, normalized.lib.sizes = TRUE, log = T)
Result_gonad <- glmLRT(fit, coef = 2) 
Result_gonad_table <- topTags(Result_gonad,n = dim(logCPM)[1], sort.by = "none")$table





#### Figure 4A ####


results <- fread("./Data/Result_O_T.txt", head = T, sep = "\t") %>%
  mutate(Gene.name = ifelse(Gene.name == "", geneSymbol, Gene.name))



O_T_data <- results %>%
  mutate(
    logFC_group = as.factor(
      case_when(
        FDR > 0.05 ~ "No significance",
        logFC > 1 & FDR <= 0.05 ~ "Male-biased gene",
        logFC < -1 & FDR <= 0.05 ~ "Female-biased gene",
        TRUE ~ "No significance"
      )
    ),
    Chromosome.scaffold.name = as.factor(Chromosome.scaffold.name),
    negative_P_value = -log10(PValue + 1e-300),
    label = ifelse(Gene.name %in% results$Gene.name[which(FDR <= 0.05 & abs(logFC) > 1)], Gene.name, NA),
    color_by_FDR = ifelse(FDR <= 0.05 & abs(logFC) > 1, "Above_Threshold", "Below_Threshold")
  ) %>%
  dplyr::select(Gene.name, Gene.type, Chromosome.scaffold.name, negative_P_value, logFC, logFC_group, label, color_by_FDR, FDR) # 컬럼 선택에 color_by_FDR 추가


O_T_data$Chromosome.scaffold.name <- factor(O_T_data$Chromosome.scaffold.name, levels = c(as.character(seq_len(28)),"W","Z"))
# 

# Define colors for each group
color_labels <- c("Male-biased gene" = "#604DA0", "Female-biased gene" = "#F7A31C", "No significance" = "#BDBDBD")

# Plot
FB_MD_plot <- ggplot(O_T_data, aes(x = Chromosome.scaffold.name, y = negative_P_value, color = logFC_group)) +
  geom_segment(aes(x = Chromosome.scaffold.name, xend = Chromosome.scaffold.name, y = 0, yend = 300), alpha = 1, size = 9, color = "#f2f2f2", lineend = "round") +
  geom_jitter(size = 5, alpha = 0.8, position = position_jitterdodge(dodge.width = 0.75, jitter.width = 0.35)) +
  scale_color_manual(values = color_labels) +
  theme_classic() +
  scale_y_continuous(expand = c(0, 0)) +
  labs(y = "-log10 (P-value) of the test for difference between reproductive organs") +
  coord_cartesian(ylim = c(0, 310)) +
  theme(axis.text.x = element_text(size = 13, hjust = 0.5),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = 12),
        legend.position = "none")

FB_MD_plot

ggsave("ManhattanPlot.pdf",width=10,height=4)
#### Figure 4B-4E ####

### Load data ###
load(file="./Data/O_T.RData")


idx <- c()


idx <- c(idx,which(Annotation$Gene.name=="AMH" | Annotation$Gene.name=="DMRT1" | Annotation$geneSymbol=="ENSCJPG00005003198"|Annotation$Gene.name=="SOX9")
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



for(i in 1:length(temp_name)){
  
  # 필터링된 데이터 생성
  filtered_data <- plot_data[plot_data$GeneSymbol %in% temp_name[i], ]
  
  #ed6a66
  # 박스플롯 및 지터 플롯 생성
  p2 <- ggplot(filtered_data, aes(x = Tissue, y = Expression, fill = Tissue)) +
    geom_boxplot(outlier.shape = NA,alpha=0.7) +  # 박스플롯 설정
    geom_jitter(aes(color = Tissue),alpha = 0.5, size = 5) +  # 지터 플롯 추가
    scale_fill_manual(values = c("Ovary" = "#F5D042", "Testis" = "#0A174E")) +  # 채우기 색상 설정
    scale_color_manual(values = c("Ovary" = "#F5D042", "Testis" = "#0A174E")) +  # 점 색상 설정
    theme_classic2()+
    theme(
      legend.position = "none",
      axis.title.x = element_blank(),
      axis.title.y = element_blank(),
      axis.text = element_text(size=12, color="black"),
      strip.text = element_text(size=12, color="black")
    )
  
  p2
  
  ggsave(paste0(temp_name[i],".pdf"),width=3,height=3)
  
  
  
} 



#### Figure 4F ####

data <- fread("./Data/Result_O_T.txt")

data <- data[data$FDR<=0.05,]

temp <- data

temp_table <- table(temp$Gene.type)
chromosome_df <- as.data.frame(temp_table)
chromosome_df$percentage <- with(chromosome_df, round(Freq / sum(Freq) * 100, 4))


chromosome_df$Var1 <- factor(chromosome_df$Var1,levels=c("protein_coding","lncRNA","miRNA","misc_RNA","pseudogene","rRNA","scaRNA","snoRNA","snRNA"))

# Color mapping

color <- c( "#00539C",  "#EEA47F",  "#ee6850",  "#eb5b6f",  "#ea7baf",  "#b198c6",  "#f9bb31",  "#8aaf6f",  "#1ba2ac","#f17430")

ggplot() +
  geom_bar(data=chromosome_df, aes(x = 1.3, y = percentage, fill = Var1), stat = "identity", width = 1,color="black",size=1) +
  coord_polar(theta = "y") +
  scale_fill_manual(values = color) +
  theme_void() +
  theme(legend.position = "right")

ggsave("Piechart_genetype.pdf",width=3,height=3)




#### Figure 4G ####

load(file="./Data/O_T.RData")


data <- merge(Annotation,Result_Gonad_table,by.x="geneSymbol",by.y="genes",sort=F)


#### Sex chromosome ###


idx <- data[which(data$Chromosome.scaffold.name=="Z"),]


temp <- idx[idx$FDR<=0.05,]


data$Chromosome.scaffold.name <- factor(data$Chromosome.scaffold.name,levels = c(1:28,"W","Z"))


idx <- which(data$Chromosome.scaffold.name=="Z")


temp <- data[idx,]



# ggplot을 사용하여 박스 플롯 그리기
ggplot(data, aes(x = Chromosome.scaffold.name, y = logFC, fill = Chromosome.scaffold.name)) +
  geom_boxplot(outliers = FALSE,alpha=0.8) +  # 박스플롯
  scale_fill_viridis(discrete = TRUE, option = "D") +  # viridis 색상 그라데이션 적용
  theme_classic() +
  geom_hline(yintercept =  mean(data$logFC), linetype="dashed",color = "red", size = 1.3,alpha=0.8) +  # logFC -1에 빨간색 점선 추가
  theme(legend.position = "None",
        axis.text.x=element_text(size=15,color="black"),
        axis.text.y=element_text(size=15,color="black"),
        axis.title.x=element_blank(),
        axis.title.y=element_blank(),
        
        
  )  # x축 레이블 회전 및 조정



ggsave("Boxplot_chromosome.pdf",width=10,height=4)

#### Figure 4H ####

options(scipen=30)
color_label <- function(data){
  label <- c()
  
  color <- c()
  #ABF200
  for(i in 1:nrow(data)){
    
    if(abs(data$logFC[i]) == 0){
      label <- c(label,"logFC=0")
    }else if(abs(data$logFC[i]) > 0  && abs(data$logFC[i]) <= 0.5){
      label <- c(label,"0<logFC<0.5")
    } else if(abs(data$logFC[i]) > 0.5  && abs(data$logFC[i]) <= 1){
      label <- c(label,"0.5<logFC<1")
    } else if(abs(data$logFC[i]) > 1 && abs(data$logFC[i]) <= 2 ){
      label <- c(label,"1<logFC<2")
    }else if(abs(data$logFC[i]) > 2 && abs(data$logFC[i]) <= 3){
      label <- c(label,"2<logFC<3")
    }else if(abs(data$logFC[i]) > 3 && abs(data$logFC[i]) <= 4){
      label <- c(label,"3<logFC<4")
    }else if(abs(data$logFC[i]) > 4 && abs(data$logFC[i]) <= 5){
      label <- c(label,"4<logFC<5")
    }else if(abs(data$logFC[i]) > 5 && abs(data$logFC[i]) <= 6){
      label <- c(label,"5<logFC<6")
    }else if(abs(data$logFC[i]) > 6){
      label <- c(label,"logFC>6")
    }
    
  }
  
  return(label)
}

data <- fread("./Data/Result_O_T.txt")


Chromosome_index <- which(data$Chromosome.scaffold.name=="Z")


data <- data[Chromosome_index]


label_list <- color_label(data)

data$label <- factor(label_list,levels = c("0<logFC<0.5","0.5<logFC<1","1<logFC<2","2<logFC<3","3<logFC<4","4<logFC<5","5<logFC<6","logFC>6"))

# Order factor
data <- data %>% arrange(desc(label))

# Replace geneSymbol
data$Gene.name[which(data$Gene.name == "")] <- data$geneSymbol[which(data$Gene.name == "")]

# Set palette
num_colors <- length(unique(data$label))
colors <- rev(brewer.pal(num_colors, "RdYlGn"))

Dosage <- ggplot(data=data,aes(x=Gene.start..bp.,y=logFC,label=Gene.name,color=label))+
  geom_jitter(alpha=0.8,size=5)+
  scale_color_manual(values=colors)+
  coord_cartesian(ylim = c(-14, 14))+
  geom_smooth(method = "loess",aes(group=1),color="black",size=2)+
  geom_hline(yintercept = mean(data$logFC),size=2,alpha=0.6,linetype='longdash',color="red")+
  geom_text_repel(fontface="italic",color="black")+
  theme_bw()+
  theme_classic()+
  theme(
    axis.title = element_blank(),
    axis.text = element_text(size=12,color="black"),
    strip.text=element_text(size=12,color="black"),
    legend.position = "None"
  )

Dosage

ggsave("ScatterPlot.pdf",width=15,height=4)
#### Figure 4I ####
# Preprocessing of anootation result #
column <- c("Category","Term","Count","%","PValue","Genes","List Total","Pop Hits","Pop Total","Fold Enrichment","Bonferroni","Benjamini","FDR")

lines <- readLines("./Data/Figure4_data/Male_up.txt")

lines <- lines[!grepl("Annotation Cluster", lines)]


lines <- lines[!grepl("Category\t", lines)]

df <- data.frame()

# 각 줄을 데이터프레임에 추가
for (line in lines) {
  # 탭으로 분리하여 데이터프레임에 추가
  row <- unlist(strsplit(line, "\t"))
  df <- rbind(df, row)
}

colnames(df) <- column

df$FDR <- as.numeric(df$FDR)
df$PValue <- as.numeric(df$PValue)


df <- df[df$FDR<=0.01,]


## Load annotation data ##

data <- fread("./Data/Figure4_data/Male_up_bar_plot_data.txt")

data <- data[1:8,]


palette <- c("#ff595e","#ff924c", "#ffca3a",  "#c5ca30",  "#8ac926",  "#52a675",  "#1982c4",  "#4267ac",  "#6a4c93","#b5a6c9")

# data.table을 사용하여 Term을 PValue에 따라 내림차순으로 정렬
setorder(data, PValue)

data$Term <- factor(data$Term, levels = unique(data$Term))

Zplot <- ggplot(data, aes(x=Term, y=-log10(data$PValue), fill=Term)) +
  geom_bar(stat="identity", colour="black", size=1,width=1) +  # geom_bar를 사용하여 막대 그래프 생성
  scale_fill_manual(values = palette) +  # 색상 팔레트 적용
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


#### Figure 4J ####
# Preprocess #
column <- c("Category","Term","Count","%","PValue","Genes","List Total","Pop Hits","Pop Total","Fold Enrichment","Bonferroni","Benjamini","FDR")

lines <- readLines("./Data/Figure4_data/Female_up.txt")

lines <- lines[!grepl("Annotation Cluster", lines)]


lines <- lines[!grepl("Category\t", lines)]


df <- data.frame()

# 각 줄을 데이터프레임에 추가
for (line in lines) {
  # 탭으로 분리하여 데이터프레임에 추가
  row <- unlist(strsplit(line, "\t"))
  df <- rbind(df, row)
}

colnames(df) <- column

df$FDR <- as.numeric(df$FDR)
df$PValue <- as.numeric(df$PValue)


df <- df[df$FDR<=0.01,]

## Load Bar plot data ##

data <- fread("./Data/Figure4_data/Female_up_bar_plot_data.txt")


palette <- c("#536e0a",  "#7da50f",  "#bdda0f",  "#ffa800",  "#ff7a00",  "#ff3d35",  "#e52b6f",  "#6226a9"  ,"#2c29a2")
# data.table을 사용하여 Term을 PValue에 따라 내림차순으로 정렬
setorder(data, PValue)


data$Term <- factor(data$Term, levels = unique(data$Term))


Zplot <- ggplot(data, aes(x=Term, y=-log10(data$PValue), fill=Term)) +
  geom_bar(stat="identity", colour="black", size=1,width=0.95) +  # geom_bar를 사용하여 막대 그래프 생성
  scale_fill_manual(values = palette) +  # 색상 팔레트 적용
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

ggsave("GOterm.pdf",width=4,height=3)