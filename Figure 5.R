#### Result 5 code for validation ####

library(ggpubr)
library(ggplot2)
library(data.table)
library(RColorBrewer)
library(ggrepel)
library(randomcoloR)
library(dplyr)
library(edgeR)
library(DESeq2)
library(data.table)
#### Main plot ####
cut <- 1
getCol <- function(data){
  
  
  palette <- c("#FCE05E","#F68D51","#AC2324","#BDD74B","#7CB5D2","#6574cd","#9561e2","#18A68D","grey30")
  
  
  plotCol <- rep("#000000",nrow(data))
  plotLab <- rep("",nrow(data))
  
  
  
  plotCol[which(data$logFC_gonad>cut & data$logFC_brain>cut)] <- palette[1]
  plotLab[which(data$logFC_gonad>cut & data$logFC_brain>cut)] <- "Section 1"
  
  
  plotCol[which((data$logFC_gonad < cut & data$logFC_gonad> -cut) & data$logFC_brain > cut)] <- palette[2]
  plotLab[which((data$logFC_gonad < cut & data$logFC_gonad> -cut) & data$logFC_brain > cut)] <- "Section 2"
  
  
  plotCol[which(data$logFC_gonad< -cut & data$logFC_brain>cut)] <- palette[3]
  plotLab[which(data$logFC_gonad< -cut & data$logFC_brain>cut)] <- "Section 3"
  
  
  
  plotCol[which((data$logFC_brain<cut & data$logFC_brain> -cut) & data$logFC_gonad < -cut)] <- palette[4]
  plotLab[which((data$logFC_brain<cut & data$logFC_brain> -cut) & data$logFC_gonad < -cut)] <- "Section 4"
  
  plotCol[which(data$logFC_gonad< -cut & data$logFC_brain< -cut)] <- palette[5]
  plotLab[which(data$logFC_gonad< -cut & data$logFC_brain< -cut)] <- "Section 5"
  
  
  plotCol[which((data$logFC_gonad<cut & data$logFC_gonad> -cut) & data$logFC_brain< -cut)] <- palette[6]
  plotLab[which((data$logFC_gonad<cut & data$logFC_gonad> -cut) & data$logFC_brain< cut)] <- "Section 6"
  
  plotCol[which(data$logFC_gonad>cut & data$logFC_brain< -cut)] <- palette[7]
  plotLab[which(data$logFC_gonad>cut & data$logFC_brain< -cut)] <- "Section 7"
  
  
  
  
  plotCol[which((data$logFC_brain<cut & data$logFC_brain> -cut) & data$logFC_gonad>cut)] <- palette[8]
  plotLab[which((data$logFC_brain<cut & data$logFC_brain> -cut) & data$logFC_gonad>cut)] <- "Section 8"
  
  
  
  plotCol[which((data$logFC_brain<cut & data$logFC_brain> -cut) & (data$logFC_gonad<cut & data$logFC_gonad> -cut))] <- "#BDBDBD"
  plotLab[which((data$logFC_brain<cut & data$logFC_brain> -cut) & (data$logFC_gonad<cut & data$logFC_gonad> -cut))] <- "Section 9"
  
  
  plotCol <- factor(plotCol,levels=c(palette,"#BDBDBD"))
  
  plotLab <- factor(plotLab,levels=c("Section 1","Section 2","Section 3","Section 4","Section 5","Section 6","Section 7","Section 8","Section 9"))
  
  return(list(col=plotCol,lab=plotLab))
}

# 데이터 불러오기
data <- fread("./Data/Figure5_data/Figure4_data.txt")


# 색상 열 추가
temp <- getCol(data)


data$group <- temp$lab
data$palette <- temp$col
data <- data[-which(data$group=="Section 9"),]


data$GeneName <- ifelse(is.na(data$GeneName) | data$GeneName == "", data$geneSymbol, data$GeneName)

ggplot(data, aes(x=logFC_gonad, y=logFC_brain, color=palette,label=GeneName)) +
  geom_point(size=5, alpha=0.5) +
  scale_color_identity() +
  geom_text_repel(color="black",fontface="italic")+
  geom_hline(yintercept = cut, color="black", size=1, alpha=0.4, linetype="dashed") +
  geom_hline(yintercept = -cut, color="black", size=1, alpha=0.4, linetype="dashed") +
  geom_vline(xintercept = cut, color="black", size=1,alpha=0.4, linetype="dashed") +
  geom_vline(xintercept = -cut, color="black", size=1, alpha=0.4, linetype="dashed") +
  theme_bw() +
  scale_y_continuous(expand = c(0, 0), limits = c(-15, 15)) +
  scale_x_continuous(expand = c(0, 0), limits = c(-15, 15)) +
  theme_classic() +
  theme()

#### Subplot ####

### Load table that have common genes between reproductive organ and brain tissue ###

Annotation <- fread("./Data/Figure5_data/Figure4_data.txt")
## Loading featureCount data
meta_data <- data.frame(fread("./Data/meta_data.txt", sep = "\t",header=T))


fileList <- paste0("./Data/Figure5_data/",meta_data$Run, '.txt')

temp <- data.frame(fread(fileList[1]))

temp <- temp[temp$Geneid %in% Annotation$geneSymbol,]
dim(temp)
str(temp)


all.equal(Annotation$geneSymbol,temp$Geneid)

geneSymbol <- temp$Geneid


count <- data.frame(temp[, 7])


for (i in 2:length(fileList)) {
  temp <- data.frame(fread(fileList[i]))
  temp <- temp[temp$Geneid %in% Annotation$geneSymbol,]
  count <- data.frame(count, temp[, 7])
}


## Defining experimental variables

Tissue <- factor(meta_data$Tissue, levels = c("Gonad", "Brain"))
Sex <- factor(meta_data$Sexual, levels = c("Female", "Male"))



design <- model.matrix( ~ Tissue*Sex)

## TMM normalization using edgeR package
Y <- DGEList(counts = count, genes = geneSymbol)
Y <- calcNormFactors(Y, method = "TMM") ## TMM Normalization
logCPM <- edgeR::cpm(Y, normalized.lib.sizes = TRUE, log = T)





gene_list <- c("CCNH","ENSCJPG00005003198","SLC13A4","VIP","DMRT1","TTR","RPS6","WNT4")




for(i in 1:length(gene_list)){
  idx <- which(Annotation$GeneName==gene_list[i] | Annotation$geneSymbol==gene_list[i])
  
  temp_name <- Annotation$GeneName[idx]
  
  if(temp_name==""){
    temp_name <- Annotation$geneSymbol[idx]
  }
  
  
  # Generate Plot data
  plot_data <- data.frame(
    Expression = as.numeric(logCPM[idx[1], ]),
    GeneSymbol = rep(temp_name, ncol(logCPM)),
    Tissue,
    Sex,
    Sex:Tissue
  )
  
  
  # Preprocessing
  plot_data_summary <- plot_data %>%
    group_by(Sex, Tissue) %>%
    summarise(
      Mean = mean(Expression),
      SE = sd(Expression) / sqrt(n()),
      .groups = 'drop'
    )
  
  # Drawing
  ggplot() +
    geom_jitter(data = plot_data, aes(x = Sex, y = Expression, color = Tissue), width = 0.1, height = 0, size = 5, alpha = 0.5) + # Jitter 추가
    geom_line(data = plot_data_summary, aes(x = Sex, y = Mean, group = Tissue, color = Tissue), size = 1.5) + # 선 그래프
    geom_errorbar(data = plot_data_summary, aes(x = Sex, ymin = Mean - SE, ymax = Mean + SE, group = Tissue, color = Tissue), width = 0.2) + # 표준 오차 바
    scale_color_manual(values = c("#FF7A01", "#70398C")) + # 사용자 정의 색상
    theme_classic() + # 클래식 테마
    theme(axis.title = element_blank(), # 축 타이틀 제거
          legend.position = "none") # 범례 제거
  
  ggsave(paste0(temp_name,".pdf"),width=3,height=1.5)
}





