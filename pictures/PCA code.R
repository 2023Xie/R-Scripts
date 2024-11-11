rm(list=ls())#clear Global Environment
library(vegan)#计算距离时需要的包
library(ggplot2)#绘图包
library(ggpubr)
library(ggsignif)
library(vegan) 
library(GUniFrac)
library(ape)
library(ade4)

# 颜色
color=c("#0E606B","#1597A5","#FEB3AE")
color=c("#1597A5","#FFC24B","#FEB3AE")
color=c("#264653","#287271","#2a9d8c","#8ab07d")
color=c("#264653","#287271","#2a9d8c","#8ab07d","#e9c46b")
color=c("#264653","#287271","#2a9d8c","#8ab07d","#e9c46b","#f3a261","#e66f51")

# 分组顺序
group_order <- c("CON", "T2DM", "MET")
group_order <- c("CON", "HFD", "LFU", "HFU")
group_order <- c("CON", "HFD", "AKK", "AKK+LFU", "AKK+HFU")
group_order <- c("CON", "HFD", "AKK", "HFU", "LFU", "AKK+HFU", "AKK+LFU")

# 比较组
comparisons <- list(c("CON", "T2DM"), c("CON", "MET"), c("T2DM", "MET"))
comparisons <- list(c("CON", "HFU"), c("CON", "LFU"), c("HFU", "LFU"))
comparisons <- list(c("CON", "HFD"), c("CON", "AKK"), c("CON", "LFU"), c("CON", "HFU"), c("CON", "AKK+HFU"), c("CON", "AKK+LFU"), c("HFD", "HFU"), c("HFD", "LFU"), c("HFD", "AKK+HFU"), c("HFD", "AKK+LFU"), c("HFD", "AKK"), c("AKK", "AKK+LFU"), c("AKK", "AKK+HFU"), c("AKK", "LFU"), c("AKK", "HFU"), c("LFU", "HFU"), c("LFU", "AKK+HFU"), c("LFU", "AKK+LFU"), c("HFU", "AKK+HFU"), c("HFU", "AKK+LFU"), c("AKK+HFU", "AKK+LFU"))

# 显著性阈值
alpha_level <- 0.05

# 导入数据
otu_table <- read.table("feature-table.csv", header = T, row.names = 1, comment.char = "", quote = "", stringsAsFactors = F, sep = ",")
name <- read.table("feature-table.csv", header = F, row.names = 1, comment.char = "", quote = "", stringsAsFactors = F, sep = ",")
names(otu_table) <- name[1, ] 
rm(name)

#以最小丰度值为依据，进行抽平
otu_table <- as.data.frame(t(rrarefy(t(otu_table), min(colSums(otu_table)))))
#查看抽平后的每个群落的总丰度
colSums(otu_table)

# 转置数据
otu <- t(otu_table)

# 导入分组
group <- read.table("group.txt", header = T, row.names = NULL, comment.char = "", sep = "\t", quote = "", stringsAsFactors = F)
names(group) <- c("sample", "group")
group$group <- as.factor(group$group)

# 导入树文件
tree <- read.tree("tree.nwk")

# 计算距离
bray <- vegdist(otu, method = 'bray')
jaccard <- vegdist(otu, method = 'jaccard')
unifracs <- GUniFrac(otu, tree, alpha = c(0, 0.5, 1))$unifracs
weight <- unifracs[, , "d_0.5"]
unweight <- unifracs[, , "d_UW"]		# Unweighted UniFrac

# 使用prcomp()函数进行PCA分析
pca <- prcomp(bray, scale. = TRUE)
data <- as.data.frame(pca$x)
data$sample <- rownames(data)
data <- merge(data, group, by = "sample")
# 绘制PCA图表
pc <- ggplot(data, aes(x = PC1, y = PC2)) +
  geom_point(size = 3) +
  theme_bw()
pc
pc1 <- ggplot(data = data, aes(x = PC1, y = PC2)) +
  theme_bw() +
  geom_point(aes(color = group), shape = 19, size = 3) +
  theme(panel.grid = element_blank(), legend.position = "bottom") +
  geom_vline(xintercept = 0, lty = "dashed", size = 1, color = 'grey50') +
  geom_hline(yintercept = 0, lty = "dashed", size = 1, color = 'grey50') +
  geom_text(aes(label = sample, y = PC2 + 0.03, x = PC1 + 0.03,
                vjust = 0, color = group), size = 3.5, show.legend = F, alpha = 0.7) +
  stat_ellipse(data = data,
               geom = "polygon", level = 0.95,
               linetype = 2, size = 0.5,
               aes(fill = group),
               alpha = 0.2) +
  scale_color_manual(values = color, breaks = group_order) +
  scale_fill_manual(values = color, breaks = group_order) +
  labs(x = "PC1", y = "PC2") +
  theme(axis.title.x = element_text(size = 12),
        axis.title.y = element_text(size = 12, angle = 90),
        axis.text.y = element_text(size = 10),
        axis.text.x = element_text(size = 10),
        panel.grid = element_blank())
pc1
#anosim
group1 <- factor(group$group)
pca_coord <- pca$x[,1:2]
anosim_pca <- anosim(dist(pca_coord), group1)
# Extract R and P values
R_value <- round(anosim_pca$statistic, 2)
P_value <- round(anosim_pca$signif, 3)
#添加P和R
pc2 <- pc1 + annotate("text", x = -4, y = 8, label = paste0("R = ", R_value, ", P = ", P_value), size = 5, fontface = "bold")
pc2
# 计算并添加贡献率标签
variance <- pca$sdev^2 / sum(pca$sdev^2)
pc3 <- pc2 + labs(x = paste0("PC1 (", round(variance[1] * 100, 1), "%)"), y = paste0("PC2 (", round(variance[2] * 100, 1), "%)"))
pc3

# 绘制y轴为PC2值的分组箱线图
pc4 <- ggplot(data,aes(x= factor(group, levels = group_order),y=PC2))+
  stat_boxplot(geom = "errorbar", width=0.1,size=0.5)+
  geom_boxplot(aes(fill=group), 
               outlier.colour="white",size=0.5)+
  theme(panel.background =element_blank(), 
        axis.line=element_line(color = "white"),
        axis.text.y = element_blank(),
        panel.grid = element_blank(),
        axis.ticks = element_blank(),
        legend.position = 'none')+
  xlab("") + ylab("")+
  scale_fill_manual(values= color, breaks = group_order)+
  geom_signif(
    comparisons = comparisons[which(sapply(comparisons, function(c) wilcox.test(data$PC2[data$group == c[1]], data$PC2[data$group == c[2]])$p.value < alpha_level))],
    map_signif_level = T, 
    test = t.test,
    y_position = c(4,5,4.5),
    tip_length = c(c(0,0),
                   c(0,0),
                   c(0,0)),
    size=0.6,color="black")
pc4
# 绘制y轴为PC1值的分组箱线图
pc5 <- ggplot(data,aes(x = factor(group, levels = group_order),y=PC1))+
  stat_boxplot(geom = "errorbar", width=0.1,size=0.5)+
  coord_flip()+
  geom_boxplot(aes(fill=group), 
               outlier.colour="white",size=0.5)+
  theme(panel.background =element_blank(), 
        axis.line=element_line(color = "white"),
        axis.text.x = element_blank(),
        panel.grid = element_blank(),
        axis.ticks = element_blank(),
        legend.position = 'none')+
  xlab("") + ylab("")+
  scale_fill_manual(values= color, breaks = group_order)+
  geom_signif(
    comparisons = comparisons[which(sapply(comparisons, function(c) wilcox.test(data$PC1[data$group == c[1]], data$PC1[data$group == c[2]])$p.value < alpha_level))],
    map_signif_level = T,
    test = t.test,
    y_position = c(7,8,7.5),
    tip_length = c(c(0,0),
                   c(0,0),
                   c(0,0)),
    size=0.6,color="black")
pc5
# ggpubr::ggarrange()函数对图进行拼接
ggarrange(pc5, NULL, pc3, pc4, widths = c(6,2), heights = c(2,6), align = "hv")
ggsave("PCA_bray.tiff", height = 7, width = 8, dpi = 600)

# 使用prcomp()函数进行PCA分析
pca <- prcomp(jaccard, scale. = TRUE)
data <- as.data.frame(pca$x)
data$sample <- rownames(data)
data <- merge(data, group, by = "sample")
# 绘制PCA图表
pc <- ggplot(data, aes(x = PC1, y = PC2)) +
  geom_point(size = 3) +
  theme_bw()
pc
pc1 <- ggplot(data = data, aes(x = PC1, y = PC2)) +
  theme_bw() +
  geom_point(aes(color = group), shape = 19, size = 3) +
  theme(panel.grid = element_blank(), legend.position = "bottom") +
  geom_vline(xintercept = 0, lty = "dashed", size = 1, color = 'grey50') +
  geom_hline(yintercept = 0, lty = "dashed", size = 1, color = 'grey50') +
  geom_text(aes(label = sample, y = PC2 + 0.03, x = PC1 + 0.03,
                vjust = 0, color = group), size = 3.5, show.legend = F, alpha = 0.7) +
  stat_ellipse(data = data,
               geom = "polygon", level = 0.95,
               linetype = 2, size = 0.5,
               aes(fill = group),
               alpha = 0.2) +
  scale_color_manual(values = color, breaks = group_order) +
  scale_fill_manual(values = color, breaks = group_order) +
  labs(x = "PC1", y = "PC2") +
  theme(axis.title.x = element_text(size = 12),
        axis.title.y = element_text(size = 12, angle = 90),
        axis.text.y = element_text(size = 10),
        axis.text.x = element_text(size = 10),
        panel.grid = element_blank())
pc1
#anosim
group1 <- factor(group$group)
pca_coord <- pca$x[,1:2]
anosim_pca <- anosim(dist(pca_coord), group1)
# Extract R and P values
R_value <- round(anosim_pca$statistic, 2)
P_value <- round(anosim_pca$signif, 3)
#添加P和R
pc2 <- pc1 + annotate("text", x = -4, y = 8, label = paste0("R = ", R_value, ", P = ", P_value), size = 5, fontface = "bold")
pc2
# 计算并添加贡献率标签
variance <- pca$sdev^2 / sum(pca$sdev^2)
pc3 <- pc2 + labs(x = paste0("PC1 (", round(variance[1] * 100, 1), "%)"), y = paste0("PC2 (", round(variance[2] * 100, 1), "%)"))
pc3

# 绘制y轴为PC2值的分组箱线图
pc4 <- ggplot(data,aes(x= factor(group, levels = group_order),y=PC2))+
  stat_boxplot(geom = "errorbar", width=0.1,size=0.5)+
  geom_boxplot(aes(fill=group), 
               outlier.colour="white",size=0.5)+
  theme(panel.background =element_blank(), 
        axis.line=element_line(color = "white"),
        axis.text.y = element_blank(),
        panel.grid = element_blank(),
        axis.ticks = element_blank(),
        legend.position = 'none')+
  xlab("") + ylab("")+
  scale_fill_manual(values= color, breaks = group_order)+
  geom_signif(
    comparisons = comparisons[which(sapply(comparisons, function(c) wilcox.test(data$PC2[data$group == c[1]], data$PC2[data$group == c[2]])$p.value < alpha_level))],
    map_signif_level = T, 
    test = t.test,
    y_position = c(4,5,4.5),
    tip_length = c(c(0,0),
                   c(0,0),
                   c(0,0)),
    size=0.6,color="black")
pc4
# 绘制y轴为PC1值的分组箱线图
pc5 <- ggplot(data,aes(x = factor(group, levels = group_order),y=PC1))+
  stat_boxplot(geom = "errorbar", width=0.1,size=0.5)+
  coord_flip()+
  geom_boxplot(aes(fill=group), 
               outlier.colour="white",size=0.5)+
  theme(panel.background =element_blank(), 
        axis.line=element_line(color = "white"),
        axis.text.x = element_blank(),
        panel.grid = element_blank(),
        axis.ticks = element_blank(),
        legend.position = 'none')+
  xlab("") + ylab("")+
  scale_fill_manual(values= color, breaks = group_order)+
  geom_signif(
    comparisons = comparisons[which(sapply(comparisons, function(c) wilcox.test(data$PC1[data$group == c[1]], data$PC1[data$group == c[2]])$p.value < alpha_level))],
    map_signif_level = T,
    test = t.test,
    y_position = c(7,8,7.5),
    tip_length = c(c(0,0),
                   c(0,0),
                   c(0,0)),
    size=0.6,color="black")
pc5
# ggpubr::ggarrange()函数对图进行拼接
ggarrange(pc5, NULL, pc3, pc4, widths = c(6,2), heights = c(2,6), align = "hv")
ggsave("PCA_jaccard.tiff", height = 7, width = 8, dpi = 600)

# 使用prcomp()函数进行PCA分析
pca <- prcomp(weight, scale. = TRUE)
data <- as.data.frame(pca$x)
data$sample <- rownames(data)
data <- merge(data, group, by = "sample")
# 绘制PCA图表
pc <- ggplot(data, aes(x = PC1, y = PC2)) +
  geom_point(size = 3) +
  theme_bw()
pc
pc1 <- ggplot(data = data, aes(x = PC1, y = PC2)) +
  theme_bw() +
  geom_point(aes(color = group), shape = 19, size = 3) +
  theme(panel.grid = element_blank(), legend.position = "bottom") +
  geom_vline(xintercept = 0, lty = "dashed", size = 1, color = 'grey50') +
  geom_hline(yintercept = 0, lty = "dashed", size = 1, color = 'grey50') +
  geom_text(aes(label = sample, y = PC2 + 0.03, x = PC1 + 0.03,
                vjust = 0, color = group), size = 3.5, show.legend = F, alpha = 0.7) +
  stat_ellipse(data = data,
               geom = "polygon", level = 0.95,
               linetype = 2, size = 0.5,
               aes(fill = group),
               alpha = 0.2) +
  scale_color_manual(values = color, breaks = group_order) +
  scale_fill_manual(values = color, breaks = group_order) +
  labs(x = "PC1", y = "PC2") +
  theme(axis.title.x = element_text(size = 12),
        axis.title.y = element_text(size = 12, angle = 90),
        axis.text.y = element_text(size = 10),
        axis.text.x = element_text(size = 10),
        panel.grid = element_blank())
pc1
#anosim
group1 <- factor(group$group)
pca_coord <- pca$x[,1:2]
anosim_pca <- anosim(dist(pca_coord), group1)
# Extract R and P values
R_value <- round(anosim_pca$statistic, 2)
P_value <- round(anosim_pca$signif, 3)
#添加P和R
pc2 <- pc1 + annotate("text", x = -4, y = 8, label = paste0("R = ", R_value, ", P = ", P_value), size = 5, fontface = "bold")
pc2
# 计算并添加贡献率标签
variance <- pca$sdev^2 / sum(pca$sdev^2)
pc3 <- pc2 + labs(x = paste0("PC1 (", round(variance[1] * 100, 1), "%)"), y = paste0("PC2 (", round(variance[2] * 100, 1), "%)"))
pc3

# 绘制y轴为PC2值的分组箱线图
pc4 <- ggplot(data,aes(x= factor(group, levels = group_order),y=PC2))+
  stat_boxplot(geom = "errorbar", width=0.1,size=0.5)+
  geom_boxplot(aes(fill=group), 
               outlier.colour="white",size=0.5)+
  theme(panel.background =element_blank(), 
        axis.line=element_line(color = "white"),
        axis.text.y = element_blank(),
        panel.grid = element_blank(),
        axis.ticks = element_blank(),
        legend.position = 'none')+
  xlab("") + ylab("")+
  scale_fill_manual(values= color, breaks = group_order)+
  geom_signif(
    comparisons = comparisons[which(sapply(comparisons, function(c) wilcox.test(data$PC2[data$group == c[1]], data$PC2[data$group == c[2]])$p.value < alpha_level))],
    map_signif_level = T, 
    test = t.test,
    y_position = c(4,5,4.5),
    tip_length = c(c(0,0),
                   c(0,0),
                   c(0,0)),
    size=0.6,color="black")
pc4
# 绘制y轴为PC1值的分组箱线图
pc5 <- ggplot(data,aes(x = factor(group, levels = group_order),y=PC1))+
  stat_boxplot(geom = "errorbar", width=0.1,size=0.5)+
  coord_flip()+
  geom_boxplot(aes(fill=group), 
               outlier.colour="white",size=0.5)+
  theme(panel.background =element_blank(), 
        axis.line=element_line(color = "white"),
        axis.text.x = element_blank(),
        panel.grid = element_blank(),
        axis.ticks = element_blank(),
        legend.position = 'none')+
  xlab("") + ylab("")+
  scale_fill_manual(values= color, breaks = group_order)+
  geom_signif(
    comparisons = comparisons[which(sapply(comparisons, function(c) wilcox.test(data$PC1[data$group == c[1]], data$PC1[data$group == c[2]])$p.value < alpha_level))],
    map_signif_level = T,
    test = t.test,
    y_position = c(7,8,7.5),
    tip_length = c(c(0,0),
                   c(0,0),
                   c(0,0)),
    size=0.6,color="black")
pc5
# ggpubr::ggarrange()函数对图进行拼接
ggarrange(pc5, NULL, pc3, pc4, widths = c(6,2), heights = c(2,6), align = "hv")
ggsave("PCA_weight.tiff", height = 7, width = 8, dpi = 600)

# 使用prcomp()函数进行PCA分析
pca <- prcomp(unweight, scale. = TRUE)
data <- as.data.frame(pca$x)
data$sample <- rownames(data)
data <- merge(data, group, by = "sample")
# 绘制PCA图表
pc <- ggplot(data, aes(x = PC1, y = PC2)) +
  geom_point(size = 3) +
  theme_bw()
pc
pc1 <- ggplot(data = data, aes(x = PC1, y = PC2)) +
  theme_bw() +
  geom_point(aes(color = group), shape = 19, size = 3) +
  theme(panel.grid = element_blank(), legend.position = "bottom") +
  geom_vline(xintercept = 0, lty = "dashed", size = 1, color = 'grey50') +
  geom_hline(yintercept = 0, lty = "dashed", size = 1, color = 'grey50') +
  geom_text(aes(label = sample, y = PC2 + 0.03, x = PC1 + 0.03,
                vjust = 0, color = group), size = 3.5, show.legend = F, alpha = 0.7) +
  stat_ellipse(data = data,
               geom = "polygon", level = 0.95,
               linetype = 2, size = 0.5,
               aes(fill = group),
               alpha = 0.2) +
  scale_color_manual(values = color, breaks = group_order) +
  scale_fill_manual(values = color, breaks = group_order) +
  labs(x = "PC1", y = "PC2") +
  theme(axis.title.x = element_text(size = 12),
        axis.title.y = element_text(size = 12, angle = 90),
        axis.text.y = element_text(size = 10),
        axis.text.x = element_text(size = 10),
        panel.grid = element_blank())
pc1
#anosim
group1 <- factor(group$group)
pca_coord <- pca$x[,1:2]
anosim_pca <- anosim(dist(pca_coord), group1)
# Extract R and P values
R_value <- round(anosim_pca$statistic, 2)
P_value <- round(anosim_pca$signif, 3)
#添加P和R
pc2 <- pc1 + annotate("text", x = -4, y = 8, label = paste0("R = ", R_value, ", P = ", P_value), size = 5, fontface = "bold")
pc2
# 计算并添加贡献率标签
variance <- pca$sdev^2 / sum(pca$sdev^2)
pc3 <- pc2 + labs(x = paste0("PC1 (", round(variance[1] * 100, 1), "%)"), y = paste0("PC2 (", round(variance[2] * 100, 1), "%)"))
pc3

# 绘制y轴为PC2值的分组箱线图
pc4 <- ggplot(data,aes(x= factor(group, levels = group_order),y=PC2))+
  stat_boxplot(geom = "errorbar", width=0.1,size=0.5)+
  geom_boxplot(aes(fill=group), 
               outlier.colour="white",size=0.5)+
  theme(panel.background =element_blank(), 
        axis.line=element_line(color = "white"),
        axis.text.y = element_blank(),
        panel.grid = element_blank(),
        axis.ticks = element_blank(),
        legend.position = 'none')+
  xlab("") + ylab("")+
  scale_fill_manual(values= color, breaks = group_order)+
  geom_signif(
    comparisons = comparisons[which(sapply(comparisons, function(c) wilcox.test(data$PC2[data$group == c[1]], data$PC2[data$group == c[2]])$p.value < alpha_level))],
    map_signif_level = T, 
    test = t.test,
    y_position = c(4,5,4.5),
    tip_length = c(c(0,0),
                   c(0,0),
                   c(0,0)),
    size=0.6,color="black")
pc4
# 绘制y轴为PC1值的分组箱线图
pc5 <- ggplot(data,aes(x = factor(group, levels = group_order),y=PC1))+
  stat_boxplot(geom = "errorbar", width=0.1,size=0.5)+
  coord_flip()+
  geom_boxplot(aes(fill=group), 
               outlier.colour="white",size=0.5)+
  theme(panel.background =element_blank(), 
        axis.line=element_line(color = "white"),
        axis.text.x = element_blank(),
        panel.grid = element_blank(),
        axis.ticks = element_blank(),
        legend.position = 'none')+
  xlab("") + ylab("")+
  scale_fill_manual(values= color, breaks = group_order)+
  geom_signif(
    comparisons = comparisons[which(sapply(comparisons, function(c) wilcox.test(data$PC1[data$group == c[1]], data$PC1[data$group == c[2]])$p.value < alpha_level))],
    map_signif_level = T,
    test = t.test,
    y_position = c(7,8,7.5),
    tip_length = c(c(0,0),
                   c(0,0),
                   c(0,0)),
    size=0.6,color="black")
pc5
# ggpubr::ggarrange()函数对图进行拼接
ggarrange(pc5, NULL, pc3, pc4, widths = c(6,2), heights = c(2,6), align = "hv")
ggsave("PCA_unweight.tiff", height = 7, width = 8, dpi = 600)

