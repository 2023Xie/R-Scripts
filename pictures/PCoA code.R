rm(list=ls())#clear Global Environment
library(vegan)#计算距离时需要的包
library(ggplot2)#绘图包
library(ggpubr)
library(ggsignif)
library(vegan) 
library(GUniFrac)
library(ape)
library(ade4)
library(readxl)



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
comparisons <- list(c("CON", "HFD"), c("CON", "LFU"), c("CON", "HFU"), c("HFD", "HFU"), c("HFD", "LFU"), c("LFU", "HFU"))
comparisons <- list(c("CON", "HFD"), c("CON", "AKK"), c("CON", "AKK+HFU"), c("CON", "AKK+LFU"), c("HFD", "AKK+HFU"), c("HFD", "AKK+LFU"), c("HFD", "AKK"), c("AKK", "AKK+LFU"), c("AKK", "AKK+HFU"), c("AKK+HFU", "AKK+LFU"))
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
weight <- as.dist(weight)
unweight <- as.dist(unweight)

# pcoa
# 使用ade4这个包中的dudi.pco ()函数做PCoA分析
pcoa = dudi.pco(bray,
                scannf = F,   # 一种逻辑值，指示是否应该显示特征值条形图
                nf=2)         # 保留几个维度的坐标信息
data = pcoa$li
data$sample = rownames(data)
data <- merge(data,group, by = "sample")
#使用ggplot2包绘图
pc <- ggplot(data,aes(x=A1,y=A2))+#指定数据、X轴、Y轴
  geom_point(size=3)+#绘制点图并设定大小
  theme_bw()#主题
pc
pc1<-ggplot(data=data,aes(x=A1,y=A2))+#指定数据、X轴、Y轴，颜色
  theme_bw()+#主题设置
  geom_point(aes(color = group), shape = 19, size=3)+#绘制点图并设定大小
  theme(panel.grid = element_blank(), legend.position = "bottom")+
  geom_vline(xintercept = 0,lty="dashed", size = 1, color = 'grey50')+
  geom_hline(yintercept = 0,lty="dashed", size = 1, color = 'grey50')+#图中虚线
  geom_text(aes(label=sample, y=A2+0.03,x=A1+0.03,
                vjust=0, color = group),size=3.5, show.legend = F, alpha = 0.7)+#添加数据点的标签
  stat_ellipse(data=data,
               geom = "polygon",level=0.95,
               linetype = 2,size=0.5,
               aes(fill=group),
               alpha=0.2)+
  scale_color_manual(values = color, breaks = group_order) +#点的颜色设置
  scale_fill_manual(values = color, breaks = group_order)+#椭圆颜色
  labs(x = "PCoA1") +
  labs(y = " PCoA2") +
  theme(axis.title.x=element_text(size=21),#修改X轴标题文本
        axis.title.y=element_text(size=21,angle=90),#修改y轴标题文本
        axis.text.y=element_text(size=15),#修改x轴刻度标签文本
        axis.text.x=element_text(size=15),#修改y轴刻度标签文本
        panel.grid=element_blank())#隐藏网格线
pc1
# 计算并添加贡献率标签
eig <- pcoa$eig # 获取每个主成分的特征值
variance <- eig / sum(eig) # 计算每个主成分的贡献率
pc2 <- pc1 + labs(x = paste0("PCoA1 (", round(variance[1] * 100, 1), "%)"), y = paste0("PCoA2 (", round(variance[2] * 100, 1), "%)"))
pc2

#anosim
anosim_result <- anosim(bray, group$group)
# Extract R and P values
R_value <- round(anosim_result $statistic, 2)
P_value <- round(anosim_result $signif, 3)
#添加P和R
pc3 <- pc2 + annotate("text", x = -0.25, y = 0.5, label = paste0("R = ", R_value, ", P = ", P_value), size = 5, fontface = "bold")
pc3


# 绘制y轴为PC2值的分组箱线图
pc4 <- ggplot(data,aes(x = factor(group, levels = group_order),y=A2))+
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
  scale_fill_manual(values = color, breaks = group_order)+
  geom_signif(
    comparisons = comparisons[which(sapply(comparisons, function(c) wilcox.test(data$A2[data$group == c[1]], data$A2[data$group == c[2]])$p.value < alpha_level))],
    map_signif_level = T, 
    test = t.test,
    y_position = c(0.3,0.4,0.35),
    tip_length = c(c(0,0),
                   c(0,0),
                   c(0,0)),
    size=0.4,color="black")
pc4
# 绘制y轴为PC1值的分组箱线图
pc5 <- ggplot(data,aes(x = factor(group, levels = group_order),y=A1))+
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
  scale_fill_manual(values = color, breaks = group_order)+
  geom_signif(
    comparisons = comparisons[which(sapply(comparisons, function(c) wilcox.test(data$A1[data$group == c[1]], data$A1[data$group == c[2]])$p.value < alpha_level))],
    map_signif_level = T,
    test = t.test,
    y_position = c(0.48,0.55,0.55),
    tip_length = c(c(0,0),
                   c(0,0),
                   c(0,0)),
    size=0.4,color="black")
pc5
# ggpubr::ggarrange()函数对图进行拼接
ggarrange(pc5, NULL, pc3, pc4, widths = c(7,2), heights = c(1.5,4), align = "hv")
ggsave("PCoA_bray.tiff", height = 8, width = 8, dpi = 600)

# pcoa
# 使用ade4这个包中的dudi.pco ()函数做PCoA分析
pcoa = dudi.pco(jaccard,
                scannf = F,   # 一种逻辑值，指示是否应该显示特征值条形图
                nf=2)         # 保留几个维度的坐标信息
data = pcoa$li
data$sample = rownames(data)
data <- merge(data,group, by = "sample")
#使用ggplot2包绘图
pc <- ggplot(data,aes(x=A1,y=A2))+#指定数据、X轴、Y轴
  geom_point(size=3)+#绘制点图并设定大小
  theme_bw()#主题
pc
pc1<-ggplot(data=data,aes(x=A1,y=A2))+#指定数据、X轴、Y轴，颜色
  theme_bw()+#主题设置
  geom_point(aes(color = group), shape = 19, size=3)+#绘制点图并设定大小
  theme(panel.grid = element_blank(), legend.position = "bottom")+
  geom_vline(xintercept = 0,lty="dashed", size = 1, color = 'grey50')+
  geom_hline(yintercept = 0,lty="dashed", size = 1, color = 'grey50')+#图中虚线
  geom_text(aes(label=sample, y=A2+0.03,x=A1+0.03,
                vjust=0, color = group),size=3.5, show.legend = F, alpha = 0.7)+#添加数据点的标签
  stat_ellipse(data=data,
               geom = "polygon",level=0.95,
               linetype = 2,size=0.5,
               aes(fill=group),
               alpha=0.2)+
  scale_color_manual(values = color, breaks = group_order) +#点的颜色设置
  scale_fill_manual(values = color, breaks = group_order)+#椭圆颜色
  labs(x = "PCoA1") +
  labs(y = " PCoA2") +
  theme(axis.title.x=element_text(size=21),#修改X轴标题文本
        axis.title.y=element_text(size=21,angle=90),#修改y轴标题文本
        axis.text.y=element_text(size=15),#修改x轴刻度标签文本
        axis.text.x=element_text(size=15),#修改y轴刻度标签文本
        panel.grid=element_blank())#隐藏网格线
pc1
# 计算并添加贡献率标签
eig <- pcoa$eig # 获取每个主成分的特征值
variance <- eig / sum(eig) # 计算每个主成分的贡献率
pc2 <- pc1 + labs(x = paste0("PCoA1 (", round(variance[1] * 100, 1), "%)"), y = paste0("PCoA2 (", round(variance[2] * 100, 1), "%)"))
pc2

#anosim
anosim_result <- anosim(jaccard, group$group)
# Extract R and P values
R_value <- round(anosim_result $statistic, 2)
P_value <- round(anosim_result $signif, 3)
#添加P和R
pc3 <- pc2 + annotate("text", x = -0.25, y = 0.5, label = paste0("R = ", R_value, ", P = ", P_value), size = 5, fontface = "bold")
pc3


# 绘制y轴为PC2值的分组箱线图
pc4 <- ggplot(data,aes(x = factor(group, levels = group_order),y=A2))+
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
  scale_fill_manual(values = color, breaks = group_order)+
  geom_signif(
    comparisons = comparisons[which(sapply(comparisons, function(c) wilcox.test(data$A2[data$group == c[1]], data$A2[data$group == c[2]])$p.value < alpha_level))],
    map_signif_level = T, 
    test = t.test,
    y_position = c(0.3,0.4,0.35),
    tip_length = c(c(0,0),
                   c(0,0),
                   c(0,0)),
    size=0.4,color="black")
pc4
# 绘制y轴为PC1值的分组箱线图
pc5 <- ggplot(data,aes(x = factor(group, levels = group_order),y=A1))+
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
  scale_fill_manual(values = color, breaks = group_order)+
  geom_signif(
    comparisons = comparisons[which(sapply(comparisons, function(c) wilcox.test(data$A1[data$group == c[1]], data$A1[data$group == c[2]])$p.value < alpha_level))],
    map_signif_level = T,
    test = t.test,
    y_position = c(0.48,0.55,0.55),
    tip_length = c(c(0,0),
                   c(0,0),
                   c(0,0)),
    size=0.4,color="black")
pc5
# ggpubr::ggarrange()函数对图进行拼接
ggarrange(pc5, NULL, pc3, pc4, widths = c(7,2), heights = c(1.5,4), align = "hv")
ggsave("PCoA_jaccard.tiff", height = 8, width = 8, dpi = 600)


# pcoa
# 使用ade4这个包中的dudi.pco ()函数做PCoA分析
pcoa = dudi.pco(weight,
                scannf = F,   # 一种逻辑值，指示是否应该显示特征值条形图
                nf=2)         # 保留几个维度的坐标信息
data = pcoa$li
data$sample = rownames(data)
data <- merge(data,group, by = "sample")
#使用ggplot2包绘图
pc <- ggplot(data,aes(x=A1,y=A2))+#指定数据、X轴、Y轴
  geom_point(size=3)+#绘制点图并设定大小
  theme_bw()#主题
pc
pc1<-ggplot(data=data,aes(x=A1,y=A2))+#指定数据、X轴、Y轴，颜色
  theme_bw()+#主题设置
  geom_point(aes(color = group), shape = 19, size=3)+#绘制点图并设定大小
  theme(panel.grid = element_blank(), legend.position = "bottom")+
  geom_vline(xintercept = 0,lty="dashed", size = 1, color = 'grey50')+
  geom_hline(yintercept = 0,lty="dashed", size = 1, color = 'grey50')+#图中虚线
  geom_text(aes(label=sample, y=A2+0.03,x=A1+0.03,
                vjust=0, color = group),size=3.5, show.legend = F, alpha = 0.7)+#添加数据点的标签
  stat_ellipse(data=data,
               geom = "polygon",level=0.95,
               linetype = 2,size=0.5,
               aes(fill=group),
               alpha=0.2)+
  scale_color_manual(values = color, breaks = group_order) +#点的颜色设置
  scale_fill_manual(values = color, breaks = group_order)+#椭圆颜色
  labs(x = "PCoA1") +
  labs(y = " PCoA2") +
  theme(axis.title.x=element_text(size=21),#修改X轴标题文本
        axis.title.y=element_text(size=21,angle=90),#修改y轴标题文本
        axis.text.y=element_text(size=15),#修改x轴刻度标签文本
        axis.text.x=element_text(size=15),#修改y轴刻度标签文本
        panel.grid=element_blank())#隐藏网格线
pc1
# 计算并添加贡献率标签
eig <- pcoa$eig # 获取每个主成分的特征值
variance <- eig / sum(eig) # 计算每个主成分的贡献率
pc2 <- pc1 + labs(x = paste0("PCoA1 (", round(variance[1] * 100, 1), "%)"), y = paste0("PCoA2 (", round(variance[2] * 100, 1), "%)"))
pc2

#anosim
anosim_result <- anosim(weight, group$group)
# Extract R and P values
R_value <- round(anosim_result $statistic, 2)
P_value <- round(anosim_result $signif, 3)
#添加P和R
pc3 <- pc2 + annotate("text", x = -0.25, y = 0.5, label = paste0("R = ", R_value, ", P = ", P_value), size = 5, fontface = "bold")
pc3


# 绘制y轴为PC2值的分组箱线图
pc4 <- ggplot(data,aes(x = factor(group, levels = group_order),y=A2))+
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
  scale_fill_manual(values = color, breaks = group_order)+
  geom_signif(
    comparisons = comparisons[which(sapply(comparisons, function(c) wilcox.test(data$A2[data$group == c[1]], data$A2[data$group == c[2]])$p.value < alpha_level))],
    map_signif_level = T, 
    test = t.test,
    y_position = c(0.3,0.4,0.35),
    tip_length = c(c(0,0),
                   c(0,0),
                   c(0,0)),
    size=0.4,color="black")
pc4
# 绘制y轴为PC1值的分组箱线图
pc5 <- ggplot(data,aes(x = factor(group, levels = group_order),y=A1))+
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
  scale_fill_manual(values = color, breaks = group_order)+
  geom_signif(
    comparisons = comparisons[which(sapply(comparisons, function(c) wilcox.test(data$A1[data$group == c[1]], data$A1[data$group == c[2]])$p.value < alpha_level))],
    map_signif_level = T,
    test = t.test,
    y_position = c(0.48,0.55,0.55),
    tip_length = c(c(0,0),
                   c(0,0),
                   c(0,0)),
    size=0.4,color="black")
pc5
# ggpubr::ggarrange()函数对图进行拼接
ggarrange(pc5, NULL, pc3, pc4, widths = c(7,2), heights = c(1.5,4), align = "hv")
ggsave("PCoA_weight.tiff", height = 8, width = 8, dpi = 600)

# pcoa
# 使用ade4这个包中的dudi.pco ()函数做PCoA分析
pcoa = dudi.pco(unweight,
                scannf = F,   # 一种逻辑值，指示是否应该显示特征值条形图
                nf=2)         # 保留几个维度的坐标信息
data = pcoa$li
data$sample = rownames(data)
data <- merge(data,group, by = "sample")
#使用ggplot2包绘图
pc <- ggplot(data,aes(x=A1,y=A2))+#指定数据、X轴、Y轴
  geom_point(size=3)+#绘制点图并设定大小
  theme_bw()#主题
pc
pc1<-ggplot(data=data,aes(x=A1,y=A2))+#指定数据、X轴、Y轴，颜色
  theme_bw()+#主题设置
  geom_point(aes(color = group), shape = 19, size=3)+#绘制点图并设定大小
  theme(panel.grid = element_blank(), legend.position = "bottom")+
  geom_vline(xintercept = 0,lty="dashed", size = 1, color = 'grey50')+
  geom_hline(yintercept = 0,lty="dashed", size = 1, color = 'grey50')+#图中虚线
  geom_text(aes(label=sample, y=A2+0.03,x=A1+0.03,
                vjust=0, color = group),size=3.5, show.legend = F, alpha = 0.7)+#添加数据点的标签
  stat_ellipse(data=data,
               geom = "polygon",level=0.95,
               linetype = 2,size=0.5,
               aes(fill=group),
               alpha=0.2)+
  scale_color_manual(values = color, breaks = group_order) +#点的颜色设置
  scale_fill_manual(values = color, breaks = group_order)+#椭圆颜色
  labs(x = "PCoA1") +
  labs(y = " PCoA2") +
  theme(axis.title.x=element_text(size=21),#修改X轴标题文本
        axis.title.y=element_text(size=21,angle=90),#修改y轴标题文本
        axis.text.y=element_text(size=15),#修改x轴刻度标签文本
        axis.text.x=element_text(size=15),#修改y轴刻度标签文本
        panel.grid=element_blank())#隐藏网格线
pc1
# 计算并添加贡献率标签
eig <- pcoa$eig # 获取每个主成分的特征值
variance <- eig / sum(eig) # 计算每个主成分的贡献率
pc2 <- pc1 + labs(x = paste0("PCoA1 (", round(variance[1] * 100, 1), "%)"), y = paste0("PCoA2 (", round(variance[2] * 100, 1), "%)"))
pc2

#anosim
anosim_result <- anosim(unweight, group$group)
# Extract R and P values
R_value <- round(anosim_result $statistic, 2)
P_value <- round(anosim_result $signif, 3)
#添加P和R
pc3 <- pc2 + annotate("text", x = -0.1, y = 0.5, label = paste0("R = ", R_value, ", P = ", P_value), size = 5, fontface = "bold")
pc3


# 绘制y轴为PC2值的分组箱线图
pc4 <- ggplot(data,aes(x = factor(group, levels = group_order),y=A2))+
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
  scale_fill_manual(values = color, breaks = group_order)+
  geom_signif(
    comparisons = comparisons[which(sapply(comparisons, function(c) wilcox.test(data$A2[data$group == c[1]], data$A2[data$group == c[2]])$p.value < alpha_level))],
    map_signif_level = T, 
    test = t.test,
    y_position = c(0.3,0.4,0.35),
    tip_length = c(c(0,0),
                   c(0,0),
                   c(0,0)),
    size=0.4,color="black")
pc4
# 绘制y轴为PC1值的分组箱线图
pc5 <- ggplot(data,aes(x = factor(group, levels = group_order),y=A1))+
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
  scale_fill_manual(values = color, breaks = group_order)+
  geom_signif(
    comparisons = comparisons[which(sapply(comparisons, function(c) wilcox.test(data$A1[data$group == c[1]], data$A1[data$group == c[2]])$p.value < alpha_level))],
    map_signif_level = T,
    test = t.test,
    y_position = c(0.48,0.55,0.55),
    tip_length = c(c(0,0),
                   c(0,0),
                   c(0,0)),
    size=0.4,color="black")
pc5
# ggpubr::ggarrange()函数对图进行拼接
ggarrange(pc5, NULL, pc3, pc4, widths = c(7,2), heights = c(1.5,4), align = "hv")
ggsave("PCoA_unweight.tiff", height = 8, width = 8, dpi = 600)
