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
weight <- as.dist(weight)
unweight <- as.dist(unweight)

# NMDS分析
# NMDS排序分析——vegan包中的metaMDS函数
df_nmds <- metaMDS(bray, k = 2)

#结果查看——关注stress、points及species三个指标
summary(df_nmds)

#应力函数值（<=0.2合理）
df_nmds_stress <- df_nmds$stress
df_nmds_stress

#检查观测值非相似性与排序距离之间的关系——没有点分布在线段较远位置表示该数据可以使用NMDS分析
stressplot(df_nmds)

#提取作图数据
df_points <- as.data.frame(df_nmds$points)

#添加samp1es变量
df_points$sample <- row.names(df_points)

#修改列名
names(df_points)[1:2] <- c('NMDS1', 'NMDS2')
head(df_points)

#绘制散点图
p <- ggplot(df_points,aes(x=NMDS1, y=NMDS2))+#指定数据、X轴、Y轴
  geom_point(size=3)+#绘制点图并设定大小
  theme_bw()#主题
p

#将绘图数据和分组合并
df <- merge(df_points,group,by="sample")
head(df)

#使用ggplot2包绘图
p1<-ggplot(data=df,aes(x=NMDS1,y=NMDS2))+#指定数据、X轴、Y轴，颜色
  theme_bw()+#主题设置
  geom_point(aes(color = group), shape = 19, size=3)+#绘制点图并设定大小
  theme(panel.grid = element_blank())+
  geom_vline(xintercept = 0,lty="dashed", size = 1, color = 'grey50')+
  geom_hline(yintercept = 0,lty="dashed", size = 1, color = 'grey50')+#图中虚线
  geom_text(aes(label=sample, y=NMDS2+0.03,x=NMDS1+0.03,
                vjust=0, color = group),size=3.5, show.legend = F, alpha =0.5)+#添加数据点的标签
  stat_ellipse(data=df,
               geom = "polygon",level=0.95,
               linetype = 2,size=0.5,
               aes(fill=group),
               alpha=0.2)+
  scale_color_manual(values = color, breaks = group_order) +#点的颜色设置
  scale_fill_manual(values = color, breaks = group_order)+#椭圆颜色
  theme(axis.title.x=element_text(size=12),#修改X轴标题文本
        axis.title.y=element_text(size=12,angle=90),#修改y轴标题文本
        axis.text.y=element_text(size=10),#修改x轴刻度标签文本
        axis.text.x=element_text(size=10),#修改y轴刻度标签文本
        panel.grid=element_blank())+#隐藏网格线
  ggtitle(paste('Stress=',round(df_nmds_stress, 3)))#添加应力函数值
p1
ggsave("NMDS_bray.tiff", height = 8, width = 8, dpi = 600)

# 绘制y轴为PC2值的分组箱线图
p2 <- ggplot(df,aes(x=factor(group, levels = group_order),y=NMDS2))+
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
  scale_fill_manual(values=color, breaks = group_order)+
  geom_signif(
    comparisons = comparisons[which(sapply(comparisons, function(c) wilcox.test(df$NMDS2[df$group == c[1]], df$NMDS2[df$group == c[2]])$p.value < alpha_level))],
    map_signif_level = T, 
    test = t.test,
    y_position = c(0.3,0.4,0.35),
    tip_length = c(c(0,0),
                   c(0,0),
                   c(0,0)),
    size=0.4,color="black",
    step_increase = 0.08,)
p2
# 绘制y轴为PC1值的分组箱线图
p3 <- ggplot(df,aes(x=factor(group, levels = group_order),y=NMDS1))+
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
  scale_fill_manual(values=color, breaks = group_order)+
  geom_signif(
    comparisons = comparisons[which(sapply(comparisons, function(c) wilcox.test(df$NMDS1[df$group == c[1]], df$NMDS1[df$group == c[2]])$p.value < alpha_level))],
    map_signif_level = T,
    test = t.test,
    y_position = c(0.78,0.90,0.86,0.82),
    tip_length = c(c(0,0),
                   c(0,0),
                   c(0,0)),
    size=0.4,color="black")
p3
# ggpubr::ggarrange()函数对图进行拼接
ggarrange(p3, NULL, p1, p2, widths = c(5,2), heights = c(2,4), align = "hv")
ggsave("NMDS_bray.tiff", height = 8, width = 8, dpi = 500)


#NMDS排序分析——vegan包中的metaMDS函数
df_nmds <- metaMDS(jaccard, k = 2)
#结果查看——关注stress、points及species三个指标
summary(df_nmds)
#应力函数值（<=0.2合理）
df_nmds_stress <- df_nmds$stress
df_nmds_stress
#检查观测值非相似性与排序距离之间的关系——没有点分布在线段较远位置表示该数据可以使用NMDS分析
stressplot(df_nmds)
#提取作图数据
df_points <- as.data.frame(df_nmds$points)
#添加samp1es变量
df_points$sample <- row.names(df_points)
#修改列名
names(df_points)[1:2] <- c('NMDS1', 'NMDS2')
head(df_points)
#绘制散点图
p <- ggplot(df_points,aes(x=NMDS1, y=NMDS2))+#指定数据、X轴、Y轴
  geom_point(size=3)+#绘制点图并设定大小
  theme_bw()#主题
p
#将绘图数据和分组合并
df <- merge(df_points,group,by="sample")
head(df)
#使用ggplot2包绘图
p1<-ggplot(data=df,aes(x=NMDS1,y=NMDS2))+#指定数据、X轴、Y轴，颜色
  theme_bw()+#主题设置
  geom_point(aes(color = group), shape = 19, size=3)+#绘制点图并设定大小
  theme(panel.grid = element_blank())+
  geom_vline(xintercept = 0,lty="dashed", size = 1, color = 'grey50')+
  geom_hline(yintercept = 0,lty="dashed", size = 1, color = 'grey50')+#图中虚线
  geom_text(aes(label=sample, y=NMDS2+0.03,x=NMDS1+0.03,
                vjust=0, color = group),size=3.5, show.legend = F, alpha = 0.7)+#添加数据点的标签
  stat_ellipse(data=df,
               geom = "polygon",level=0.95,
               linetype = 2,size=0.5,
               aes(fill=group),
               alpha=0.2)+
  scale_color_manual(values = color, breaks = group_order) +#点的颜色设置
  scale_fill_manual(values = color, breaks = group_order)+#椭圆颜色
  theme(axis.title.x=element_text(size=12),#修改X轴标题文本
        axis.title.y=element_text(size=12,angle=90),#修改y轴标题文本
        axis.text.y=element_text(size=10),#修改x轴刻度标签文本
        axis.text.x=element_text(size=10),#修改y轴刻度标签文本
        panel.grid=element_blank())+#隐藏网格线
  ggtitle(paste('Stress=',round(df_nmds_stress, 3)))#添加应力函数值
p1
ggsave("NMDS_ jaccard.tiff", height = 8, width = 8, dpi = 600)

# 绘制y轴为PC2值的分组箱线图
p2 <- ggplot(df,aes(x=factor(group, levels = group_order),y=NMDS2))+
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
  scale_fill_manual(values=color, breaks = group_order)+
  geom_signif(
    comparisons = comparisons[which(sapply(comparisons, function(c) wilcox.test(df$NMDS2[df$group == c[1]], df$NMDS2[df$group == c[2]])$p.value < alpha_level))],
    map_signif_level = T, 
    test = t.test,
    y_position = c(0.3,0.4,0.35),
    tip_length = c(c(0,0),
                   c(0,0),
                   c(0,0)),
    size=0.4,color="black")
p2
# 绘制y轴为PC1值的分组箱线图
p3 <- ggplot(df,aes(x=group,y=NMDS1))+
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
  scale_fill_manual(values=color, breaks = group_order)+
  geom_signif(
    comparisons = comparisons[which(sapply(comparisons, function(c) wilcox.test(df$NMDS1[df$group == c[1]], df$NMDS1[df$group == c[2]])$p.value < alpha_level))],
    map_signif_level = T,
    test = t.test,
    y_position = c(0.78,0.92,0.85),
    tip_length = c(c(0,0),
                   c(0,0),
                   c(0,0)),
    size=0.4,color="black")
p3
# ggpubr::ggarrange()函数对图进行拼接
ggarrange(p3, NULL, p1, p2, widths = c(5,2), heights = c(2,4), align = "hv")
ggsave("NMDS_ jaccard.tiff", height = 8, width = 8, dpi = 500)


#weight_unifrac
#NMDS排序分析——vegan包中的metaMDS函数
df_nmds <- metaMDS(weight, k = 2)
#结果查看——关注stress、points及species三个指标
summary(df_nmds)
#应力函数值（<=0.2合理）
df_nmds_stress <- df_nmds$stress
df_nmds_stress
#检查观测值非相似性与排序距离之间的关系——没有点分布在线段较远位置表示该数据可以使用NMDS分析
stressplot(df_nmds)
#提取作图数据
df_points <- as.data.frame(df_nmds$points)
#添加samp1es变量
df_points$sample <- row.names(df_points)
#修改列名
names(df_points)[1:2] <- c('NMDS1', 'NMDS2')
head(df_points)
#绘制散点图
p <- ggplot(df_points,aes(x=NMDS1, y=NMDS2))+#指定数据、X轴、Y轴
  geom_point(size=3)+#绘制点图并设定大小
  theme_bw()#主题
p
#将绘图数据和分组合并
df <- merge(df_points,group,by="sample")
head(df)
#使用ggplot2包绘图
p1<-ggplot(data=df,aes(x=NMDS1,y=NMDS2))+#指定数据、X轴、Y轴，颜色
  theme_bw()+#主题设置
  geom_point(aes(color = group), shape = 19, size=3)+#绘制点图并设定大小
  theme(panel.grid = element_blank())+
  geom_vline(xintercept = 0,lty="dashed", size = 1, color = 'grey50')+
  geom_hline(yintercept = 0,lty="dashed", size = 1, color = 'grey50')+#图中虚线
  geom_text(aes(label=sample, y=NMDS2+0.03,x=NMDS1+0.03,
                vjust=0, color = group),size=3.5, show.legend = F, alpha=0.7)+#添加数据点的标签
  stat_ellipse(data=df,
               geom = "polygon",level=0.95,
               linetype = 2,size=0.5,
               aes(fill=group),
               alpha=0.2)+
  scale_color_manual(values = color, breaks = group_order) +#点的颜色设置
  scale_fill_manual(values = color, breaks = group_order)+#椭圆颜色
  theme(axis.title.x=element_text(size=12),#修改X轴标题文本
        axis.title.y=element_text(size=12,angle=90),#修改y轴标题文本
        axis.text.y=element_text(size=10),#修改x轴刻度标签文本
        axis.text.x=element_text(size=10),#修改y轴刻度标签文本
        panel.grid=element_blank())+#隐藏网格线
  ggtitle(paste('Stress=',round(df_nmds_stress, 3)))#添加应力函数值
p1
ggsave("NMDS_ weight_unifrac.tiff", height = 8, width = 8, dpi = 600)

# 绘制y轴为PC2值的分组箱线图
p2 <- ggplot(df,aes(x=factor(group, levels = group_order),y=NMDS2))+
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
  scale_fill_manual(values=color, breaks = group_order)+
  geom_signif(
    comparisons = comparisons[which(sapply(comparisons, function(c) wilcox.test(df$NMDS2[df$group == c[1]], df$NMDS2[df$group == c[2]])$p.value < alpha_level))],
    map_signif_level = T, 
    test = t.test,
    y_position = c(0.3,0.4,0.35),
    tip_length = c(c(0,0),
                   c(0,0),
                   c(0,0)),
    size=0.4,color="black")
p2
# 绘制y轴为PC1值的分组箱线图
p3 <- ggplot(df,aes(x=group,y=NMDS1))+
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
  scale_fill_manual(values=color, breaks = group_order)+
  geom_signif(
    comparisons = comparisons[which(sapply(comparisons, function(c) wilcox.test(df$NMDS1[df$group == c[1]], df$NMDS1[df$group == c[2]])$p.value < alpha_level))],
    map_signif_level = T,
    test = t.test,
    y_position = c(0.78,0.92,0.85),
    tip_length = c(c(0,0),
                   c(0,0),
                   c(0,0)),
    size=0.4,color="black")
p3
# ggpubr::ggarrange()函数对图进行拼接
ggarrange(p3, NULL, p1, p2, widths = c(5,2), heights = c(2,4), align = "hv")
ggsave("NMDS_ weight_unifrac.tiff", height = 8, width = 8, dpi = 500)


#unweight_unifrac
#NMDS排序分析——vegan包中的metaMDS函数
df_nmds <- metaMDS(unweight, k = 2)
#结果查看——关注stress、points及species三个指标
summary(df_nmds)
#应力函数值（<=0.2合理）
df_nmds_stress <- df_nmds$stress
df_nmds_stress
#检查观测值非相似性与排序距离之间的关系——没有点分布在线段较远位置表示该数据可以使用NMDS分析
stressplot(df_nmds)
#提取作图数据
df_points <- as.data.frame(df_nmds$points)
#添加samp1es变量
df_points$sample <- row.names(df_points)
#修改列名
names(df_points)[1:2] <- c('NMDS1', 'NMDS2')
head(df_points)
#绘制散点图
p <- ggplot(df_points,aes(x=NMDS1, y=NMDS2))+#指定数据、X轴、Y轴
  geom_point(size=3)+#绘制点图并设定大小
  theme_bw()#主题
p
#将绘图数据和分组合并
df <- merge(df_points,group,by="sample")
head(df)
#使用ggplot2包绘图
p1<-ggplot(data=df,aes(x=NMDS1,y=NMDS2))+#指定数据、X轴、Y轴，颜色
  theme_bw()+#主题设置
  geom_point(aes(color = group), shape = 19, size=3)+#绘制点图并设定大小
  theme(panel.grid = element_blank())+
  geom_vline(xintercept = 0,lty="dashed", size = 1, color = 'grey50')+
  geom_hline(yintercept = 0,lty="dashed", size = 1, color = 'grey50')+#图中虚线
  geom_text(aes(label=sample, y=NMDS2+0.03,x=NMDS1+0.03,
                vjust=0, color = group),size=3.5, show.legend = F, alpha = 0.7)+#添加数据点的标签
  stat_ellipse(data=df,
               geom = "polygon",level=0.95,
               linetype = 2,size=0.5,
               aes(fill=group),
               alpha=0.2)+
  scale_color_manual(values = color, breaks = group_order) +#点的颜色设置
  scale_fill_manual(values = color, breaks = group_order)+#椭圆颜色
  theme(axis.title.x=element_text(size=12),#修改X轴标题文本
        axis.title.y=element_text(size=12,angle=90),#修改y轴标题文本
        axis.text.y=element_text(size=10),#修改x轴刻度标签文本
        axis.text.x=element_text(size=10),#修改y轴刻度标签文本
        panel.grid=element_blank())+#隐藏网格线
  ggtitle(paste('Stress=',round(df_nmds_stress, 3)))#添加应力函数值
p1
ggsave("NMDS_ unweight_unifrac.tiff", height = 8, width = 8, dpi = 600)

# 绘制y轴为PC2值的分组箱线图
p2 <- ggplot(df,aes(x=factor(group,levels = group_order),y=NMDS2))+
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
  scale_fill_manual(values=color, break = group_order)+
  geom_signif(
    comparisons = comparisons[which(sapply(comparisons, function(c) wilcox.test(df$NMDS2[df$group == c[1]], df$NMDS2[df$group == c[2]])$p.value < alpha_level))],
    map_signif_level = T, 
    test = t.test,
    y_position = c(0.3,0.4,0.35),
    tip_length = c(c(0,0),
                   c(0,0),
                   c(0,0)),
    size=0.4,color="black")
p2
# 绘制y轴为PC1值的分组箱线图
p3 <- ggplot(df,aes(x=group,y=NMDS1))+
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
  scale_fill_manual(values=color, break = group_order)+
  geom_signif(
    comparisons = comparisons[which(sapply(comparisons, function(c) wilcox.test(df$NMDS1[df$group == c[1]], df$NMDS1[df$group == c[2]])$p.value < alpha_level))],
    map_signif_level = T,
    test = t.test,
    y_position = c(0.48,0.55,0.55),
    tip_length = c(c(0,0),
                   c(0,0),
                   c(0,0)),
    size=0.4,color="black")
p3
# ggpubr::ggarrange()函数对图进行拼接
ggarrange(p3, NULL, p1, p2, widths = c(5,2), heights = c(2,4), align = "hv")
ggsave("NMDS_ unweight_unifrac.tiff", height = 8, width = 8, dpi = 500)

