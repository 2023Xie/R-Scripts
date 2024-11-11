# 加载R包
library(ggplot2) #用于 ggplot2作图
library(ggsignif)
library(ggpubr)
library(tidyverse)
library(vegan)

# 导入数据
otu_table <- read.table(" otu_Flattening.csv", header = T, row.names = 1, sep = ",", quote = "", stringsAsFactors = F)
name <- read.table(" otu_Flattening.csv", header = F, row.names = 1, sep = ",", quote = "", stringsAsFactors = F)
names(otu_table) <- name[1, ]
rm(name)

# 转置
otu <- t(otu_table)

# 导入分组文件
group <- read.table('group.txt', header=T, sep="\t", quote = "",comment.char="",stringsAsFactors = FALSE)
names(group) <- c("sample", "group")
group$group <- as.factor(group$group)

# 分组顺序和颜色
group_order <- c("CON", "HFD", "AKK", "HFU", "LFU", "AKK+HFU", "AKK+LFU")
color=c("#264653","#287271","#2a9d8c","#8ab07d","#e9c46b","#f3a261","#e66f51")

# 计算Alpha多样性函数
Alpha_diversity_index <- function(x, tree = NULL, base = exp(1)) {
  est <- estimateR(x)
  Obs <-  est[1, ]
  Shannon <- diversity(x, index = 'shannon', base = base)
  Simpson <- diversity(x, index = 'simpson')
  Pielou <- Shannon / log(Obs, base)
  goods_coverage <- 1 - rowSums(x == 1) / rowSums(x)
  
  result <- rbind(est, Shannon, Simpson,
                  Pielou, goods_coverage)
  if (!is.null(tree)) {
    Pd <- pd(x, tree, include.root = FALSE)[1]
    Pd <- t(Pd)
    result <- rbind(result, Pd)
  }
  result <- as.data.frame(t(result))
  return(result)
}

# 计算alpha多样性
alpha_diversity <- Alpha_diversity_index(otu)
alpha_diversity[c(1:3), ]

# 比较信息
alpha_level <- 0.05
comparisons <- list(c("CON", "HFD"), c("CON", "AKK"), c("CON", "LFU"), c("CON", "HFU"), c("CON", "AKK+HFU"), c("CON", "AKK+LFU"), c("HFD", "HFU"), c("HFD", "LFU"), c("HFD", "AKK+HFU"), c("HFD", "AKK+LFU"), c("HFD", "AKK"), c("AKK", "AKK+LFU"), c("AKK", "AKK+HFU"), c("AKK", "LFU"), c("AKK", "HFU"), c("LFU", "HFU"), c("LFU", "AKK+HFU"), c("LFU", "AKK+LFU"), c("HFU", "AKK+HFU"), c("HFU", "AKK+LFU"), c("AKK+HFU", "AKK+LFU"))

# 添加sample列
alpha_diversity$sample <- rownames(alpha_diversity)

# 绘图文件
plot_data <- merge(group, alpha_diversity, by = "sample")

# 设置主题
custom_theme <- theme(panel.background = element_blank(),
                      axis.line.x = element_line(colour = "black"),
                      axis.line.y = element_line(colour = "black"),
                      axis.title.x = element_blank(),
                      axis.title.y = element_text(size = 30),
                      axis.text.x = element_text(size = 21),
                      axis.text.y = element_text(size = 21),
)

# SHANNON
shannon <- ggplot(data = plot_data, aes(x = group, y = Shannon, fill = factor(group, levels = group_order)))+ 
  stat_boxplot(geom ='errorbar', width = 0.1, size = 0.5)+
  geom_boxplot(fill = color, color = "black")+ 
  geom_signif(
    comparisons = comparisons[which(sapply(comparisons, function(c) wilcox.test(plot_data$Shannon[plot_data$group == c[1]], plot_data$Shannon[plot_data$group == c[2]])$p.value < alpha_level))],
    map_signif_level = TRUE, textsize = 9,
    step_increase = 0.15,
  ) +
  ylim(NA, 6)+
  xlab("")+
  scale_x_discrete(limits = group_order) +
  ylab("Shannon")+
  custom_theme
shannon
ggsave("shannon.tiff", plot=shannon, width =6, height =6, dpi = 600)

#ACE
ace <- ggplot(data = plot_data, aes(x = group, y = S.ACE, fill = factor(group, levels = group_order)))+ 
  stat_boxplot(geom ='errorbar', width = 0.1, size = 0.5)+
  geom_boxplot(fill = color, color = "black")+ 
  geom_signif(
    comparisons = comparisons[which(sapply(comparisons, function(c) wilcox.test(plot_data$S.ACE[plot_data$group == c[1]], plot_data$S.ACE[plot_data$group == c[2]])$p.value < alpha_level))],
    map_signif_level = TRUE, textsize = 9,
    step_increase = 0.15,
  ) +
  ylim(NA, 320)+
  xlab("")+
  scale_x_discrete(limits = group_order) +
  ylab("ACE")+
  custom_theme
ace
ggsave("ACE.tiff", plot=ace, width =6, height =6, dpi = 600)

#CHAO1
chao1 <- ggplot(data = plot_data, aes(x = group, y = S.chao1, fill = factor(group, levels = group_order)))+ 
  stat_boxplot(geom ='errorbar', width = 0.1, size = 0.5)+
  geom_boxplot(fill = color, color = "black")+ 
  geom_signif(
    comparisons = comparisons[which(sapply(comparisons, function(c) wilcox.test(plot_data$S.chao1[plot_data$group == c[1]], plot_data$S.chao1[plot_data$group == c[2]])$p.value < alpha_level))],
    map_signif_level = TRUE, textsize = 9,
    step_increase = 0.15,
  ) +
  ylim(NA, 320)+
  xlab("")+
  scale_x_discrete(limits = group_order) +
  ylab("Chao1")+
  custom_theme
chao1
ggsave("chao1.tiff", plot=chao1, width =6, height =6, dpi = 600)

#SIMPSON
Simpson <- ggplot(data = plot_data, aes(x = group, y = Simpson, fill = factor(group, levels = group_order)))+ 
  stat_boxplot(geom ='errorbar', width = 0.1, size = 0.5)+
  geom_boxplot(fill = color, color = "black")+ 
  geom_signif(
    comparisons = comparisons[which(sapply(comparisons, function(c) wilcox.test(plot_data$Simpson[plot_data$group == c[1]], plot_data$Simpson[plot_data$group == c[2]])$p.value < alpha_level))],
    map_signif_level = TRUE, textsize = 9,
    step_increase = 0.15,
  ) +
  ylim(NA, 1)+
  xlab("")+
  scale_x_discrete(limits = group_order) +
  ylab("Simpson")+
  custom_theme
Simpson
ggsave("Simpson.tiff", plot=Simpson, width =6, height =6, dpi = 600)

# 稀释曲线
alpha_index <- function(x, method = 'richness', tree = NULL, base = exp(1)) {
  if (method == 'richness') result <- rowSums(x > 0)    #丰富度指数
  else if (method == 'chao1') result <- estimateR(x)[2, ]    #Chao1 指数
  else if (method == 'ace') result <- estimateR(x)[4, ]    #ACE 指数
  else if (method == 'shannon') result <- diversity(x, index = 'shannon', base = base)    #Shannon 指数
  else if (method == 'simpson') result <- diversity(x, index = 'simpson')    #Gini-Simpson 指数
  else if (method == 'pielou') result <- diversity(x, index = 'shannon', base = base) / log(estimateR(x)[1, ], base)    #Pielou 均匀度
  else if (method == 'gc') result <- 1 - rowSums(x == 1) / rowSums(x)    #goods_coverage
  else if (method == 'pd' & !is.null(tree)) {    #PD_whole_tree
    pd <- pd(x, tree, include.root = FALSE)
    result <- pd[ ,1]
    names(result) <- rownames(pd)
  }
  result
}


alpha_curves <- function(x, step, method = 'richness', rare = NULL, tree = NULL, base = exp(1)) {
  x_nrow <- nrow(x)
  if (is.null(rare)) rare <- rowSums(x) else rare <- rep(rare, x_nrow)
  alpha_rare <- list()
  
  for (i in 1:x_nrow) {
    step_num <- seq(0, rare[i], step)
    if (max(step_num) < rare[i]) step_num <- c(step_num, rare[i])
    
    alpha_rare_i <- NULL
    for (step_num_n in step_num) alpha_rare_i <- c(alpha_rare_i, alpha_index(x = rrarefy(x[i, ], step_num_n), method = method, tree = tree, base = base))
    names(alpha_rare_i) <- step_num
    alpha_rare <- c(alpha_rare, list(alpha_rare_i))
  }
  names(alpha_rare) <- rownames(x)
  alpha_rare
}


# richness稀释曲线

richness_curves <- alpha_curves(otu, step = 150, method = 'richness')

plot_richness <- data.frame()

for (i in names(richness_curves)) {
  
  richness_curves_i <- (richness_curves[[i]])
  
  richness_curves_i <- data.frame(depth = names(richness_curves_i), richness = richness_curves_i, sample = i, stringsAsFactors = FALSE)
  
  plot_richness <- rbind(plot_richness, richness_curves_i)
  
}

rownames(plot_richness) <- NULL

plot_richness$depth <- as.numeric(plot_richness$depth)

plot_richness$richness <- as.numeric(plot_richness$richness)

plot_richness <- plot_richness %>% mutate(depth = ifelse(depth == 0, 1, depth))

plot_richness <- plot_richness %>% mutate(richness = ifelse(richness == 0, 1, richness))

plot_richness <- merge(plot_richness, group, by = "sample")

S1 <- ggplot(plot_richness, aes(depth, richness, color = factor(group, levels = group_order), group = `sample`))+
  geom_smooth(se = FALSE, method = "lm", formula = y~log(x), linewidth = 0.5, alpha =0.7)+
  scale_color_manual(values = color)+ # 设置分组颜色
  labs(color = "group")
S1
ggsave("richness_rarefaction_samples.tiff", height = 5, width = 8, dpi = 300)

S2 <- ggplot(plot_richness, aes(depth, richness, color = factor(group, levels = group_order), group = `group`))+
  geom_smooth(se = FALSE, method = "lm", formula = y~log(x) , linewidth = 0.5, alpha =0.7)+
  scale_color_manual(values = color)+ # 设置分组颜色
  labs(color = "group")
S2
ggsave("richness_rarefaction_groups.tiff", height = 5, width = 8, dpi = 300)

# shannon稀释曲线

shannon_curves <- alpha_curves(otu, step = 150, method = 'shannon')

plot_shannon <- data.frame()

for (i in names(shannon_curves)) {
  
  shannon_curves_i <- (shannon_curves[[i]])
  
  shannon_curves_i <- data.frame(depth = names(shannon_curves_i), shannon = shannon_curves_i, sample = i, stringsAsFactors = FALSE)
  
  plot_shannon <- rbind(plot_shannon, shannon_curves_i)
  
}

rownames(plot_shannon) <- NULL

plot_shannon$depth <- as.numeric(plot_shannon$depth)

plot_shannon$shannon <- as.numeric(plot_shannon$shannon)

plot_shannon <- plot_shannon %>% mutate(depth = ifelse(depth == 0, 1, depth))

plot_shannon <- plot_shannon %>% mutate(shannon = ifelse(shannon == 0, 1, shannon))

plot_shannon <- merge(plot_shannon, group, by = "sample")

S1 <- ggplot(plot_shannon, aes(depth, shannon, color = factor(group, levels = group_order), group = `sample`))+
  geom_smooth(se = FALSE, method = "lm", formula = y~log(x), linewidth = 0.5, alpha =0.7)+
  scale_color_manual(values = color)+ # 设置分组颜色
  labs(color = "group")
S1
ggsave("shannon_rarefaction_samples.tiff", height = 5, width = 8, dpi = 300)


S2 <- ggplot(plot_shannon, aes(depth, shannon, color = factor(group, levels = group_order), group = `group`))+
  geom_smooth(se = FALSE, method = "lm", formula = y~log(x) , linewidth = 0.5, alpha =0.7)+
  scale_color_manual(values = color)+ # 设置分组颜色
  labs(color = "group")
S2
ggsave("shannon_rarefaction_groups.tiff", height = 5, width = 8, dpi = 300)


# simpson稀释曲线

simpson_curves <- alpha_curves(otu, step = 150, method = 'simpson')

plot_simpson <- data.frame()

for (i in names(simpson_curves)) {
  
  simpson_curves_i <- (simpson_curves[[i]])
  
  simpson_curves_i <- data.frame(depth = names(simpson_curves_i), simpson = simpson_curves_i, sample = i, stringsAsFactors = FALSE)
  
  plot_simpson <- rbind(plot_simpson, simpson_curves_i)
  
}

rownames(plot_simpson) <- NULL

plot_simpson$depth <- as.numeric(plot_simpson$depth)

plot_simpson$simpson <- as.numeric(plot_simpson$simpson)

plot_simpson <- plot_simpson %>% mutate(depth = ifelse(depth == 0, 1, depth))

plot_simpson <- plot_simpson %>% mutate(simpson = ifelse(simpson == 0, 1, simpson))

plot_simpson <- merge(plot_simpson, group, by = "sample")

S1 <- ggplot(plot_simpson, aes(depth, simpson, color = factor(group, levels = group_order), group = `sample`))+
  geom_smooth(se = FALSE, method = "lm", formula = y~log(x), linewidth = 0.5, alpha =0.7)+
  scale_color_manual(values = color)+ # 设置分组颜色
  labs(color = "group")
S1
ggsave("simpson_rarefaction_samples.tiff", height = 5, width = 8, dpi = 300)


S2 <- ggplot(plot_simpson, aes(depth, simpson, color = factor(group, levels = group_order), group = `group`))+
  geom_smooth(se = FALSE, method = "lm", formula = y~log(x) , linewidth = 0.5, alpha =0.7)+
  scale_color_manual(values = color)+ # 设置分组颜色
  labs(color = "group")
S2
ggsave("simpson_rarefaction_groups.tiff", height = 5, width = 8, dpi = 300)



# chao1稀释曲线

chao1_curves <- alpha_curves(otu, step = 150, method = 'chao1')

plot_chao1 <- data.frame()

for (i in names(chao1_curves)) {
  
  chao1_curves_i <- (chao1_curves[[i]])
  
  chao1_curves_i <- data.frame(depth = names(chao1_curves_i), chao1 = chao1_curves_i, sample = i, stringsAsFactors = FALSE)
  
  plot_chao1 <- rbind(plot_chao1, chao1_curves_i)
  
}

rownames(plot_chao1) <- NULL

plot_chao1$depth <- as.numeric(plot_chao1$depth)

plot_chao1$chao1 <- as.numeric(plot_chao1$chao1)

plot_chao1 <- plot_chao1 %>% mutate(depth = ifelse(depth == 0, 1, depth))

plot_chao1 <- plot_chao1 %>% mutate(chao1 = ifelse(chao1 == 0, 1, chao1))

plot_chao1 <- merge(plot_chao1, group, by = "sample")

S1 <- ggplot(plot_chao1, aes(depth, chao1, color = factor(group, levels = group_order), group = `sample`))+
  geom_smooth(se = FALSE, method = "lm", formula = y~log(x), linewidth = 0.5, alpha =0.7)+
  scale_color_manual(values = color)+ # 设置分组颜色
  labs(color = "group")
S1
ggsave("chao1_rarefaction_samples.tiff", height = 5, width = 8, dpi = 300)


S2 <- ggplot(plot_chao1, aes(depth, chao1, color = factor(group, levels = group_order), group = `group`))+
  geom_smooth(se = FALSE, method = "lm", formula = y~log(x) , linewidth = 0.5, alpha =0.7)+
  scale_color_manual(values = color)+ # 设置分组颜色
  labs(color = "group")
S2
ggsave("chao1_rarefaction_groups.tiff", height = 5, width = 8, dpi = 300)


# ace稀释曲线

ace_curves <- alpha_curves(otu, step = 150, method = 'ace')

plot_ace <- data.frame()

for (i in names(ace_curves)) {
  
  ace_curves_i <- (ace_curves[[i]])
  
  ace_curves_i <- data.frame(depth = names(ace_curves_i), ace = ace_curves_i, sample = i, stringsAsFactors = FALSE)
  
  plot_ace <- rbind(plot_ace, ace_curves_i)
  
}

rownames(plot_ace) <- NULL

plot_ace$depth <- as.numeric(plot_ace$depth)

plot_ace$ace <- as.numeric(plot_ace$ace)

plot_ace <- plot_ace %>% mutate(depth = ifelse(depth == 0, 1, depth))

plot_ace <- plot_ace %>% mutate(ace = ifelse(ace == 0, 1, ace))

plot_ace <- merge(plot_ace, group, by = "sample")

S1 <- ggplot(plot_ace, aes(depth, ace, color = factor(group, levels = group_order), group = `sample`))+
  geom_smooth(se = FALSE, method = "lm", formula = y~log(x), linewidth = 0.5, alpha =0.7)+
  scale_color_manual(values = color)+ # 设置分组颜色
  labs(color = "group")
S1
ggsave("ace_rarefaction_samples.tiff", height = 5, width = 8, dpi = 300)


S2 <- ggplot(plot_ace, aes(depth, ace, color = factor(group, levels = group_order), group = `group`))+
  geom_smooth(se = FALSE, method = "lm", formula = y~log(x) , linewidth = 0.5, alpha =0.7)+
  scale_color_manual(values = color)+ # 设置分组颜色
  labs(color = "group")
S2
ggsave("ace_rarefaction_groups.tiff", height = 5, width = 8, dpi = 300)


