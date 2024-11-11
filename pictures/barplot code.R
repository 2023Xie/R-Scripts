# 导入R包
library(tidyverse)
library(ggsci)

## 物种组成数据
spetax <- read.table("feature-table.txt", header<-T, row.names=NULL, sep="\t", comment.char="",stringsAsFactors = F,quote = "", skip = 1)
spetax1 <- read.table("feature-table.txt", header=F, row.names=NULL, sep="\t", comment.char="",stringsAsFactors = F,quote = "", skip = 1)
names(spetax) <- spetax1[1,]
rm(spetax1)

taxnommy <- read.table("taxonomy.tsv", header=T, row.names=NULL, sep="\t", comment.char="",stringsAsFactors = F,quote = "")
taxnommy <- taxnommy[1:2]
names(taxnommy) <- c("#OTU ID", "Taxon")

otu_table <- merge(spetax, taxnommy, by = "#OTU ID")

domain <- c()
delete_d <- c()
phylum <- c()
delete_p <-c()
class <- c()
delete_c <- c()
order <- c()
delete_o <- c()
family <- c()
delete_f <- c()
genus <- c()
delete_g <- c()
species <- c()
delete_s <- c()


count = 1
for(i in otu_table[,58]){
        tmp <- unlist(strsplit(as.character(i), split="; |__"))
        domain <- c(domain, tmp[2])
        phylum <- c(phylum, tmp[4])
        class <- c(class, tmp[6])
        order <- c(order, tmp[8])
        family <- c(family, tmp[10])
        genus <- c(genus, tmp[12])
        species <- c(species, tmp[14])
        
        if(is.na(tmp[2]) | tmp[2] == ""){
                delete_d <- c(delete_d, count)
        }
        if(is.na(tmp[4]) | tmp[4] == ""){
                delete_p <- c(delete_p, count)
        }
        if(is.na(tmp[6]) | tmp[6] == ""){
                delete_c <- c(delete_c, count)
        }
        if(is.na(tmp[8]) | tmp[8] == ""){
                delete_o <- c(delete_o, count)
        }
        if(is.na(tmp[10]) | tmp[10] == ""){
                delete_f <- c(delete_f, count)
        }
        if(is.na(tmp[12]) | tmp[12] == ""){
                delete_g <- c(delete_g, count)
        }
        if(is.na(tmp[14]) | tmp[14] == ""){
                delete_s <- c(delete_s, count)
        }
        count = count + 1
}
otu_table$Domain <- domain
otu_table$Phylum <- phylum
otu_table$Class <- class
otu_table$Order <- order
otu_table$Family <- family
otu_table$Genus <- genus
otu_table$species <- species

for (i in 59:65){
        for (j in 1:nrow(otu_table)){
               if (otu_table[j, i] == "" | is.na(otu_table[j, i])) {
                        if (endsWith(as.character(otu_table[j, i-1]), "_unclassified") | endsWith(as.character(otu_table[j, i-1]), "_unclassified")){
                                  otu_table[j, i] <- otu_table[j, i-1]
                        }
                        if (!endsWith(as.character(otu_table[j, i-1]), "_unclassified") & !endsWith(as.character(otu_table[j, i-1]), "_unclassified")){
                                  otu_table[j, i] <- paste(otu_table[j, i-1], "_unclassified")
                        }
               }
               if (otu_table[j, i] == "uncultured"){
                 otu_table[j, i] <- paste(otu_table[j, i-1], "_unclassified")
               }
        
        }
}

# 颜色
color<-c("#0E606B","#1597A5","#FEB3AE")
color<-c("#1597A5","#FFC24B","#FEB3AE")
color<-c("#264653","#287271","#2a9d8c","#8ab07d")
color<-c("#264653","#287271","#2a9d8c","#8ab07d","#e9c46b")
color<-c("#264653","#287271","#2a9d8c","#8ab07d","#e9c46b","#f3a261","#e66f51")

col=pal_d3("category20")(20)
col2 = pal_d3("category20",alpha = 0.5)(20)
mypal=c(col,col2[-8])

# 分组顺序
group_order <- c("CON", "T2DM", "MET")
group_order <- c("CON", "HFD", "LFU", "HFU")
group_order <- c("CON", "HFD", "AKK", "AKK+LFU", "AKK+HFU")
group_order <- c("CON", "HFD", "AKK", "HFU", "LFU", "AKK+HFU", "AKK+LFU")

# 比较组
comparisons <- list(c("CON", "T2DM"), c("CON", "MET"), c("T2DM", "MET"))
comparisons <- list(c("CON", "HFU"), c("CON", "LFU"), c("HFU", "LFU"))
comparisons <- list(c("CON", "HFD"), c("CON", "AKK"), c("CON", "LFU"), c("CON", "HFU"), c("CON", "AKK+HFU"), c("CON", "AKK+LFU"), c("HFD", "HFU"), c("HFD", "LFU"), c("HFD", "AKK+HFU"), c("HFD", "AKK+LFU"), c("HFD", "AKK"), c("AKK", "AKK+LFU"), c("AKK", "AKK+HFU"), c("AKK", "LFU"), c("AKK", "HFU"), c("LFU", "HFU"), c("LFU", "AKK+HFU"), c("LFU", "AKK+LFU"), c("HFU", "AKK+HFU"), c("HFU", "AKK+LFU"), c("AKK+HFU", "AKK+LFU"))

spe <- otu_table[2:57]
tax <- otu_table[, c(59:65)]

# 导入分组数据
group <- read.table("group.txt", header=T, row.names=NULL, sep="\t", comment.char="",stringsAsFactors = F,quote = "")
names(group) <- c("sample", "group")
group$group <- as.factor(group$group)
group$group <- factor(group$group, levels = group_order)

phy <- spe %>%
  group_by(tax$Phylum) %>% # 使用tax中的门水平进行分类
  summarise_all(sum) %>%
  rename(Phylum = `tax$Phylum`) %>%
  gather(key="sample",value = "abun",-Phylum) %>% # 数据形式转换：“宽”转“长”
  left_join(group,by=c("sample"="sample")) %>%
  select(group,Phylum,abun) %>%
  group_by(group,Phylum) %>% # 求均值
  summarise_all(mean)

gen <- spe %>%
  group_by(tax$Genus) %>% # 使用tax中的门水平进行分类
  summarise_all(sum) %>%
  rename(Genus = `tax$Genus`) %>%
  gather(key="sample",value = "abun",-Genus) %>% # 数据形式转换：“宽”转“长”
  left_join(group,by=c("sample"="sample")) %>%
  select(group,Genus,abun) %>%
  group_by(group,Genus) %>% # 求均值
  summarise_all(mean)

pdf("Phyum_stack.pdf",width = 12,height = 10,family="Times")
ggplot()+
  geom_bar(data=phy,
           aes(x=factor(group, levels = group_order),
               weight=abun,
               fill=reorder(Phylum,-abun)),
           position = "fill", # ggplot2自行计算相对丰度
           width=0.6,
           )+
  scale_fill_manual(values = mypal[-18])+ # 颜色与绝对丰度堆叠柱形图保持一致
  theme_bw()+
  scale_y_continuous(#expand=c(0,0), #设置横坐标轴紧挨柱形图
    name = "Relative abundance (%)",
    limits = c(0,1),
    breaks = seq(0,1,0.25),
    labels = paste(seq(0,100,25),"%")
  )+
  guides(fill=guide_legend(title = "Phylum",ncol = 1))+
  labs(x=NULL)+
  theme(legend.position="right",
        axis.title = element_text(face = "bold", 
                                  size = 30,colour = "black"))+
  theme(axis.text = element_text(face = "bold", 
                                 size = 21,color="black"),
        strip.text.x = element_text(face = "bold", 
                                    size =21,color="black", angle = 45 ))+
  theme(panel.grid=element_blank())+
  theme(legend.title = element_text(face = "bold", 
                                    size =27,color="black"),
        legend.text = element_text(face = "bold", 
                                   size =24,color="black"))
dev.off()

# 计算每个属的总丰度
gen_total <- gen %>%
  group_by(Genus) %>%
  summarise(total_abun = sum(abun))

# 对总丰度进行排序，并选择前二十的属
top_20_Genus <- gen_total %>%
  arrange(desc(total_abun)) %>%
  head(20) %>%
  pull(Genus)

# 将数据集中的属分为前二十和其他
gen$Genus <- as.character(gen$Genus)
gen$Genus <- ifelse(gen$Genus %in% top_20_Genus, gen$Genus, "Others")

gen$Genus <- as.factor(gen$Genus)
gen$group <- as.factor(gen$group)

genu <- gen %>%
  group_by(Genus, group) %>%
  summarise(abun_sum = sum(abun))

pdf("genus_stack.pdf",width = 12,height = 10,family="Times")
ggplot()+
  geom_bar(data=genu,
           aes(x=factor(group, levels = group_order),
               weight=abun_sum,
               fill=reorder(Genus,-abun_sum)),
           position = "fill", # ggplot2会自行计算相对丰度，无需提前计算。
           width=0.6)+
  scale_fill_manual(values = mypal[-18])+ # 颜色与绝对丰度堆叠柱形图保持一致
  theme_bw()+
  scale_y_continuous(#expand=c(0,0), #设置横坐标轴紧挨柱形图
    name = "Relative abundance (%)",
    limits = c(0,1),
    breaks = seq(0,1,0.25),
    labels = paste(seq(0,100,25),"%")
  )+
  guides(fill=guide_legend(title = "Genus",ncol = 1))+
  labs(x=NULL)+
  theme(legend.position="right",
        axis.title = element_text(face = "bold", 
                                  size = 30,colour = "black"))+
  theme(axis.text = element_text(face = "bold", 
                                 size = 21,color="black"),
        strip.text.x = element_text(face = "bold", 
                                    size =21,color="black"))+
  theme(panel.grid=element_blank())+
  theme(legend.title = element_text(face = "bold", 
                                    size =27,color="black"),
        legend.text = element_text(face = "bold", 
                                   size =24,color="black"))
dev.off()

