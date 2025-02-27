#Loading packages
library(biomformat)
library(BiocManager)
library(file2meco)
library(MicrobiomeStat)
library(WGCNA)
library(ggtree)
library(metagenomeSeq)
library(ALDEx2)
library(ANCOMBC)
library(microeco)
library(ape)
library(plyr)
library(magrittr)
library(tidygraph)
library(ggcor)
library(ggplot2)
library(ggpubr)
library(dplyr)
library(vegan)
library(tidyverse)
library(ggplot2)
library(tidyr)
library(dplyr)
library(paletteer)
library(colorspace)
# Load the required packages
library(lme4)
library(lmerTest)
library(multcomp)

# Developmental version
devtools::install_github("EdwinTh/dutchmasters")

bet = paletteer_d("dichromat::BluetoOrange_10")

bet = paletteer_d("rcartocolor::Magenta")

bet = paletteer_d("calecopal::chaparral3")

install.packages("ggsignif")
library(ggsignif)


##Importing microbiome data into R and converting this data into microeco package
##Genotype X Drought_Stage -- Phylum abundance plots

abund_phylum = read.csv("gen-drought-phylum.csv")

# Reshape the data to long format
abund_d_phylum <- abund_phylum %>%
  pivot_longer(
    cols = -c(Drought_Stage, Genotype),    # Select columns after the first four factors
    names_to = "Phylum",
    values_to = "Abundance"
  )

total_abundance_phylum <- abund_d_phylum %>%
  group_by(Phylum) %>%
  summarize(Total_Abundance = sum(Abundance, na.rm = TRUE)) %>%
  arrange(desc(Total_Abundance))

# Select the top 10 microbial groups
top_15_microbes_phylum <- total_abundance_phylum %>%
  top_n(15, Total_Abundance) %>%
  pull(Phylum)

# Filter the data to keep only the top 20 microbial groups
abund_d_phylum <- abund_d_phylum %>%
  filter(Phylum %in% top_15_microbes_phylum)

# Reorder Microbial_Group factor levels based on total abundance
abund_d_phylum$Phylum <- factor(abund_d_phylum$Phylum, levels = top_15_microbes_phylum)


# Reorder Drought_Stage factor levels
abund_d_phylum$Drought_Stage <- factor(abund_d_phylum$Drought_Stage, levels = c("V2-Stage", "V4-Stage", "R1-Stage", "R4-Stage"))

farrowAndBall_palette <- c(
  "#4d5b6a" #Stiffkey Blue
  ,"#6a90b4" #Cook's Blue
  ,"#599ec4" #ST Giles Blue
  ,"#a1c5c8" #Blue Ground
  ,"#7997a1" #Stone Blue
  ,"#427e83" #Vardo
  ,"#84b59c" #Arsenic
  ,"#919f70" #Yeabridge Green
  ,"#686a47" #Bancha
  ,"#c8bd83" #Churlish Green
  ,"#cb9e59" #India Yellow
  ,"#ecc363" #Babouche
  ,"#c57b67" #Red Earth
  ,"#d65f3d" #Charlotte's Locks
  ,"#a04344" #Incarnadine
  ,"#bf7a8f" #Rangwali
  ,"#8d8089" #Brassica
  ,"#50414c" #Pelt
  ,"#e5e0db" #Strong White
  ,"#444546" #Off-Black
)

top15_palette <- c("#290AD8", "#264DFF", 
                   "#3FA0FF", "#AAF7FF", 
                   "#B2FFB2", "#FFFFBF", 
                   "#FFE099", "#FFAD72", 
                   "#F76D5E", "#D82632", 
                   "#A50021", '#F3B79C',
                   '#D64D72', '#841859',
                   '#312A56', '#6B6100',
                   '#004F78', '#0096B5',
                   "#427e83", "#686a47")

alpha_pallete = c("dodgerblue", "#f54260")

# Inspect the transformed data
head(abund)

# Plot the data
abundnace_plot_phylum = ggplot(abund_d_phylum, aes(x = Drought_Stage, y = Abundance, fill = Phylum)) +
  geom_bar(stat = "identity", position = "stack") +
  facet_wrap(~ Genotype) +
  theme_minimal() +
  labs(title = "Phylum",
       x = "Drought Stage",
       y = "Abundance",
       fill = "Phylum") +
  theme(axis.text.x = element_text(size = 14, angle = 45, hjust = 1),  # Adjust x-axis text size
        axis.text.y = element_text(size = 14),  # Adjust y-axis text size
        axis.title.x = element_text(size = 14), # Adjust x-axis label size
        axis.title.y = element_text(size = 14), # Adjust y-axis label size
        strip.text = element_text(size = 14),   # Adjust facet labels text size
        legend.title = element_text(size = 14), # Adjust legend title text size
        legend.text = element_text(size = 12),  # Adjust legend text size
        plot.title = element_text(size = 14, hjust = 0.5)) +  # Adjust title size and center it
  scale_fill_manual(values = top15_palette)  # Apply the custom color palette

abundnace_plot_phylum

ggsave("abundance_plot_phylum.jpeg", plot = abundnace_plot_phylum, 
       width = 8, height = 6, units = "in", 
       dpi = 1000)

##########################################################################################
##########Gentotype X Stages  - Class ##########################

abund_class = read.csv("gen-drought-class.csv")

# Reshape the data to long format
abund_d_class <- abund_class %>%
  pivot_longer(
    cols = -c(Drought_Stage, Genotype),    # Select columns after the first four factors
    names_to = "Class",
    values_to = "Abundance"
  )

total_abundance_class <- abund_d_class %>%
  group_by(Class) %>%
  summarize(Total_Abundance = sum(Abundance, na.rm = TRUE)) %>%
  arrange(desc(Total_Abundance))

# Select the top 10 microbial groups
top_15_microbes_class <- total_abundance_class %>%
  top_n(15, Total_Abundance) %>%
  pull(Class)

# Filter the data to keep only the top 20 microbial groups
abund_d_class<- abund_d_class %>%
  filter(Class %in% top_15_microbes_class)

# Reorder Microbial_Group factor levels based on total abundance
abund_d_class$Class <- factor(abund_d_class$Class, levels = top_15_microbes_class)


# Reorder Drought_Stage factor levels
abund_d_class$Drought_Stage <- factor(abund_d_class$Drought_Stage, levels = c("V2-Stage", "V4-Stage", "R1-Stage", "R4-Stage"))

# Plot the data
abundnace_plot_class = ggplot(abund_d_class, aes(x = Drought_Stage, y = Abundance, fill = Class)) +
  geom_bar(stat = "identity", position = "stack") +
  facet_wrap(~ Genotype) +
  theme_minimal() +
  labs(title = "Class",
       x = "Drought Stage",
       y = "Abundance",
       fill = "Class") +
  theme(axis.text.x = element_text(size = 14, angle = 45, hjust = 1),  # Adjust x-axis text size
        axis.text.y = element_text(size = 14),  # Adjust y-axis text size
        axis.title.x = element_text(size = 14), # Adjust x-axis label size
        axis.title.y = element_text(size = 14), # Adjust y-axis label size
        strip.text = element_text(size = 14),   # Adjust facet labels text size
        legend.title = element_text(size = 14), # Adjust legend title text size
        legend.text = element_text(size = 12),  # Adjust legend text size
        plot.title = element_text(size = 14, hjust = 0.5)) +  # Adjust title size and center it
  scale_fill_manual(values = top15_palette)  # Apply the custom color palette

abundnace_plot_class

ggsave("abundance_plot_class.jpeg", plot = abundnace_plot_class, 
       width = 12, height = 10, units = "in", 
       dpi = 1000)

#########################################################################################
#######################Genotype X stages  - order ####################################

abund_order = read.csv("gen-drought-order.csv")

# Reshape the data to long format
abund_d_order <- abund_order %>%
  pivot_longer(
    cols = -c(Drought_Stage, Genotype),    # Select columns after the first four factors
    names_to = "Order",
    values_to = "Abundance"
  )

total_abundance_order <- abund_d_order %>%
  group_by(Order) %>%
  summarize(Total_Abundance = sum(Abundance, na.rm = TRUE)) %>%
  arrange(desc(Total_Abundance))

# Select the top 10 microbial groups
top_15_microbes_order <- total_abundance_order %>%
  top_n(15, Total_Abundance) %>%
  pull(Order)

# Filter the data to keep only the top 20 microbial groups
abund_d_order<- abund_d_order %>%
  filter(Order %in% top_15_microbes_order)

# Reorder Microbial_Group factor levels based on total abundance
abund_d_order$Order <- factor(abund_d_order$Order, levels = top_15_microbes_order)


# Reorder Drought_Stage factor levels
abund_d_order$Drought_Stage <- factor(abund_d_order$Drought_Stage, levels = c("V2-Stage", "V4-Stage", "R1-Stage", "R4-Stage"))

# Plot the data
abundnace_plot_order = ggplot(abund_d_order, aes(x = Drought_Stage, y = Abundance, fill = Order)) +
  geom_bar(stat = "identity", position = "stack") +
  facet_wrap(~ Genotype) +
  theme_minimal() +
  labs(title = "Order",
       x = "Drought Stage",
       y = "Abundance",
       fill = "Order") +
  theme(axis.text.x = element_text(size = 14, angle = 45, hjust = 1),  # Adjust x-axis text size
        axis.text.y = element_text(size = 14),  # Adjust y-axis text size
        axis.title.x = element_text(size = 14), # Adjust x-axis label size
        axis.title.y = element_text(size = 14), # Adjust y-axis label size
        strip.text = element_text(size = 14),   # Adjust facet labels text size
        legend.title = element_text(size = 14), # Adjust legend title text size
        legend.text = element_text(size = 12),  # Adjust legend text size
        plot.title = element_text(size = 14, hjust = 0.5)) +  # Adjust title size and center it
  scale_fill_manual(values = top15_palette)  # Apply the custom color palette

abundnace_plot_order

ggsave("abundance_plot_order.jpeg", plot = abundnace_plot_order, 
       width = 12, height = 10, units = "in", 
       dpi = 1000)

###########################################################################################
############################ genotype X stage -- family #####################

abund_family = read.csv("gen-drought-family.csv")

# Reshape the data to long format
abund_d_family <- abund_family %>%
  pivot_longer(
    cols = -c(Drought_Stage, Genotype),    # Select columns after the first four factors
    names_to = "Family",
    values_to = "Abundance"
  )

total_abundance_family <- abund_d_family %>%
  group_by(Family) %>%
  summarize(Total_Abundance = sum(Abundance, na.rm = TRUE)) %>%
  arrange(desc(Total_Abundance))

# Select the top 10 microbial groups
top_15_microbes_family <- total_abundance_family %>%
  top_n(20, Total_Abundance) %>%
  pull(Family)

# Filter the data to keep only the top 20 microbial groups
abund_d_family<- abund_d_family %>%
  filter(Family %in% top_15_microbes_family)

# Reorder Microbial_Group factor levels based on total abundance
abund_d_family$Family <- factor(abund_d_family$Family, levels = top_15_microbes_family)


# Reorder Drought_Stage factor levels
abund_d_family$Drought_Stage <- factor(abund_d_family$Drought_Stage, levels = c("V2-Stage", "V4-Stage", "R1-Stage", "R4-Stage"))

# Plot the data
abundnace_plot_family = ggplot(abund_d_family, aes(x = Drought_Stage, y = Abundance, fill = Family)) +
  geom_bar(stat = "identity", position = "stack") +
  facet_wrap(~ Genotype) +
  theme_minimal() +
  labs(title = "Family",
       x = "Drought Stage",
       y = "Abundance",
       fill = "Family") +
  theme(axis.text.x = element_text(size = 14, angle = 45, hjust = 1),  # Adjust x-axis text size
        axis.text.y = element_text(size = 14),  # Adjust y-axis text size
        axis.title.x = element_text(size = 14), # Adjust x-axis label size
        axis.title.y = element_text(size = 14), # Adjust y-axis label size
        strip.text = element_text(size = 14),   # Adjust facet labels text size
        legend.title = element_text(size = 14), # Adjust legend title text size
        legend.text = element_text(size = 12),  # Adjust legend text size
        plot.title = element_text(size = 14, hjust = 0.5)) +  # Adjust title size and center it
  scale_fill_manual(values = top15_palette)  # Apply the custom color palette

abundnace_plot_family

ggsave("abundance_plot_family.jpeg", plot = abundnace_plot_family, 
       width = 12, height = 10, units = "in", 
       dpi = 1000)

###########################################################################################
###################genotype X stage -- genus


abund_genus = read.csv("gen-drought-genus.csv")

# Reshape the data to long format
abund_d_genus <- abund_genus %>%
  pivot_longer(
    cols = -c(Drought_Stage, Genotype),    # Select columns after the first four factors
    names_to = "Genus",
    values_to = "Abundance"
  )

total_abundance_genus <- abund_d_genus %>%
  group_by(Genus) %>%
  summarize(Total_Abundance = sum(Abundance, na.rm = TRUE)) %>%
  arrange(desc(Total_Abundance))

# Select the top 10 microbial groups
top_15_microbes_genus <- total_abundance_genus %>%
  top_n(20, Total_Abundance) %>%
  pull(Genus)

# Filter the data to keep only the top 20 microbial groups
abund_d_genus<- abund_d_genus %>%
  filter(Genus %in% top_15_microbes_genus)

# Reorder Microbial_Group factor levels based on total abundance
abund_d_genus$Genus <- factor(abund_d_genus$Genus, levels = top_15_microbes_genus)

# Reorder Genus factor levels with "Others" first, followed by the top 20 genera
abund_d_genus$Genus <- factor(abund_d_genus$Genus, levels = c("Others", top_15_microbes_genus))


# Reorder Drought_Stage factor levels
abund_d_genus$Drought_Stage <- factor(abund_d_genus$Drought_Stage, levels = c("V2-Stage", "V4-Stage", "R1-Stage", "R4-Stage"))

# Plot the data
abundnace_plot_genus = ggplot(abund_d_genus, aes(x = Drought_Stage, y = Abundance, fill = Genus)) +
  geom_bar(stat = "identity", position = "stack") +
  facet_wrap(~ Genotype) +
  theme_minimal() +
  labs(title = "Genus",
       x = "Drought Stage",
       y = "Abundance",
       fill = "Genus") +
  theme(axis.text.x = element_text(size = 14, angle = 45, hjust = 1),  # Adjust x-axis text size
        axis.text.y = element_text(size = 14),  # Adjust y-axis text size
        axis.title.x = element_text(size = 14), # Adjust x-axis label size
        axis.title.y = element_text(size = 14), # Adjust y-axis label size
        strip.text = element_text(size = 14),   # Adjust facet labels text size
        legend.title = element_text(size = 14), # Adjust legend title text size
        legend.text = element_text(size = 12),  # Adjust legend text size
        plot.title = element_text(size = 14, hjust = 0.5)) +  # Adjust title size and center it
  scale_fill_manual(values = top15_palette)  # Apply the custom color palette

abundnace_plot_genus

ggsave("abundance_plot_genus.jpeg", plot = abundnace_plot_genus, 
       width = 12, height = 10, units = "in", 
       dpi = 1000)


abundance_stages = egg::ggarrange(
  abundnace_plot_phylum, 
  abundnace_plot_class,
  abundnace_plot_order, 
  abundnace_plot_family, 
  abundnace_plot_genus, 
  ncol = 2,
  heights = c(1, 1, 1, 1, 1), # Adjust these values as needed
  widths = c(1, 1)            # Adjust these values as needed
)

ggsave("abundance_stages.jpeg", plot = abundance_stages, 
       width = 16, height = 10, units = "in", 
       dpi = 1000)

###############################################################################################
###############################################################################################
################################# treatment abundance plots

ab_trt_phylum = readxl::read_xlsx("trt-gen-phylum.xlsx")

head(ab_trt_phylum)

# Reshape the data to long format
ab_trt_d_phylum <- ab_trt_phylum %>%
  pivot_longer(
    cols = -c(Treatment, Genotype),    # Select columns after the first four factors
    names_to = "Phylum",
    values_to = "Abundance"
  )

total_abundance_phylum_trt <- ab_trt_d_phylum %>%
  group_by(Phylum) %>%
  summarize(Total_Abundance = sum(Abundance, na.rm = TRUE)) %>%
  arrange(desc(Total_Abundance))

# Select the top 10 microbial groups
top_15_microbes_phylum_trt <- total_abundance_phylum_trt %>%
  top_n(20, Total_Abundance) %>%
  pull(Phylum)

# Filter the data to keep only the top 20 microbial groups
ab_trt_d_phylum <- ab_trt_d_phylum %>%
  filter(Phylum %in% top_15_microbes_phylum_trt)

# Reorder Microbial_Group factor levels based on total abundance
ab_trt_d_phylum$Phylum <- factor(ab_trt_d_phylum$Phylum, levels = top_15_microbes_phylum_trt)

# Plot the data
abundance_plot_phylum_trt = ggplot(ab_trt_d_phylum, aes(x = Treatment, y = Abundance, fill = Phylum)) +
  geom_bar(stat = "identity", position = "stack") +  # Use position = "stack" for stacked bars
  facet_wrap(~ Genotype) +
  theme_minimal() +
  labs(title = "Phylum Abundance - Treatment x Genotype",
       x = "Treatment",
       y = "Abundance",
       fill = "Phylum") +
  theme(axis.text.x = element_text(size = 14, angle = 45, hjust = 1),  # Adjust x-axis text size
        axis.text.y = element_text(size = 14),  # Adjust y-axis text size
        axis.title.x = element_text(size = 14), # Adjust x-axis label size
        axis.title.y = element_text(size = 14), # Adjust y-axis label size
        strip.text = element_text(size = 14),   # Adjust facet labels text size
        legend.title = element_text(size = 14), # Adjust legend title text size
        legend.text = element_text(size = 12),  # Adjust legend text size
        plot.title = element_text(size = 14, hjust = 0.5)) +  # Adjust title size and center it
  scale_fill_manual(values = top15_palette)   # Apply the custom color palette

# Display the plot
print(abundance_plot_phylum_trt)


ggsave("abundance_plot_phylum_trt.jpeg", plot = abundance_plot_phylum_trt, 
       width = 8, height = 6, units = "in", 
       dpi = 1000)



################################################################################################
################################################################################################
################## Beta Diversity Plots ###############################################

library(ggplot2)

beta = readxl::read_xlsx("PCOA-Bray.xlsx")

beta_plot_genotype <- ggplot(beta, aes(x = Axis.1, y = Axis.2, color = Genotype)) +
  geom_point(size = 3) +   # Adjust point size as needed
  scale_color_manual(values = c("#599ec4", "#d65f3d")) +  # Customize colors
  scale_shape_manual(values = c(16, 17)) +
  scale_fill_manual(values = c("#599ec4", "#d65f3d")) +
  stat_ellipse(aes(fill = Genotype), level = 0.97, alpha = 0.2, geom = "polygon") +  # 95% confidence ellipses
  theme(legend.position = "right") + # Customize shapes
  geom_hline(yintercept = 0, linetype = "dashed", color = "black", size = 0.8) +  # Horizontal axis line
  geom_vline(xintercept = 0, linetype = "dashed", color = "black", size = 0.8) +  # Vertical axis line
  theme_minimal() +
  theme(
    axis.line = element_line(color = "black", size = 1),  # Ensure x and y axis lines are visible
    axis.ticks = element_line(color = "black", size = 1)  # Ensure axis ticks are visible
  ) +
  labs(x = "PC1 (37.2 %)",
       y = "PC2 (9.8 %)",
       title = "")

beta_plot_genotype

ggsave("beta_plot_genotype.jpeg", plot = beta_plot_genotype, 
       width = 6, height = 4, units = "in", 
       dpi = 1000)

beta_plot_trt <- ggplot(beta, aes(x = Axis.1, y = Axis.2, color = Treatment)) +
  geom_point(size = 3) +   # Adjust point size as needed
  scale_color_manual(values = c("#312A56", "#D82632")) +  # Customize colors
  scale_shape_manual(values = c(16, 17)) +
  scale_fill_manual(values = c("#312A56", "#D82632")) +
  stat_ellipse(aes(fill = Treatment), level = 0.97, alpha = 0.2, geom = "polygon") +  # 95% confidence ellipses
  theme(legend.position = "right") + # Customize shapes
  geom_hline(yintercept = 0, linetype = "dashed", color = "black", size = 0.8) +  # Horizontal axis line
  geom_vline(xintercept = 0, linetype = "dashed", color = "black", size = 0.8) +  # Vertical axis line
  theme_minimal() +
  theme(
    axis.line = element_line(color = "black", size = 1),  # Ensure x and y axis lines are visible
    axis.ticks = element_line(color = "black", size = 1)  # Ensure axis ticks are visible
  ) +
  labs(x = "PC1 (37.2 %)",
       y = "PC2 (9.8 %)",
       title = "")

beta_plot_trt

ggsave("beta_plot_trt.jpeg", plot = beta_plot_trt, 
       width = 6, height = 4, units = "in", 
       dpi = 1000)

beta_stage = readxl::read_xlsx("PCOA-Wg.xlsx")
# Ensure Drought_Stage is a factor with the correct order
beta_stage$Drought_Stage <- factor(beta_stage$Drought_Stage, levels = c("V2-Stage", "V4-Stage", "R1-Stage", "R4-Stage"))

beta_plot_stage <- ggplot(beta_stage, aes(x = Axis.1, y = Axis.2, color = Drought_Stage, fill = Drought_Stage)) +
  geom_point(size = 3, shape = 21) +  # Use shape = 21 for filled points
  scale_color_manual(values = c("lightblue", "skyblue", "royalblue", "midnightblue")) +  # Gradual color change
  scale_fill_manual(values = c("lightblue", "skyblue", "royalblue", "midnightblue")) +  # Matching fill for ellipses
  stat_ellipse(level = 0.95, alpha = 0.2, geom = "polygon") +  # 97% confidence ellipses with color fill
  theme(legend.position = "right") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black", size = 0.8) +  # Horizontal axis line
  geom_vline(xintercept = 0, linetype = "dashed", color = "black", size = 0.8) +  # Vertical axis line
  theme_minimal() +
  theme(
    axis.line = element_line(color = "black", size = 1),  # Ensure x and y axis lines are visible
    axis.ticks = element_line(color = "black", size = 1)  # Ensure axis ticks are visible
  ) +
  labs(x = "PC1 (37.2 %)",
       y = "PC2 (9. 8%)",
       title = "")

ggsave("beta_plot_stage.jpeg", plot = beta_plot_stage, 
       width = 6, height = 4, units = "in", 
       dpi = 1000)

############################################################################################
############################# alpha diversity plots ################################

chao1_16s = read.csv("chao1.csv")
shannon_16s = read.csv("shannon.csv")

# Set the order of the Drought_Stage factor
shannon_16s$Drought_Stage <- factor(shannon_16s$Drought_Stage, levels = c("V2-Stage", "V4-Stage", "R1-Stage", "R4-Stage"))

chao1_16s$Drought_Stage <- factor(chao1_16s$Drought_Stage, levels = c("V2-Stage", "V4-Stage", "R1-Stage", "R4-Stage"))

stages_chao1_16s = ggplot(chao1_16s, aes(x = Drought_Stage, y = value, fill = Treatment)) +
  geom_boxplot(position = position_dodge(0.8)) + # Adjust box plot positions
  geom_jitter(position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.8), alpha = 0.5) + # Add jittered points
  facet_wrap(~ Genotype) +
  labs(title = "",
       x = "Drought Stage",
       y = "Chao1 Index") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 14), # Increase x-axis text size
        axis.text.y = element_text(size = 14), # Increase y-axis text size
        axis.title.x = element_text(size = 16), # Increase x-axis label size
        axis.title.y = element_text(size = 16), # Increase y-axis label size
        strip.text = element_text(size = 14), # Increase facet label size
        panel.border = element_rect(colour = "black", fill = NA, size = 1)) +# Add border
  scale_fill_manual(values = alpha_pallete)

stages_chao1_16s

ggsave("stages chao1.pdf", plot = stages_chao1_16s, device = "pdf", width = 8, height = 6, units = "in", dpi = 1000)

stages_shannon_16s = ggplot(shannon_16s, aes(x = Drought_Stage, y = value, fill = Treatment)) +
  geom_boxplot(position = position_dodge(0.8)) + # Adjust box plot positions
  geom_jitter(position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.8), alpha = 0.5) + # Add jittered points
  facet_wrap(~ Genotype) +
  labs(title = "",
       x = "Drought Stage",
       y = "Shannon Diversity Index") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 14), # Increase x-axis text size
        axis.text.y = element_text(size = 14), # Increase y-axis text size
        axis.title.x = element_text(size = 16), # Increase x-axis label size
        axis.title.y = element_text(size = 16), # Increase y-axis label size
        strip.text = element_text(size = 14), # Increase facet label size
        panel.border = element_rect(colour = "black", fill = NA, size = 1)) +# Add border
  scale_fill_manual(values = alpha_pallete)

stages_shannon_16s

ggsave("stages shannon 16s.pdf", plot = stages_shannon_16s, device = "pdf", width = 8, height = 6, units = "in", dpi = 1000)


alpha_metrics = readxl::read_xlsx("Alpha_Diversity_16s.xlsx")

shapiro.test(alpha_metrics$Chao1)

# Shapiro-Wilk test for Shannon
shapiro.test(alpha_metrics$Shannon)

# GLMM for Chao1
glmm_chao1 <- lmer(Chao1 ~ Treatment * Drought_Stage * Genotype + (1|Replication), data = alpha_metrics)
glmm_shannon <- lmer(Shannon ~ Treatment * Drought_Stage * Genotype + (1|Replication), data = alpha_metrics)



summary(glmm_chao1)

anova_results_chao1_16s <- anova(glmm_chao1)
anova_results_shannon_16s <- anova(glmm_shannon)

anova_results_chao1_16s 
anova_results_shannon_16s
# GLMM for Shannon
glmm_shannon <- lmer(Shannon ~ Group + (1|RandomFactor), data = data)
summary(glmm_shannon)

##### Soil Data box plot

setwd("~/OneDrive - Mississippi State University/Soybean project_Mahesh/Sujan_Drought/Bioinformatics/8-24-24/16S")

soil = readxl::read_xlsx("soil_data.xlsx")

# Set the order of the Drought_Stage factor
soil$Drought_Stage <- factor(soil$Drought_Stage, levels = c("V2-Stage", "V4-Stage", "R1-Stage", "R4-Stage"))

soil_pH = ggplot(soil, aes(x = Drought_Stage, y = pH, fill = Treatment)) +
  geom_boxplot(position = position_dodge(0.8)) + # Adjust box plot positions
  geom_jitter(position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.8), alpha = 0.5) + # Add jittered points
  facet_wrap(~ Genotype) +
  labs(title = "",
       x = "Drought Stage",
       y = "soil pH") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 14), # Increase x-axis text size
        axis.text.y = element_text(size = 14), # Increase y-axis text size
        axis.title.x = element_text(size = 16), # Increase x-axis label size
        axis.title.y = element_text(size = 16), # Increase y-axis label size
        strip.text = element_text(size = 14), # Increase facet label size
        panel.border = element_rect(colour = "black", fill = NA, size = 1)) +# Add border
  scale_fill_manual(values = alpha_pallete)

soil_pH

NAG = ggplot(soil, aes(x = Drought_Stage, y = NAG, fill = Treatment)) +
  geom_boxplot(position = position_dodge(0.8)) + # Adjust box plot positions
  geom_jitter(position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.8), alpha = 0.5) + # Add jittered points
  facet_wrap(~ Genotype) +
  labs(title = "",
       x = "Drought Stage",
       y = "β-glucosaminidase activity") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 14), # Increase x-axis text size
        axis.text.y = element_text(size = 14), # Increase y-axis text size
        axis.title.x = element_text(size = 16), # Increase x-axis label size
        axis.title.y = element_text(size = 16), # Increase y-axis label size
        strip.text = element_text(size = 14), # Increase facet label size
        panel.border = element_rect(colour = "black", fill = NA, size = 1)) +# Add border
  scale_fill_manual(values = alpha_pallete)

NAG
