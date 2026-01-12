# Bioassay data analysis of the Gambia phenotypic data

# load libraries
library(pacman)
p_load(readxl, tidyr, tidyverse, ggplot2, forcats, grid, dplyr, viridis, extrafont)


# load data
pheno_data = read_excel("./Mortality rate_Saad_20211012_BSA.xlsx", sheet = 1)
pheno_data

# clean out data
pheno_clean <- pheno_data %>%
  select(!c(Upper, Lower)) %>%
  drop_na()


# map villages to administrative region 1
pheno_clean <- pheno_clean %>%
  mutate(admin1 = case_when(
    VILLAGE %in% c("TALLINDING", "EBOTOWN", "BAKAU") ~ "Kanifing",
    VILLAGE %in% c("BRIKAMA", "BULLOCK") ~ "Brikama",
    VILLAGE %in% c("FARAFENNI", "ESSAU", "YALLAL BA") ~ "Kerewan",
    VILLAGE %in% c("MADINA", "DONGOROBA") ~ "Mansakonko",
    VILLAGE %in% c("KUMBIJA", "KATAMIMA", "SARE SILLI", "JAHALI", "KEREWAN", "BANSANG", "SAPU", "JARENG", "SARUJA", "BOIRAM", "MADINA UMFALLY", "SANKULAY KUNDA") ~ "Janjanbureh",
    VILLAGE %in% c("NYANGA BANTANG", "KUNTAUR") ~ "Kuntaur",
    TRUE ~ "Basse"
  ))


# map villages to administrative region 2
pheno_clean <- pheno_clean %>%
  mutate(admin2 = case_when(
    VILLAGE %in% c("TALLINDING", "EBOTOWN") ~ "Jeshwang",
    VILLAGE %in% c("BAKAU") ~ "Bakau",
    VILLAGE %in% c("BULLOCK") ~ "Foni Brefet",
    VILLAGE %in% c("BRIKAMA") ~ "Kombo Central",
    VILLAGE %in% c("FARAFENNI") ~ "Illiasa",
    VILLAGE %in% c("ESSAU") ~ "Lower Niumi",
    VILLAGE %in% c("YALLAL BA") ~ "Upper Baddibu",
    VILLAGE %in% c("DONGOROBA") ~ "Jarra East",
    VILLAGE %in% c("MADINA") ~ "Jarra West",
    VILLAGE %in% c("KUMBIJA", "KATAMIMA", "SARE SILLI", "JAHALI", "KEREWAN") ~ "Niamina East",
    VILLAGE %in% c("BANSANG") ~ "Upper Fulladu West",
    VILLAGE %in% c("SAPU") ~ "Lower Fulladu West",
    VILLAGE %in% c("JARENG", "SARUJA", "BOIRAM") ~ "Niamina West",
    VILLAGE %in% c("MADINA UMFALLY", "SANKULAY KUNDA") ~ "Niamina Dankunku",
    VILLAGE %in% c("NYANGA BANTANG", "KUNTAUR") ~ "Niani",
    TRUE ~ "Basse"
  ))


# create region2 column (inland and coastal)
pheno_clean <- pheno_clean %>%
  mutate(CC_IC = case_when(
    admin1 %in% c("Kanifing", "Kerewan", "Brikama") ~ "CC",
    TRUE ~ "IC"
  ))

# remove some insecticides and add % to the insecticide column
pheno_clean <- pheno_clean %>%
  filter(!INSECTICIDE %in% c("Bendiocarb 0.5", "Bendiocarb 1", "Deltamethrin 0.5")) %>%
  mutate(INSECTICIDE = case_when(
    INSECTICIDE == "Bendiocarb 0.1" ~ "Bendiocarb 0.1%",
    INSECTICIDE == "DDT 4" ~ "DDT 4.0%",
    INSECTICIDE == "Deltamethrin 0.05" ~ "Deltamethrin 0.05%",
    INSECTICIDE == "Deltamethrin 0.25" ~ "Deltamethrin 0.25%",
    INSECTICIDE == "Permethrin 0.75" ~ "Permethrin 0.75%",
    INSECTICIDE == "Permethrin 3.75" ~ "Permethrin 3.75%",
    INSECTICIDE == "Permethrin 7.5" ~ "Permethrin 7.5%",
    INSECTICIDE == "pirimiphos Methyl 0.25" ~ "Pirimiphosmethyl 0.25%"
  ))
  

# classify resistance levels based on the mortality percentage
df <- pheno_clean %>%
  mutate(resistance_level = case_when(
    `MORTALITY %` >= 98 ~ "Susceptible",
    `MORTALITY %` >= 90 & `MORTALITY %` < 98 ~ "Possible resistance",
    `MORTALITY %` < 90 & `N TOTAL` >= 100 ~ "Confirmed resistance", 
    `MORTALITY %` < 90 & `N TOTAL` < 100 ~ "Possible resistance",
    TRUE ~ NA_character_
  ))


# aggregate data by admin1 and INSECTICIDE - (alternative way?)
df_aggregated <- pheno_clean %>%
  group_by(CC_IC, admin1, INSECTICIDE) %>%
  summarise(
    total_exposed = sum(`N TOTAL`),
    total_dead = sum(DEAD),
    mortality_rate = (total_dead / total_exposed) * 100,  # calculate mortality percentage
    # retain original 'N TOTAL', 'DEAD' and 'MORTALITY %'values as concatenated strings - maybe not a good idea, need to find a better way
    original_n_total = paste(`N TOTAL`, collapse = ", "),  
    original_dead = paste(DEAD, collapse = ", "),
    original_mortality_pct = paste(`MORTALITY %`, collapse = ", ")
  ) %>%
  ungroup()

# classify resistance levels based on the corrected mortality percentage
df_aggregated <- df_aggregated %>%
  mutate(resistance_level = case_when(
    mortality_rate >= 98 ~ "Susceptible",
    mortality_rate >= 90 & mortality_rate < 98 ~ "Possible resistance",
    mortality_rate < 90  ~ "Confirmed resistance",
    TRUE ~ NA_character_
  ))

# # calculate resistance levels based on the aggregated mortality rate
# df_aggregated <- df_aggregated %>%
#   mutate(resistance_level = case_when(
#     mortality_rate == 100 ~ "Susceptible",
#     mortality_rate >= 98 & mortality_rate < 100 ~ "Confirmed resistance",
#     mortality_rate >= 90 & mortality_rate < 98 ~ "Moderate-high intensity resistance",
#     mortality_rate < 90 ~ "Low intensity resistance",
#     TRUE ~ NA_character_
#   ))

# calculate 95% confidence intervals for mortality using binomial exact test
df_aggregated <- df_aggregated %>%
  rowwise() %>%
  mutate(
    binom_test = list(binom.test(total_dead, total_exposed, conf.level = 0.95)),
    ci_lower = binom_test$conf.int[1] * 100,
    ci_upper = binom_test$conf.int[2] * 100
  ) %>%
  ungroup() %>%
  select(-binom_test)

# convert 'admin1', 'region2', and 'resistance_level' to factors for plotting
df_aggregated <- df_aggregated %>%
  mutate(admin1 = factor(admin1, levels = c("Kanifing", "Brikama", "Kerewan", "Mansakonko", "Janjanbureh", "Kuntaur", "Basse")),
         # admin2 = factor(admin2, levels = c("Bakau", "Jeshwang", "Foni Brefet", "Kombo Central", "Illiasa", "Lower Niumi", 
         #                                    "Upper Baddibu", "Jarra East", "Jarra West", "Lower Fulladu West", "Niamina Dankunku", 
         #                                    "Niamina East", "Niamina West", "Upper Fulladu West", "Niani", "Basse")), 
         #resistance_level = factor(resistance_level, levels = c("Susceptible", "Low intensity resistance", "Moderate-high intensity resistance", "Confirmed resistance")))
         resistance_level = factor(resistance_level, levels = c( "Confirmed resistance", "Possible resistance", "Susceptible")),
         INSECTICIDE = factor(INSECTICIDE, levels = c("Deltamethrin 0.05%", "Permethrin 0.75%", "DDT 4.0%","Bendiocarb 0.1%", "Pirimiphosmethyl 0.25%",
                                                      "Deltamethrin 0.25%", "Permethrin 3.75%", "Permethrin 7.5%")))


#======================================Figure 1: Bar plot=======================================

fig1_bar <- ggplot(df_aggregated, aes(x = admin1, y = mortality_rate, fill = resistance_level)) +
  geom_bar(stat = "identity", position = "stack", width = 0.8) +
  geom_errorbar(aes(ymin = ci_lower, ymax = ci_upper), 
                width = 0.2, position = position_dodge(width = 0.8), color = "black") +
  labs(
    x = "Administrative Region", 
    y = "Mortality Rate (%) with 95% CI", 
    fill = "Resistance Level",
    #title = "Insecticide Resistance by Region"
  ) +
  scale_fill_manual(
    values = c("Confirmed resistance" = "#E69F00", 
               "Possible resistance" = "#56B4E9",
               "Susceptible" = "#009E73")
  ) +
  theme_classic() +
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 10),
    axis.text.y = element_text(size = 10),
    text = element_text(size = 12, family = "sans"),
    strip.text = element_text(size = 10, face = "bold"),
    panel.border = element_rect(colour = "black", fill = NA),
    legend.position = "bottom",
    legend.text = element_text(size = 10),
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold")
  ) +
  facet_wrap(~INSECTICIDE, nrow=2)


# Save plot to a PNG file
ggsave("./figures/Fig1_barplot_insecticide_resistance.png", plot = fig1_bar,
       width = 8, height = 6, dpi = 300, units = "in")


#======================================Figure 1: dot + error bar =======================================

fig1_dot <- ggplot(df_aggregated, aes(x = admin1, y = mortality_rate, color = resistance_level)) +
  geom_point(size = 3) +
  geom_errorbar(aes(ymin = ci_lower, ymax = ci_upper), width = 0.2) +
  labs(
    x = "Administrative Region", 
    y = "Mortality Rate (%) with 95% CI", 
    color = "Resistance Level",
    #title = "Insecticide Resistance by Region"
  ) +
  scale_color_manual(values = c(
    "Confirmed resistance" = "#E69F00", 
    "Possible resistance" = "#56B4E9",
    "Susceptible" = "#009E73"
  )) + 
  theme_classic() +
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 10),
    axis.text.y = element_text(size = 10),
    axis.title.x = element_blank(),
    text = element_text(size = 12, family = "sans"),
    strip.text = element_text(size = 10, face = "bold"),
    panel.border = element_rect(colour = "black", fill = NA),
    legend.position = "bottom",
    legend.text = element_text(size = 10),
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold")
  ) +
  facet_wrap(~INSECTICIDE, nrow=2)


# Save plot to a PNG file
ggsave("./figures/Fig1_dotplot_insecticide_resistance.png", plot = fig1_dot,
       width = 8, height = 6, dpi = 300, units = "in")

#=================================Supplementary Figure 1: heatmap (option 1)=====================================

# Optional: order regions by mean mortality
df_aggregated2 <- df_aggregated %>%
  group_by(admin1) %>%
  mutate(mean_mortality = mean(mortality_rate, na.rm = TRUE)) %>%
  ungroup() %>%
  mutate(admin1 = reorder(admin1, mean_mortality))

# Heatmap plot
supp_plot1_a <- ggplot(df_aggregated, aes(x = INSECTICIDE, y = admin1, fill = mortality_rate)) +
  geom_tile(color = "white", width = 0.9, height = 0.9) +
  geom_text(aes(label = round(mortality_rate, 1)), size = 3, color = "white") +
  scale_fill_viridis(
    name = "Mortality Rate (%)", 
    option = "C", 
    direction = -1, 
    limits = c(0, 100)
  ) +
  labs(x = "Insecticide", y = "Region") +
  theme_minimal(base_family = "sans") +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
    axis.text.y = element_text(size = 10),
    axis.title = element_text(size = 12, face = "bold"),
    legend.position = "right",
    legend.title = element_text(size = 11),
    legend.text = element_text(size = 10),
    panel.grid = element_blank(),
    plot.title = element_text(size = 14, face = "bold", hjust = 0.5)
  )

# Save plot to a PNG file
ggsave("./figures/Supp_Fig1_heatmap_insecticide_resistance.png", plot = supp_plot1_a,
       width = 8, height = 6, dpi = 300, units = "in")

#=================================Supplementary Figure 1: heatmap (option 2)=====================================

# Combine region and classification for row labels
df_aggregated3 <- df_aggregated %>%
  mutate(region_type = paste(admin1, CC_IC, sep = " - ")) %>%
  mutate(region_type = factor(region_type, levels = unique(region_type)))

# Plot
supp_plot1_b <- ggplot(df_aggregated3, aes(x = INSECTICIDE, y = region_type, fill = mortality_rate)) +
  geom_tile(color = "white", width = 0.9, height = 0.9) +
  geom_text(aes(label = round(mortality_rate, 1)), size = 3, color = "white") +
  scale_fill_viridis(
    name = "Mortality Rate (%)", 
    option = "C", 
    direction = -1, 
    limits = c(0, 100)
  ) +
  labs(
    x = "Insecticide", 
    y = "Region", 
    #title = "Insecticide Mortality Rates by Region and Classification"
  ) +
  theme_minimal(base_family = "sans") +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
    axis.text.y = element_text(size = 9),
    axis.title = element_text(size = 12, face = "bold"),
    legend.position = "right",
    legend.title = element_text(size = 11),
    legend.text = element_text(size = 10),
    panel.grid = element_blank(),
    plot.title = element_text(size = 14, face = "bold", hjust = 0.5)
  )

# Save plot to a PNG file
ggsave("./figures/Supp_Fig1_insecticide_resistance.png", plot = supp_plot1_b,
       width = 8, height = 6, dpi = 300, units = "in")

