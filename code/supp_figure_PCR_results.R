library(tidyverse)
library(ggtext)

############
# all tests
#########

test_data <- read_csv("input/test_data.csv") %>%
  
  # Only TaqPath N-protein
  filter(
    pathogen != "SARS-CoV-2 (TaqPath ORF1ab)",
    pathogen != "SARS-CoV-2 (TaqPath S-protein)",
    pathogen != "SARS-CoV-2 (respiratory panel)"
  ) %>%
  
  group_by(pathogen) %>%
  summarise(
    n=n(),
    n_pos=sum(detected),
    n_pos/n
  ) %>%
  mutate(pathogen=fct_reorder(pathogen,n_pos))

ggplot(test_data, aes(n_pos, pathogen)) +
  geom_col() + 
  theme_bw() +
  theme(
    axis.text.y = element_markdown()
  ) +
  labs(
    y="Pathogen",
    x="Positive results"
  )

#################
# Per location
##############@@

test_data <- read_csv("input/test_data.csv") %>%
  
  # Only TaqPath N-protein
  filter(
    pathogen != "SARS-CoV-2 (TaqPath ORF1ab)",
    pathogen != "SARS-CoV-2 (TaqPath S-protein)",
    pathogen != "SARS-CoV-2 (respiratory panel)"
  ) %>%
  
  group_by(pathogen) %>%
  mutate(n_detected=sum(detected)) %>%
  ungroup %>%
  arrange(n_detected) %>%
  mutate(pathogen=fct_reorder(pathogen,n_detected))

n = nrow(distinct(test_data,samplingDate,location))

ggplot(test_data, aes(y=pathogen, alpha=detected, fill=location)) +
  scale_alpha_manual(values=c(0,1)) +
  scale_x_continuous(
   # breaks=c(0,n*0.25,n*0.5,n*0.75,n),
    #labels=c("0","0.25","0.5","0.75","1")
  ) +
  geom_bar() +
  theme_bw() +
  theme(axis.text.y = element_markdown()) +
  guides(alpha = "none") +
  labs(y="Pathogen",x="Number of positive samples",fill="Sampler location")

ggsave("output/supp_figure_PCR_results.png", width=6, height=5)
