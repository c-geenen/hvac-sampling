library(tidyverse)
library(ggtext)

# Function for determining exact binomial confidence interval
binom_CI <- function(x, n, low_high) {
  if (n>0) {
    CI <- binom.test(x, n, 0/100, alternative =  "two.sided", conf.level = 0.95)$conf.int[low_high]
    return(CI)
  } else {
    return(c(0,1)[l])
  }
}

chi_squared <- function(x1, x2, y1, y2) {
  p <- chisq.test( matrix(c(x1, x2,y1, y2), ncol = 2), correct=TRUE)$p.value
  return(p)
}

data <- read_csv("input/test_data.csv") %>%
  select(-Ct_value) %>%
  
  # only use SARS-CoV-2 N gene
  filter(
    pathogen != "SARS-CoV-2 (respiratory panel)",
    pathogen != "SARS-CoV-2 (TaqPath ORF1ab)",
    pathogen != "SARS-CoV-2 (TaqPath S-protein)"
  ) %>%
  
  pivot_wider(
    names_from = "location",
    values_from = detected
  ) %>%
  mutate(
    detected_anywhere = replace_na(`1` | `2` | `3` | HVAC, F),
    detected_room = replace_na(`1` | `2` | `3`, F)
  )

############
# Sensitivity / specificity of air samples for detecting positive nose wupes
###############

sensitivity <- data %>%
  # Add "all pathogens" category
  mutate(pathogen="-All pathogens-") %>%
  bind_rows(data) %>%
  
  # Determine sensitivity for each pathogen
  filter(detected_anywhere) %>%
  group_by(pathogen) %>%
  summarise(
    sens_hvac_pos = sum(HVAC, na.rm=T),
    sens_hvac_n = sum(!is.na(HVAC)),
    sens_hvac = sens_hvac_pos / sens_hvac_n,
    sens_hvac_low = binom_CI(sens_hvac_pos, sens_hvac_n, 1),
    sens_hvac_high = binom_CI(sens_hvac_pos, sens_hvac_n, 2),
    sens_room1_pos = sum(`1`, na.rm=T),
    sens_room1_n = sum(!is.na(`1`)),
    sens_room1 = sens_room1_pos / sens_room1_n,
    sens_room1_low = binom_CI(sens_room1_pos, sens_room1_n, 1),
    sens_room1_high = binom_CI(sens_room1_pos, sens_room1_n, 2),
    sens_room1_p = suppressWarnings(chi_squared(sens_room1_pos,sens_room1_n-sens_room1_pos,sens_hvac_pos,sens_hvac_n-sens_hvac_pos)),
    sens_room2_pos = sum(`2`, na.rm=T),
    sens_room2_n = sum(!is.na(`2`)),
    sens_room2 = sens_room2_pos / sens_room2_n,
    sens_room2_low = binom_CI(sens_room2_pos, sens_room2_n, 1),
    sens_room2_high = binom_CI(sens_room2_pos, sens_room2_n, 2),
    sens_room2_p = suppressWarnings(chi_squared(sens_room2_pos,sens_room2_n-sens_room2_pos,sens_hvac_pos,sens_hvac_n-sens_hvac_pos)),
    sens_room3_pos = sum(`3`, na.rm=T),
    sens_room3_n = sum(!is.na(`3`)),
    sens_room3 = sens_room3_pos / sens_room3_n,
    sens_room3_low = binom_CI(sens_room3_pos, sens_room3_n, 1),
    sens_room3_high = binom_CI(sens_room3_pos, sens_room3_n, 2),
    sens_room3_p = suppressWarnings(chi_squared(sens_room3_pos,sens_room3_n-sens_room3_pos,sens_hvac_pos,sens_hvac_n-sens_hvac_pos)),
  ) %>%
  select(-starts_with("n_")) %>%
  
  pivot_longer(
    cols = sens_hvac_pos:sens_room3_p,
    names_to = c("location","ci"),
    names_prefix = "sens_",
    names_sep = "_"
  ) %>%
  mutate(ci=replace_na(ci,"est")) %>%
  pivot_wider(
    names_from = "ci",
    values_from = "value"
  ) %>%
  
  mutate(
    colour=case_match(location, c("room3","room1","room2")~"Room","hvac"~"HVAC"),
    location =case_match(
      location,
      "hvac" ~ "HVAC",
      "room1" ~ "1",
      "room2" ~ "2",
      "room3" ~ "3"
    ),
    location = relevel(as.factor(location),"HVAC"),
    significant=if_else(p<0.05,"*","")
  )

write_csv(sensitivity,"output/sensitivity.csv")
  
###################
# Plot sensitivity results
####################

ggplot(sensitivity, aes(location, est, ymin=low, ymax=high, colour=colour)) +
  geom_point() +
  geom_errorbar() +
  facet_wrap(~pathogen, ncol=4) +
  labs(
    colour="Sampling location",
    y="Sensitivity",
    x="Sampling location"
  ) +
  scale_colour_manual(values=c("red","darkgrey")) +
  scale_y_continuous(limits=c(0,1)) +
  geom_text(aes(x=location, y=high, label=significant), nudge_y=0.05, size=8, colour="black") +
  geom_text(data=filter(sensitivity, location=="1"),aes(label=paste0("n=",n)), x=4, y=0.07, colour="black") +
  theme_bw() + 
  theme(
    legend.position = c(0.87,0.05),
    strip.text = ggtext::element_markdown()
  )

ggsave("output/supp_figure_sensitivity.png", width=10, height=10)
ggsave("output/supp_figure_sensitivity.pdf", width=10, height=10)

# only frequently detected pathogens
ggplot(filter(sensitivity,n>=10), aes(location, est, ymin=low, ymax=high, colour=colour)) +
  geom_point() +
  geom_errorbar() +
  facet_wrap(~pathogen, ncol=4) +
  labs(
    colour="Sampling location",
    y="Sensitivity",
    x="Sampling location"
  ) +
  scale_colour_manual(values=c("red","darkgrey")) +
  scale_y_continuous(limits=c(0,1)) +
  scale_x_discrete(breaks=c("HVAC","1","2","3")) +
  geom_text(aes(x=location, y=high, label=significant), nudge_y=0.05, size=8, colour="black") +
  geom_text(data=filter(sensitivity, n>=10, location=="1"),aes(label=paste0("n=",n)), x=4, y=0.07, colour="black") +
  theme_bw() + 
  theme(
    legend.position = c(0.87,0.1),
    strip.text = ggtext::element_markdown()
  )

ggsave("output/figure_3.png", width=9, height=5)
ggsave("output/figure_3.pdf", width=9, height=5)

