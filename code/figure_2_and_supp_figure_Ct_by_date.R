library(tidyverse)
library(ggtext)

# get environmental variables
envir_data <- read_csv("input/environmental_data.csv") %>%
  select(-sensor_samples) %>%
  pivot_wider(names_from="variable",values_from="mean") %>%
  group_by(samplingDate) %>%
  mutate(
    occupancy=(mean(occupancy, na.rm=T)),
    airSamplingStart=(mean(airSamplingStart, na.rm=T)),
    airSamplingEnd=(mean(airSamplingEnd, na.rm=T)),
    samplingDuration=(mean(samplingDuration, na.rm=T))
  ) %>%
  ungroup %>%
  filter(location!="roomAll")
  
#################
# Ct in air vs nose wipes
####################

PCR_params <- read_csv("input/PCR_parameters.csv")

test_data <- read_csv("input/test_data.csv") %>%
  
  # Only TaqPath N
  filter(
    pathogen != "SARS-CoV-2 (TaqPath ORF1ab)",
    pathogen != "SARS-CoV-2 (TaqPath S-protein)",
    pathogen != "SARS-CoV-2 (respiratory panel)"
  ) %>%
  
  # Remove pathogens which were never positive
  group_by(pathogen) %>%
  filter(sum(detected)>0) %>%
  ungroup() %>%
  
  # Add Ct at LOD
  left_join(PCR_params,by="pathogen") %>%
  
  arrange(samplingDate) %>%
  mutate(
    sample_type=if_else(
      location=="HVAC",
      "HVAC",
      "Room"
    ),
    samplingDate_str=paste0(day(samplingDate),"/",month(samplingDate))
  )

###########################
# Plot ct values
########################
    
test_data_recoded <- test_data %>%    
  mutate(
    samplingDate_str=fct_relevel(samplingDate_str,unique(test_data$samplingDate_str))
  ) %>%
  
  arrange(desc(sample_type))

add_not_detected_label <- function(label) {
  l = as.character(label)
  return(
    if_else(
      as.numeric(l)==47.5,
      "Not detected",
      if_else(
        as.numeric(l)==45,
        "",
        l
      )
    )
  )
}

test_data_not_detected <- test_data_recoded %>%
  filter(is.na(Ct_value)) %>%
  mutate(Ct_value=47.5)

ct_plot <- ggplot(filter(test_data_recoded,!is.na(Ct_value)), aes(samplingDate_str,Ct_value,colour=sample_type)) +
  geom_jitter(height=0,width=0.3,size=1.5, alpha=0.5) +
  facet_wrap(~pathogen, ncol=2) +
  theme_bw() +
  scale_y_reverse(
    labels = add_not_detected_label,
    breaks=c(seq(0, 45, by = 5),47.5)
  ) +
  scale_x_discrete(expand=expansion(add=0.5)) +
  scale_colour_manual(values=c("red","darkgrey")) +
  theme(
    panel.grid.major.x = element_blank(),
    strip.text = ggtext::element_markdown(),
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
    legend.position = c(0.8,0.02)
  ) +
  geom_jitter(data=test_data_not_detected, aes(samplingDate_str,Ct_value), height=1.5, width=0.3, size=1.5, alpha=0.5) +
  geom_hline(yintercept = 45) +
  geom_hline(aes(yintercept = Ct_at_LOD),colour="black",alpha=0.5,linetype="dashed") +
  geom_vline(xintercept=seq(0.5,30.5,by=1)) +
  labs(colour="Sampling location",y="Air sample Ct value",x="Sampling date")
ct_plot

ggsave("output/supp_figure_ct_by_date.png", width=10, height=10)

##################
# Heatmap
##############

heatmap_data <- test_data_recoded %>%
  group_by(pathogen,samplingDate_str) %>%
  summarise(
    detected_hvac = sum(sample_type=="HVAC" & detected)>0,
    detected_room = sum(grepl("Room",sample_type,ignore.case=T) & detected),
    n_room = sum(grepl("Room",sample_type,ignore.case=T)),
  ) %>%
  mutate(total_detected = sum(detected_room>0) + sum(detected_hvac) ) %>%
  ungroup %>%
  
  # arrange pathogens according to number of times detected
  arrange(total_detected) %>%
  mutate(pathogen = factor(pathogen, levels=unique(pathogen))) %>%
  
  mutate(
    detected_room = fct_recode(
      as.factor(detected_room),
      "Not detected" = "0",
      "1 sampler" = "1",
      "2 samplers" = "2",
      "3 samplers" = "3",
    )
  )
  
ggplot(select(heatmap_data, -detected_hvac), aes(samplingDate_str, pathogen, fill=detected_room)) +
  geom_tile() +
  scale_fill_manual(values=c("#ffffff","#ffd7cf","#ff9782","#ff2f05")) +
  labs(
    y="Pathogen",
    x="Sampling date",
    fill = "Detected in room"
  ) +
  ggnewscale::new_scale_fill() + 
  geom_point(data=filter(heatmap_data,detected_hvac), aes(shape=detected_hvac), size=2) +
  labs(
    shape = "Detected in HVAC"
  ) +
  coord_fixed() +
  theme_bw() +
  theme(
    axis.text.y=element_markdown(),
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)
  )

ggsave("output/figure_2.png", width=8, height=4)

