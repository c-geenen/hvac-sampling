library(tidyverse)

#################
# plot environmental data
###############################

envir_data <- read_csv("input/environmental_data.csv") %>%
  mutate(
    sample_type=if_else(
      location=="HVAC",
      "HVAC",
      "Room"
    ),
    variable=case_when(
      variable=="co2" ~ "CO2 concentration [ppm]",
      variable=="samplingDuration" ~ "Sampling duration [hours]",
      variable=="humidity" ~ "Relative humidity [%]",
      variable=="light" ~ "Light intensity [lux]",
      variable=="occupancy" ~ "Occupancy",
      variable=="temperature" ~ "Temperature [°C]",
      variable=="velocity" ~ "Air velocity [m/s]",
      variable=="airSamplingStart" ~ "Sampling start time [24h clock]",
      variable=="airSamplingEnd" ~ "Sampling end time [24h clock]",
      variable=="velocity" ~ "Air velocity [m/s]",
      .default = variable
    ),
    location=if_else(location=="roomAll","Room",location)
  ) %>%
  
  mutate(location = relevel(as.factor(location),"HVAC"))

# averages in room
envir_median_room <- envir_data %>%
  filter(location!="HVAC") %>%
  group_by(variable) %>%
  summarise(
    sd = sd(mean, na.rm=T),
    mean = mean(mean, na.rm=T)
  ) %>%
  print(n=50)

# averages in HVAC
envir_median_room <- envir_data %>%
  filter(location=="HVAC") %>%
  group_by(variable) %>%
  summarise(
    sd = sd(mean, na.rm=T),
    mean = mean(mean, na.rm=T)
  ) %>%
  print(n=50)

# mean by sensor
envir_by_sensor_ci <- envir_data %>%
  filter(variable!="Sampling duration [hours]") %>%
  group_by(location, variable) %>%
  summarise(
    ci_low = t.test(na.omit(mean))$conf.int[1],
    ci_high = t.test(na.omit(mean))$conf.int[2],
    iqr = IQR(mean, na.rm=T),
    median = median(mean, na.rm=T),
    mean = mean(mean, na.rm=T)
  ) %>%
  print(n=50)

ggplot(envir_data, aes(location, mean, colour=sample_type)) +
  geom_jitter(height=0, width=0.2, alpha=0.4, size=1) +
  scale_colour_manual(values=c("red","darkgrey")) +
  facet_wrap(~variable, scales="free") +
  geom_point(data=envir_by_sensor_ci, aes(location, mean), colour="black") +
  geom_errorbar(data=envir_by_sensor_ci, aes(location, mean, ymin=ci_low, ymax=ci_high), width=0.4, colour="black") +
  theme_bw() +
  labs(
    colour="Location: ",
    x="Location",
    y="Measured value"
  ) +
  theme(
    legend.position = "bottom"
  )

ggsave("output/supp_figure_envir_summary.png", width=7, height=8)
ggsave("output/supp_figure_envir_summary.pdf", width=7, height=8)
