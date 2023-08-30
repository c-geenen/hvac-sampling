library(tidyverse)
library(ggtext)
library(rjags)
library(patchwork)

PCR_parameters <- read_csv("input/PCR_parameters.csv")

####################
# Organise environmental variables
##########################

# get mean environmental variables from three room sensors by sampling day.
envir_room_data <- read_csv("input/environmental_data.csv") %>%
  filter(location!="HVAC") %>%
  group_by(samplingDate,variable) %>%
  summarise(
    mean=mean(mean, na.rm=T),
  ) %>%
  pivot_wider(names_from="variable",values_from="mean")

# get HVAC environmental variables
envir_hvac_data <- read_csv("input/environmental_data.csv") %>%
  filter(location=="HVAC",variable!="flow",variable!="velocity") %>%
  mutate(variable=paste0(variable,"_hvac")) %>%
  pivot_wider(names_from="variable",values_from="mean") %>%
  select(samplingDate,co2_hvac,humidity_hvac,temperature_hvac)
  
# combine environmental variables
envir_data <- envir_room_data %>%
  left_join( envir_hvac_data, by="samplingDate" ) %>%
  mutate(
    diff_temperature = temperature_hvac - temperature,
    diff_co2 = co2_hvac - co2,
    diff_humidity = humidity_hvac - humidity
  ) %>%
  select(-temperature_hvac,-co2_hvac,-humidity_hvac)

##############################
# Organise test data
#############################

diff_ct_data <- read_csv("input/test_data.csv") %>%
  
  # only include positive tests (although we will correct for censoring)
  filter(detected) %>%
  
  # only use SARS-CoV-2 N gene
  filter(
    pathogen != "SARS-CoV-2 (respiratory panel)",
    pathogen != "SARS-CoV-2 (TaqPath ORF1ab)",
    pathogen != "SARS-CoV-2 (TaqPath S-protein)"
  ) %>%
    
  pivot_wider(
    names_from = "location",
    values_from = Ct_value
  ) %>%
  
  left_join(envir_data,by="samplingDate") %>%
  
  mutate(
    hvac_position_changed = as_date(samplingDate)>=ymd("2023-03-23"),
    sampling_6h = as_date(samplingDate)>=ymd("2023-03-16")
  )

####################
# Determine censored intervals
####################

# Function to determine lower end of possible Ct value range
# (concentrations below the limit of detection being represented by the equivalent non-observable Ct value)
# If a pathogen is not detected, we assume that the Ct value is somewhere between Ct_LOD and infinity.
Ct_interval_low <- function(Ct, LOD) {
  if_else(
    !is.na(Ct),
    Ct,
    LOD # not detected: Ct is higher than Ct_LOD
  )
}
# Function to determine upper end of Ct value range
Ct_interval_high <- function(Ct, LOD) {
  if_else(
    !is.na(Ct),
    Ct,
    Inf # not detected
  )
}

# Determine intervals of Ct_diff
censored_data <- diff_ct_data %>%
  left_join(PCR_parameters,by="pathogen") %>%

  pivot_longer(c("1","2","3"), names_to = "location", values_to = "room") %>%
  
  filter(
    !is.na(HVAC) | !is.na(room), # detected in any of the two locations
  ) %>%
  
  # calculate range of possible Ct_diff
  mutate(
    censored = is.na(HVAC) | is.na(room),
    left_censored  = if_else(censored, is.na(room), NA), # left censored: not detected in room
    
    Ct_diff = if_else(
      !censored,
      HVAC - room,
      if_else(
        left_censored,
        HVAC - Ct_at_LOD,
        Ct_at_LOD - room
      )
    ),
    
    left_censored = as.integer(left_censored)
  )

#######################
# plot data points with limits
###################

plot_data <- censored_data %>%
  select(pathogen, censored, left_censored, Ct_diff, location) %>%
  mutate(
    limit = if_else(!censored, "known", if_else(left_censored==1, "high", "low")),
    limit = fct_recode(limit, "Left censoring"="high","Right censoring"="low","Observed value"="known"),
    location = as.factor(parse_number(location)),
    Pathogen=pathogen,
    `Room sampler location`=paste0("Location ",location),
    `Overall`= 1
  )
xmin = min(plot_data$Ct_diff, na.rm=T)
xmax = max(plot_data$Ct_diff, na.rm=T)

censored_plot_base <- ggplot(plot_data, aes(Ct_diff, shape=limit, size=limit, colour=location)) +
  theme_bw() +
  coord_cartesian(xlim=c(xmin,xmax)) +
  scale_x_continuous(expand=expansion(add=0)) +
  scale_y_discrete(expand=expansion(add=0.5)) +
  scale_shape_manual(values=c(41,20,40)) +
  scale_size_manual(values=c(3,1,3)) +
  theme(
    axis.text.y = element_markdown(),
    axis.title.x = element_markdown(),
    panel.grid.minor.y = element_blank(),
    panel.grid.major.y = element_blank()
  ) +
  geom_vline(xintercept = 0) +
  labs(x="&Delta;C<sub>t, HVAC - room</sub>", colour="Room sampler location",
       shape="Censoring", size="Censoring")

# by pathogen
(censored_plot_ptn <- censored_plot_base +
  geom_hline(yintercept = seq(0.5,20.5,by=1), alpha=0.2) +
  geom_jitter(aes(y = Pathogen), height=0.2, width=0, alpha=1))

# by location
(censored_plot_loc <- censored_plot_base +
  geom_hline(yintercept = seq(0.5,20.5,by=1), alpha=0.2) +
  geom_jitter(aes(y = `Room sampler location`), height=0.2, width=0, alpha=1))

# no differentiation (all locations and pathogens)
(censored_plot_all <- censored_plot_base +
  geom_jitter(aes(y = `Overall`), height=0.5, width=0, alpha=1)) +
  theme(axis.title.y=element_blank())
  
# combined
(censored_plot_all / censored_plot_ptn / censored_plot_loc) + plot_layout(guides = 'collect', heights=c(1,4,2))

######################
# Alternative visualisation
#####################

plot_data_observed <- filter(plot_data, !censored)
xmin_alt = min(plot_data_observed$Ct_diff, na.rm=T)
xmax_alt = max(plot_data_observed$Ct_diff, na.rm=T)

plot_data_censored <- filter(plot_data, censored) %>%
  mutate(Ct_diff = if_else(
    left_censored==1,
    xmin_alt - 1.25,
    xmax_alt + 1.25
  ))

add_not_detected_label <- function(l) {
  if_else(
    as.numeric(l)==xmin_alt - 1.25,
    str_wrap("Only detected in HVAC", width = 10),
    if_else(
      as.numeric(l)==xmax_alt + 1.25,
      str_wrap("Only detected in room", width = 10),
      as.character(l)
    )
  )
}

censored_plot_base_alt <- ggplot(plot_data_observed, aes(Ct_diff, colour=location)) +
  theme_bw() +
  coord_cartesian(xlim=c(xmin_alt-1.5,xmax_alt+1.5)) +
  scale_y_discrete(expand=expansion(add=0.5)) +
  scale_x_continuous(
    expand=expansion(add=0.5),
    breaks = c(xmin_alt-1.25,xmax_alt+1.25,0,4,-4),
    labels = add_not_detected_label
  ) +
  theme(
    axis.text.y = element_markdown(),
    axis.title.x = element_markdown(),
    panel.grid.minor.y = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.grid.minor.x = element_blank(),
  ) +
  geom_vline(xintercept = 0, alpha=0.5) +
  geom_vline(xintercept = c(xmin_alt-0.5,xmax_alt+0.5)) +
  labs(x="&Delta;C<sub>t, HVAC - room</sub>",colour="Room sampler location",
       shape="Censoring", size="Censoring")

# by pathogen
( censored_plot_ptn_alt <- censored_plot_base_alt +
  geom_hline(yintercept = seq(0.5,20.5,by=1), alpha=1) +
  geom_jitter(aes(y = Pathogen), height=0.2, width=0, alpha=0.5) +
  geom_jitter(data = plot_data_censored, aes(y = Pathogen), height=0.2, width=0.5, alpha=0.5) )
# by location
( censored_plot_loc_alt <- censored_plot_base_alt +
    geom_hline(yintercept = seq(0.5,20.5,by=1), alpha=1) +
    geom_jitter(aes(y = `Room sampler location`), height=0.2, width=0, alpha=0.5) +
    geom_jitter(data = plot_data_censored, aes(y = `Room sampler location`), height=0.2, width=0.5, alpha=0.5) )
# all
( censored_plot_all_alt <- censored_plot_base_alt +
    geom_jitter(aes(y = `Overall`), height=0.5, width=0, alpha=0.5) +
    geom_jitter(data = plot_data_censored, aes(y = `Overall`), height=0.5, width=0.5, alpha=0.5) )
# combined
(censored_plot_all_alt / censored_plot_ptn_alt / censored_plot_loc_alt) + plot_layout(guides = 'collect', heights=c(1,4,2))

#######################
# Model Ct value difference
# Fixed effects: -
# Random effects: pathogen, sampling day, room sampler location
######################

# based on:
# https://bmcbioinformatics.biomedcentral.com/articles/10.1186/s12859-021-04496-8
# https://onlinelibrary.wiley.com/doi/full/10.1111/oik.05985
# https://www.ling.uni-potsdam.de/~vasishth/JAGSStanTutorial/SorensenVasishthMay12014.pdf

jags_data <- censored_data %>%
  filter(!is.na(Ct_diff)) %>%
  arrange(censored)

pathogen_list <- unique(jags_data$pathogen)

# Model with formula:
# Ct_diff ~ intercept + (1|pathogen) + (1|day) + (1|loc)
model_code_all <- 
  "model {
    # fully observed Ct_diff
    for (o in 1:O) {
      Ct_diff[o] ~ dnorm(mu[o], tau)
      mu[o] <- intercept + pn_coeff[pathogen[o]] + day_coeff[day[o]] + loc_coeff[loc[o]]
    }
    # left or right censored Ct_diff
    for (c in 1:C) {
      Z[O+c] ~ dbern(p[c])
      p[c] <- pnorm(Ct_diff[O+c], mu[O+c], tau)
      mu[O+c] <- intercept + pn_coeff[pathogen[O+c]] + day_coeff[day[O+c]] + loc_coeff[loc[O+c]]
    }
    
    # Random effect of day, location, pathogen
    for (d in 1:N_day) {
      day_coeff[d] ~ dnorm(0, tau_day)
    }
    for (j in 1:N_pathogen) {
      pn_coeff[j] ~ dnorm(0, tau_pn)
    }
    for (l in 1:3) {
      loc_coeff[l] ~ dnorm(0, tau_loc)
    }
    
    # Define weakly informative priors for fixed parts
    intercept ~ dnorm(0,0.01)
    sigma ~ dunif(0,100)
    tau <- 1 / (sigma * sigma)
    
    # Define priors for random parts
    tau_day ~ dgamma(1.5, 1.0E-4)
    sigma_day <- pow(tau_day,-1/2)
    tau_pn ~ dgamma(1.5, 1.0E-4)
    sigma_pn <- pow(tau_pn,-1/2)
    tau_loc ~ dgamma(1.5, 1.0E-4)
    sigma_loc <- pow(tau_loc,-1/2)
    
    # Predict the mean Ct_diff
    mu_all <- intercept
    Ct_diff_all <- mu_all
  }"

set.seed(1)
jags_all <- jags.model(
  textConnection(model_code_all),
  data = list(
    'Ct_diff' = jags_data$Ct_diff,
    'pathogen' = match(jags_data$pathogen, pathogen_list),
    'N_pathogen' = length(pathogen_list),
    'O' = sum(jags_data$censored==F),
    'C' = sum(jags_data$censored==T),
    'Z' = jags_data$left_censored,
    'day' = match(jags_data$samplingDate, unique(jags_data$samplingDate)),
    'N_day' = length(unique(jags_data$samplingDate)),
    'loc' = match(jags_data$location, c("1","2","3"))
  ),
  n.chains = 4,
  n.adapt = 10000
)
update(jags_all, 10^4)

# Take samples from the posterior distribution of Ct_diff
samples_Ct_diff_all <- coda.samples( jags_all, c('Ct_diff_all'), 10^5 )
plot(samples_Ct_diff_all)

(s_Ct_diff_all <- summary(samples_Ct_diff_all))
post_dist_all <- as_tibble(t(s_Ct_diff_all$quantiles)) %>%
  mutate(
    `Overall`=1,
    ratio_low = 2^`2.5%`,
    ratio_est = 2^`50%`,
    ratio_high = 2^`97.5%`
  )
post_dist_all
write_csv(post_dist_all,"output/model_all/post_dist.csv")
post_dist_all <- read_csv("output/model_all/post_dist.csv")

plot_all <- censored_plot_all +
  geom_errorbar(
    data=post_dist_all,
    aes(y=`Overall`, x=`50%`, xmin=`2.5%`,xmax=`97.5%`),
    width=1.5, linewidth=0.5, colour="black", alpha=0.7, inherit.aes = F
  ) +
  geom_point(
    data=post_dist_all,
    aes(y=`Overall`, x=`50%`),
    colour="black",inherit.aes = F) +
  theme(
    axis.title.x = element_blank()
  )
plot_all

plot_all_alt <- censored_plot_all_alt +
  geom_errorbar(
    data=post_dist_all,
    aes(y=`Overall`, x=`50%`, xmin=`2.5%`,xmax=`97.5%`),
    width=1.5, linewidth=0.5, colour="black", alpha=0.7, inherit.aes = F
  ) +
  geom_point(
    data=post_dist_all,
    aes(y=`Overall`, x=`50%`),
    colour="black",inherit.aes = F) +
  theme(
    axis.title.x = element_blank(),
    axis.text.x = element_blank()
  )
plot_all_alt

#######################
# Model Ct value difference
# Fixed effects: pathogen
# Random effects: sampling day, room sampler location
######################

# Model with formula:
# Ct_diff ~ intercept + pathogen + (1|day) + (1|loc)
model_code_ptn <- 
  "model {
    # fully observed Ct_diff
    for (o in 1:O) {
      Ct_diff[o] ~ dnorm(mu[o], tau)
      mu[o] <- intercept + pn_coeff[pathogen[o]] + day_coeff[day[o]] + loc_coeff[loc[o]]
    }
    # left or right censored Ct_diff
    for (c in 1:C) {
      Z[O+c] ~ dbern(p[c])
      p[c] <- pnorm(Ct_diff[O+c], mu[O+c], tau)
      mu[O+c] <- intercept + pn_coeff[pathogen[O+c]] + day_coeff[day[O+c]] + loc_coeff[loc[O+c]]
    }
    
    # Random effect of day, loc, occ and dur
    for (d in 1:N_day) {
      day_coeff[d] ~ dnorm(0, tau_day)
    }
    for (l in 1:3) {
      loc_coeff[l] ~ dnorm(0, tau_loc)
    }
    
    # Define weakly informative priors for fixed parts
    pn_coeff[1] <- 0
    for (j in 2:N_pathogen) {
      pn_coeff[j] ~ dnorm(0,0.01) # categorical coefficient prior: SD=10
    }
    intercept ~ dnorm(0,0.01)
    sigma ~ dunif(0,100)
    tau <- 1 / (sigma * sigma)
    
    # Define priors for random parts
    tau_day ~ dgamma(1.5, 1.0E-4)
    sigma_day <- pow(tau_day,-1/2)
    tau_loc ~ dgamma(1.5, 1.0E-4)
    sigma_loc <- pow(tau_loc,-1/2)
    
    ## Predict the mean Ct_diff for each pathogen at 6h sampling duration
    for (pn in 1:N_pathogen) {
      mu_pn[pn] <- intercept + pn_coeff[pn]
      Ct_diff_pn[pn] <- mu_pn[pn] # for actual posterior distribution rather than distrubution of mean: ratio_pn[pn] ~ dnorm(mu_pn[pn],tau)
    }
  }"

# Define and update the model in JAGS
set.seed(2)
jags_ptn <- jags.model(
  textConnection(model_code_ptn),
  data = list(
    'Ct_diff' = jags_data$Ct_diff,
    'pathogen' = match(jags_data$pathogen, pathogen_list),
    'N_pathogen' = length(pathogen_list),
    'O' = sum(jags_data$censored==F),
    'C' = sum(jags_data$censored==T),
    'Z' = jags_data$left_censored,
    'day' = match(jags_data$samplingDate, unique(jags_data$samplingDate)),
    'N_day' = length(unique(jags_data$samplingDate)),
    'loc' = match(jags_data$location, c("1","2","3"))
  ),
  n.chains = 4,
  n.adapt = 10000
)
update(jags_ptn, 10^4)

# Take samples from the posterior distribution of Ct_diff for each pathogen
samples_Ct_diff_pn <- coda.samples( jags_ptn, c('Ct_diff_pn','intercept'), 10^5 )
plot(samples_Ct_diff_pn)

# summarise posterior distributions
(s_Ct_diff_pn <- summary(samples_Ct_diff_pn))
post_dist_ptn <- as_tibble(s_Ct_diff_pn$quantiles, rownames="parameter") %>%
  mutate( Pathogen = c(pathogen_list,NA) )
post_dist_ptn
write_csv(post_dist_ptn,"output/model_pathogen/post_dist.csv")
post_dist_ptn <- read_csv("output/model_pathogen/post_dist.csv")

plot_ptn <- censored_plot_ptn +
  geom_errorbar(
    data=post_dist_ptn,
    aes(y=Pathogen, x=`50%`, xmin=`2.5%`,xmax=`97.5%`),
    width=0.5, linewidth=0.5, colour="black", alpha=0.7, inherit.aes = F
  ) +
  geom_point(data=post_dist_ptn, aes(y=Pathogen, x=`50%`), colour="black",inherit.aes = F) +
  theme(
    axis.title.x = element_blank()
  )
plot_ptn

post_dist_ptn_alt <- post_dist_ptn %>% filter(`2.5%`>xmin_alt, `97.5%`<xmax_alt)
plot_ptn_alt <- censored_plot_ptn_alt +
  geom_errorbar(
    data=post_dist_ptn_alt,
    aes(y=Pathogen, x=`50%`, xmin=`2.5%`,xmax=`97.5%`),
    width=0.5, linewidth=0.5, colour="black", alpha=0.7, inherit.aes = F
  ) +
  geom_point(data=post_dist_ptn_alt, aes(y=Pathogen, x=`50%`), colour="black",inherit.aes = F) +
  theme(
    axis.title.x = element_blank(),
    axis.text.x = element_blank()
  )
plot_ptn_alt

#######################
# Model Ct value difference
# Fixed effect: room sampler location
# Random effects: sampling day, pathogen
######################

# Model with formula:
# Ct_diff ~ intercept + loc + (1|pathogen) + (1|day)
model_code_loc <- 
  "model {
    # fully observed Ct_diff
    for (o in 1:O) {
      Ct_diff[o] ~ dnorm(mu[o], tau)
      mu[o] <- intercept + pn_coeff[pathogen[o]] + day_coeff[day[o]] + loc_coeff[loc[o]]
    }
    # left or right censored Ct_diff
    for (c in 1:C) {
      Z[O+c] ~ dbern(p[c])
      p[c] <- pnorm(Ct_diff[O+c], mu[O+c], tau)
      mu[O+c] <- intercept + pn_coeff[pathogen[O+c]] + day_coeff[day[O+c]] + loc_coeff[loc[O+c]]
    }
    
    # Random effect of day and pathogen
    for (d in 1:N_day) {
      day_coeff[d] ~ dnorm(0, tau_day)
    }
    for (j in 1:N_pathogen) {
      pn_coeff[j] ~ dnorm(0, tau_pn)
    }
    
    # Define weakly informative priors for fixed parts
    loc_coeff[1] <- 0
    for (l in 2:3) {
      loc_coeff[l] ~ dnorm(0,0.01) # categorical coefficient prior: SD=10
    }
    intercept ~ dnorm(0,0.01)
    sigma ~ dunif(0,100)
    tau <- 1 / (sigma * sigma)
    
    # Define priors for random parts
    tau_day ~ dgamma(1.5, 1.0E-4)
    sigma_day <- pow(tau_day,-1/2)
    tau_pn ~ dgamma(1.5, 1.0E-4)
    sigma_pn <- pow(tau_pn,-1/2)
    
    # Predict the mean Ct_diff for each location
    for (k in 1:3) {
      mu_loc[k] <- intercept + loc_coeff[k]
      Ct_diff_loc[k] <- mu_loc[k]
    }
  }"

set.seed(3)
jags_loc <- jags.model(
  textConnection(model_code_loc),
  data = list(
    'Ct_diff' = jags_data$Ct_diff,
    'pathogen' = match(jags_data$pathogen, pathogen_list),
    'N_pathogen' = length(pathogen_list),
    'O' = sum(jags_data$censored==F),
    'C' = sum(jags_data$censored==T),
    'Z' = jags_data$left_censored,
    'day' = match(jags_data$samplingDate, unique(jags_data$samplingDate)),
    'N_day' = length(unique(jags_data$samplingDate)),
    'loc' = match(jags_data$location, c("1","2","3"))
  ),
  n.chains = 4,
  n.adapt = 10000
)
update(jags_loc, 10^4)

# Take samples from the posterior distribution of Ct_diff for each location
samples_Ct_diff_loc <- coda.samples( jags_loc, c('Ct_diff_loc','intercept'), 10^5 )
plot(samples_Ct_diff_loc)

# Summarise posterior distributions
(s_Ct_diff_loc <- summary(samples_Ct_diff_loc))
post_dist_loc <- as_tibble(s_Ct_diff_loc$quantiles, rownames="parameter") %>%
  mutate( `Room sampler location` = c(1:3,NA) )
post_dist_loc
write_csv(post_dist_loc,"output/model_location/post_dist.csv")
post_dist_loc <- read_csv("output/model_location/post_dist.csv")

plot_loc <- censored_plot_loc +
  geom_errorbar(
    data=post_dist_loc,
    aes(y=`Room sampler location`, x=`50%`, xmin=`2.5%`,xmax=`97.5%`),
    width=0.7, linewidth=0.5, colour="black", alpha=0.7, inherit.aes = F
  ) +
  geom_point(
    data=post_dist_loc,
    aes(y=`Room sampler location`, x=`50%`),
    colour="black",inherit.aes = F)
plot_loc

plot_loc_alt <- censored_plot_loc_alt +
  geom_errorbar(
    data=post_dist_loc,
    aes(y=`Room sampler location`, x=`50%`, xmin=`2.5%`,xmax=`97.5%`),
    width=0.7, linewidth=0.5, colour="black", alpha=0.7, inherit.aes = F
  ) +
  geom_point(
    data=post_dist_loc,
    aes(y=`Room sampler location`, x=`50%`),
    colour="black",inherit.aes = F)
plot_loc_alt

############################
# Model Ct value difference, including all environmental variables
# Fixed effects: co2, occupancy, humidity, sampling duration, temperature, HVAC sampler moved
# Random effects: pathogen, sampling day, room sampler location
#########################@

# Model with formula:
# Ct_diff ~ intercept + co2 + occupancy + humidity + samplingDuration + hvacSamplerMoved + (1|location) + (1|pathogen) + (1|day)
model_code_envir <- 
  "model {
    # fully observed Ct_diff
    for (o in 1:O) {
      Ct_diff[o] ~ dnorm(mu[o], tau)
      mu[o] <- intercept + pn_coeff[pathogen[o]] + day_coeff[day[o]] + loc_coeff[loc[o]] + co2_coeff * co2[o] + occ_coeff * occ[o] + hum_coeff * hum[o] + temp_coeff * temp[o] + dur_coeff * dur[o] + moved_coeff * moved[o]
    }
    # left or right censored Ct_diff
    for (c in 1:C) {
      Z[O+c] ~ dbern(p[c])
      p[c] <- pnorm(Ct_diff[O+c], mu[O+c], tau)
      mu[O+c] <- intercept + pn_coeff[pathogen[O+c]] + day_coeff[day[O+c]] + loc_coeff[loc[O+c]] + co2_coeff * co2[O+c] + occ_coeff * occ[O+c] + hum_coeff * hum[O+c] + temp_coeff * temp[O+c] + dur_coeff * dur[O+c] + moved_coeff * moved[O+c]
    }
    
    # Random effect of day
    for (d in 1:N_day) {
      day_coeff[d] ~ dnorm(0, tau_day)
    }
    
    # Define weakly informative priors for fixed parts
    # (ignore some variables, without an assumed causal link, to help convergence)
    pn_coeff[1] <- 0
    for (j in 2:N_pathogen) {
      pn_coeff[j] ~ dnorm(0,0.01) # categorical coefficient prior: SD=10
    }
    loc_coeff[1] <- 0
    for (l in 2:3) {
      loc_coeff[l] ~ dnorm(0,0.01) # categorical coefficient prior: SD=10
    }
    intercept ~ dnorm(0,0.01)
    co2_coeff ~ dnorm(0,0.01)
    occ_coeff ~ dnorm(0,0.01)
    hum_coeff ~ dnorm(0,0.01)
    temp_coeff ~ dnorm(0,0.01)
    dur_coeff ~ dnorm(0,0.01)
    moved_coeff ~ dnorm(0,0.01)
    sigma ~ dunif(0,100)
    tau <- 1 / (sigma * sigma)
    
    # Define priors for random parts
    tau_day ~ dgamma(1.5, 1.0E-4)
    sigma_day <- pow(tau_day,-1/2)
  }"

set.seed(4)
jags_envir <- jags.model(
  textConnection(model_code_envir),
  data = list(
    'Ct_diff' = jags_data$Ct_diff,
    'pathogen' = match(jags_data$pathogen, pathogen_list),
    'N_pathogen' = length(pathogen_list),
    'O' = sum(jags_data$censored==F),
    'C' = sum(jags_data$censored==T),
    'Z' = jags_data$left_censored,
    'co2' = jags_data$co2,
    'occ' = jags_data$occupancy,
    'dur' = jags_data$samplingDuration,
    'hum' = jags_data$humidity,
    'temp' = jags_data$temperature, 
    'moved' = jags_data$hvac_position_changed,
    'day' = match(jags_data$samplingDate, unique(jags_data$samplingDate)),
    'N_day' = length(unique(jags_data$samplingDate)),
    'loc' = match(jags_data$location, c("1","2","3"))
  ),
  n.chains = 4,
  n.adapt = 10000,
)

update(jags_envir, 5*10^5)

# Posterior distribution of each coefficient
samples_coeff_envir <-
  coda.samples( jags_envir,
                c('intercept','occ_coeff','co2_coeff','temp_coeff',
                  'hum_coeff','moved_coeff','dur_coeff'
                  ,'loc_coeff','pn_coeff'
                ), 5*10^6 )

#save(samples_coeff_envir, file="output/model_envir/mcmc_output_envir_5mil.rda")
#load("output/model_envir/mcmc_output_envir_5mil.rda")

# create tibble of model output
samples_envir_tibble <- NULL
for (i in 1:4) {
  samples_envir_tibble <- bind_rows(
    samples_envir_tibble,
    samples_coeff_envir[[i]] %>%
      as_tibble %>%
      mutate(
        chain=as.factor(i),
        iteration=seq(510001:5510000)
      )
  )
}
#write_csv(samples_envir_tibble,"output/model_envir/mcmc_output_envir_5mil.csv")

trace_plot <- ggplot(filter(samples_envir_tibble,iteration%%10==0), aes(x=iteration, y=intercept, colour=chain)) +
  geom_line() +
  theme_bw() +
  labs(colour="Chain", x="Iteration", y="Modelled intercept")
density_plot <- ggplot(filter(samples_envir_tibble,iteration%%10==0), aes(x=intercept, colour=chain)) +
  geom_density() +
  theme_bw() +
  theme(legend.position="none") +
  labs(x="Modelled intercept", y="Density")
trace_plot / density_plot + plot_layout(guides = 'collect')
ggsave("output/model_envir/mcmc_output_envir.png", width=8, height=10)

quantiles <- samples_envir_tibble %>%
  select(-chain,-iteration) %>%
  sapply(quantile, probs = c(0.025, 0.5, 0.975)) %>%
  t
post_dist_envir <- as_tibble(quantiles, rownames="parameter") %>%
  #filter(parameter!="intercept") %>%
  mutate(
    change = case_when(
      parameter=="occ_coeff" ~ "20 more attendees",
      parameter=="hum_coeff" ~ "Relative humidity 20% higher (absolute)",
      parameter=="dur_coeff" ~ "Sampling duration 6h instead of 2h",
      parameter=="co2_coeff" ~ "CO2 concentration 200ppm higher",
      parameter=="temp_coeff" ~ "Temperature 1°C higher",
      parameter=="moved_coeff" ~ "HVAC sampler more central in duct",
      grepl("pn_coeff",parameter) ~ pathogen_list[suppressWarnings(parse_number(parameter))],
      grepl("loc_coeff",parameter) ~ paste0("Location ",suppressWarnings(parse_number(parameter))),
      .default = parameter
    ),
    multiplier = case_match(
      parameter,
      "occ_coeff" ~ 20,
      "hum_coeff" ~ 20,
      "dur_coeff" ~ 4,
      "co2_coeff" ~ 200,
      "temp_coeff" ~ 1,
      "moved_coeff" ~ 1,
      .default = 1
    )
  ) %>%
  mutate_at(vars(matches("%")), ~.*multiplier) %>%
  mutate(change=fct_reorder(change,!grepl("pn_coeff",parameter)))
print(post_dist_envir,n=50)
write_csv(post_dist_envir,"output/model_envir/post_dist.csv")
post_dist_envir <- read_csv("output/model_envir/post_dist.csv")

plot_envir_all <- ggplot(post_dist_envir, aes(`50%`,change, xmin=`2.5%`,xmax=`97.5%`)) +
  geom_errorbar( colour="black") +
  geom_point( colour="black" ) +
  labs(x="Effect on &Delta;Ct", y="Change in variable") +
  theme_bw() +
  theme(
    axis.title.x=element_markdown()
  ) +
  geom_vline(xintercept=0, linetype="dotted") +
  coord_cartesian(xlim=c(-6,6)) 
plot_envir_all

ggsave("output/model_envir/Ct_diff_envir_coefficients_all.png", width=6, height=6)

summary_plot_post_dist <- post_dist_envir %>%
  filter(
    !grepl("pn_coeff",parameter),
    !grepl("loc_coeff",parameter),
    parameter!="intercept"
  )
plot_envir <- ggplot(summary_plot_post_dist, aes(`50%`,change, xmin=`2.5%`,xmax=`97.5%`)) +
  geom_errorbar( colour="black") +
  geom_point( colour="black" ) +
  labs(x="Effect on &Delta;Ct", y="Change in variable") +
  theme_bw() +
  theme(
    axis.title.x=element_markdown()
  ) +
  geom_vline(xintercept=0, linetype="dotted") +
  coord_cartesian(xlim=c(-6,6)) 
plot_envir

ggsave("output/model_envir/Ct_diff_envir_coefficients_summary.png", width=6, height=2)


#############################
# combine plots
#####################@@

(plot_all / plot_ptn / plot_loc / plot_envir)  +
  plot_layout(guides = 'collect', heights=c(1,5,2,2)) & plot_annotation(tag_levels = 'a', tag_suffix=".")

ggsave("output/supp_figure_diffCt.png", width=10, height=10)

(plot_all_alt / plot_ptn_alt / plot_loc_alt / plot_envir)  +
  plot_layout(guides = 'collect', heights=c(1,6,2,2)) & plot_annotation(tag_levels = 'a', tag_suffix=".") & theme(legend.position = 'top')

ggsave("output/figure_4.png", width=8, height=10)

