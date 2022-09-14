#Code for all analyses from Teitelbaum et al.
#Author: Claire Teitelbaum (claire.teitelbaum@gmail.com)
  
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##                      PACKAGES AND SETUP                                  ----
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#load packages
library(dplyr); library(lubridate); library(tidyr); library(stringr)
library(sf)
library(purrr)
library(ggplot2)
library(amt)
library(glmmTMB); library(DHARMa); library(emmeans)
crs <- raster::crs


#ggplot global options
theme_set(theme_classic())
options(ggplot2.continuous.colour="viridis")
options(ggplot2.continuous.fill = "viridis")


##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##                            LOAD DATA                                     ----
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

move_dats_meta <- read.csv("individual_data_attributes.csv") #species, sex, infection, etc.
mig_tracks <- read.csv("spring_migration_data.csv") #movement data for spring migration
local_tracks <- read.csv("local_movement_data.csv") #movement data for local movements

#for birds sampled more than once, add sampling number
move_dats_meta <- mutate(move_dats_meta, sampling_date = ymd(sampling_date)) %>%
  arrange(ID,sampling_date) %>% group_by(ID) %>%
  mutate(sampling_event = 1:length(sampling_date)) %>% ungroup()

#convert to date and get local timestamp for both data sets
mig_tracks <- mutate(mig_tracks, timestamp = ymd_hms(timestamp, tz = "UTC"), 
                     time_local = with_tz(timestamp, tzone = "US/Pacific"))
local_tracks <- mutate(local_tracks, timestamp = ymd_hms(timestamp, tz = "UTC"), 
                       time_local = with_tz(timestamp, tzone = "US/Pacific"))

#add metadata to migration and local movement data
mig_tracks <- left_join(mig_tracks, move_dats_meta) %>% mutate(season = str_replace_all(season," ","_"))






##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##                          MIGRATION TIMING                                ----
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##                ~extract migration departure date                         ----
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#get first and last date of each range/movement season by calendar year
start_dates <- mig_tracks %>%
  #winters stretch across years so date conversion needed: include winter of "2018-2019" as "2019" (ending year)
  mutate(year = ifelse(season=="winter" & month(time_local) > 6, year(time_local) + 1, year(time_local))) %>% 
  group_by(ID,species,season, year) %>% 
  summarize(start = min(time_local), end = max(time_local), segment_id = min(segment_id)) %>% 
  mutate(doy_start = yday(start), doy_end = yday(end)) 

#convert to wide format
#and remove duplicates where a bird was tracked in the previous year too
start_dates_wide <- pivot_wider(start_dates, id_cols=c(ID,species,year),
                                names_from = season, values_from = c(start,end,doy_start,doy_end)) %>%
  arrange(ID, desc(year)) %>% group_by(ID) %>% slice(1) %>% ungroup() 


#calculate spring migration start date from stopover, migration, and breeding start dates
mod_dats_spring_migration <- start_dates_wide %>%
  left_join(move_dats_meta) %>%
  rowwise() %>%
  mutate(start_used =  min(c(start_spring_stopover, start_breed, start_spring_migration), na.rm=T), #earliest start date
         doy_start = yday(start_used), year_start = year(start_used),
         time_since_sampling = difftime(start_used, sampling_date, units = "days"),
         time_since_sampling = as.numeric(time_since_sampling)) %>%
  ungroup() %>%
  mutate(species = factor(species, levels = c("Anser albifrons","Aythya valisineria","Anas acuta")))


#range of start dates and time since sampling
range(mod_dats_spring_migration$doy_start)
range(mod_dats_spring_migration$time_since_sampling)
range(mod_dats_spring_migration$time_since_sampling[mod_dats_spring_migration$species=="Aythya valisineria"])
mean(mod_dats_spring_migration$time_since_sampling)



##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##                    ~extract stopover duration                            ----
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#calculate duration of first stopover for each individual
range_locs <- mig_tracks %>% group_by(segment_id,ID,season,species) %>%
  summarize(duration = difftime(max(time_local), min(time_local), units = "days"),
            year = year(min(time_local)), net_displacement = median(net_displacement),
            .groups = "drop") %>%
  filter(season=="spring_stopover") %>%
  mutate(first_stop = !duplicated(paste(ID,year))) %>%
  left_join(mod_dats_spring_migration %>% select(ID,year_start)) %>%
  filter(year==year_start)

#exclude birds where this stopover is the last point monitored
#because then we are not confident that we have the full duration
last_point <- mig_tracks %>% group_by(ID) %>% summarize(maxRange = max(segment_id)) 
range_locs <- left_join(range_locs, last_point) %>%
  mutate(last_range = segment_id==maxRange)

#get first stopover site info
first_stops <- filter(range_locs, first_stop & !last_range) %>%
  mutate(duration=as.numeric(duration))

#birds with no detected stopovers who should be noted as stopover duration=0: sampling frequency high enough to have detected stopovers if they occurred
#candidates
ids <- mod_dats_spring_migration %>%
  filter(!ID %in% first_stops$ID
         & !ID %in% (range_locs %>% filter(first_stop & last_range) %>% pull(ID))) %>%
  select(ID, end_winter, start_spring_migration, end_spring_migration, start_breed)
#check monitoring gaps
missing_tracks <- mig_tracks %>% filter(ID %in% ids$ID) %>%
  left_join(ids) %>%
  left_join(mod_dats_spring_migration %>% select(ID,year_start)) %>%
  filter(year(time_local)==year_start) %>%
  mutate(dt = difftime(time_local, lag(time_local), units = "days"), dt = ifelse(ID==lag(ID), dt, NA))
ids2 <- missing_tracks %>%
  group_by(ID) %>%
  summarize(dt = max(dt[season=="spring_migration"]),
            end_spring_migration = season[length(season)]=="spring_migration") %>%
  filter(dt < 4 & !end_spring_migration)
missing_tracks %>% filter(ID %in% ids2$ID) %>%
  ggplot(aes(x = time_local, y = net_displacement)) +
  geom_point(aes(color=season)) +
  geom_path() +
  facet_wrap(~ID, scales = "free")


add_zeroes = cbind.data.frame(ID = c("171259_61112","177129_61432","171267_61162","180958.1","180671.1"),
                              species = c("Aythya valisineria","Aythya valisineria","Aythya valisineria","Anser albifrons","Anas acuta"),
                              duration = c(0,0,0,0,0), ND = c(Inf,Inf,Inf,Inf,Inf),
                              year = c(2018,2019,2018,2019,2019))
first_stops = bind_rows(first_stops, add_zeroes)
#"180636.1", "180637.1" excluded because first stop ends up being very near breeding site; unclear if it is a stopover or early arrival
#"177124_41451" excluded: no spring migration points but gap between winter and breed is 9 days

#add this information to already-filtered migration and AIV data. this pre-filters out birds as above
stop_dats <- mod_dats_spring_migration %>% select(year=year_start, ID, 
                                                  active_infection, antibody_detection, species,
                                                  doy_start_mig = doy_start,
                                                  doy_start_breed, doy_start_spring_migration,
                                                  time_since_sampling) %>%
  mutate(year = as.numeric(as.character(year))) %>%
  filter(ID %in% first_stops$ID) %>%
  left_join(first_stops) %>%
  mutate(species = factor(species, levels = levels(mod_dats_spring_migration$species)))



##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##                       ~migration models                                  ----
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#### active infection vs migration start date
#run LMM and LM
mod_dats_pcr_mig <- filter(mod_dats_spring_migration, !is.na(active_infection) & species %in% c("Anas acuta","Aythya valisineria")) 
glmm_pcr_mig <- glmmTMB(doy_start ~ active_infection*species + (1|year_start), data = mod_dats_pcr_mig, na.action = "na.fail", REML = F)
lm_pcr_mig <- lm(doy_start ~ active_infection*species, data = mod_dats_pcr_mig, na.action = "na.fail")
AIC(glmm_pcr_mig); AIC(lm_pcr_mig)
#use LM; lower AIC
mod_pcr_mig <- lm_pcr_mig

#model evaluation
resids <- simulateResiduals(mod_pcr_mig)
testDispersion(mod_pcr_mig)
plot(resids)

#summary tables: coefficients and post-hoc means
summ_mod_pcr <- summary(mod_pcr_mig)

#post-hoc means
means_mod_pcr <- emmeans(mod_pcr_mig, ~ active_infection + species, type="response", cov.reduce = min)
means_mod_pcr %>% emmeans::contrast("pairwise",by="species")
#pairwise comparisons
pairwise_mig_pcr <- means_mod_pcr %>% emmeans::contrast("pairwise",by="species") %>%
  data.frame() %>%
  mutate(CI.low = estimate-1.96*SE, CI.high = estimate+1.96*SE)
#R2
summary(mod_pcr_mig)$r.squared
summary(update(mod_pcr_mig, . ~ . -active_infection - species:active_infection))$r.squared


#### antibodies vs migration start date models
#antibody status only
mod_dats_ab_mig <- filter(mod_dats_spring_migration, !is.na(antibody_detection)) 

glmm_ab_mig <- glmmTMB(doy_start ~ antibody_detection*species + (1|year_start), data = mod_dats_ab_mig, na.action = "na.fail", REML=F)
lm_ab_mig <- lm(doy_start ~ antibody_detection*species, data = mod_dats_ab_mig, na.action = "na.fail")
AIC(glmm_ab_mig); AIC(lm_ab_mig)
mod_ab_mig <- lm_ab_mig

#model evaluation
resids <- simulateResiduals(mod_ab_mig)
testDispersion(mod_ab_mig)
plot(resids)

summ_mod_ab <- summary(mod_ab_mig)
#summary tables: coefficients and post-hoc means
summ_mod_ab <- summ_mod_ab$coefficients %>% round(3) 
means_mod_ab <- emmeans(mod_ab_mig, ~ antibody_detection + species, type="response")
pairwise_mig_ab <- means_mod_ab %>% emmeans::contrast("pairwise",by="species") %>% 
  data.frame() %>%
  mutate(CI.low = estimate-1.96*SE, CI.high = estimate+1.96*SE)
#R2
summary(mod_ab_mig)$r.squared
summary(update(mod_ab_mig, . ~ . -antibody_detection - species:antibody_detection))$r.squared


#### active infection vs stopover models

stop_dats_pcr <- filter(stop_dats, !is.na(active_infection)) 

glmm_stop_duration_pcr <- glmmTMB(duration ~ active_infection*species+doy_start_mig + (1|year), data = stop_dats_pcr, na.action = "na.fail")
lm_stop_duration_pcr <- lm(duration ~ active_infection*species+doy_start_mig, data = stop_dats_pcr, na.action = "na.fail")
AIC(glmm_stop_duration_pcr); AIC(lm_stop_duration_pcr)

mod_stop_duration_pcr <- lm_stop_duration_pcr
summ_stop_pcr <- summary(mod_stop_duration_pcr)$coefficients %>% round(3) 

resids <- simulateResiduals(mod_stop_duration_pcr)
plotResiduals(resids)
means_stop_pcr <- emmeans(mod_stop_duration_pcr, ~ active_infection + species, type="response")
pairwise_stop_pcr <- means_stop_pcr %>% emmeans::contrast("pairwise",by="species") %>% 
  data.frame() %>%
  mutate(CI.low = estimate-1.96*SE, CI.high = estimate+1.96*SE)
as.data.frame(means_stop_pcr) %>% mutate_at(-c(1:2), round, 3) 


#### antibody status vs stopover models
#antibody status
stop_dats_ab <- stop_dats %>% filter(!is.na(antibody_detection))
glmm_stop_duration_ab <- glmmTMB(duration ~ antibody_detection*species+doy_start_mig + (1|year), data = stop_dats_ab)
lm_stop_duration_ab <- lm(duration ~ antibody_detection*species+doy_start_mig, data = stop_dats_ab, na.action = "na.fail")
AIC(glmm_stop_duration_ab); AIC(lm_stop_duration_ab)

mod_stop_duration_ab <- lm_stop_duration_ab

resids <- simulateResiduals(mod_stop_duration_ab)
plotResiduals(resids)

summ_stop_ab <- summary(mod_stop_duration_ab)$coefficients %>% round(3) 

means_stop_ab <- emmeans(mod_stop_duration_ab, ~ antibody_detection + species, type="response")
pairwise_stop_ab <- means_stop_ab %>% emmeans::contrast("pairwise",by="species") %>% 
  data.frame() %>%
  mutate(CI.low = estimate-1.96*SE, CI.high = estimate+1.96*SE)

as.data.frame(means_stop_ab) %>% mutate_if(is.numeric, round, 3) 

#predictions for each group
expand.grid(antibody_detection = c(T,F), species = unique(stop_dats_ab$species), doy_start_mig = 73) %>%
  mutate(pred = predict(lm_stop_duration_ab, .),
         se = predict(lm_stop_duration_ab, ., se.fit = T)$se.fit,
         lower = pred-1.96*se, upper = pred+1.96*se)



##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##                ~migration model summaries and plots                      ----
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#### table of pairwise comparisons from all models
pairwise_full <- list("mig_pcr" = pairwise_mig_pcr, "mig_ab" = pairwise_mig_ab,
                      "stop_pcr" = pairwise_stop_pcr, "stop_ab" = pairwise_stop_ab) %>%
  bind_rows(.id = "model") %>%
  mutate(influenza_type = str_split_fixed(model,"_",2)[,2],
         influenza_type = ifelse(influenza_type == "pcr", "Active infection (rRT-PCR)", "Antibody (bELISA)"),
         model = str_split_fixed(model,"_",2)[,1],
         model = ifelse(model == "mig", "Migration start date", "Stopover duration")) %>%
  mutate_at(vars(estimate,t.ratio,CI.low,CI.high), ~ifelse(contrast == "FALSE - TRUE", -1*.x, .x)) %>%
  mutate_if(is.numeric, round, 3) %>%
  mutate(estimate = paste0(estimate," [",CI.low,", ",CI.high,"]")) %>%
  left_join(distinct(mod_dats_spring_migration,species,species)) %>%
  select(Model=model,species=species,`Influenza metric`=influenza_type,`Difference (positive-negative)`=estimate,SE,df,t=t.ratio,p=p.value) 



#### migration results plot, fig.width=7.5, fig.height = 4.5}
species_labels <- c("Anas acuta"="Northern pintail",
                    "Anser albifrons"="Gr. white-fronted goose",
                    "Anas platyrhynchos"="Mallard",
                    "Aythya valisineria"="Canvasback")

plot_dats_mig <- bind_rows(
  mod_dats_ab_mig %>% mutate(status = antibody_detection, type = "ab"),
  mod_dats_pcr_mig %>% mutate(status = active_infection, type = "pcr")) %>% 
  mutate(type_name = ifelse(type == "ab", "Antibodies", "Active infection"),
         species = factor(species, levels = c("Aythya valisineria","Anas acuta","Anser albifrons")),
         species = forcats::fct_relabel(species, function(x) species_labels[x]))
mod_dats_mig <- bind_rows(
  data.frame(means_mod_pcr) %>% mutate(status = active_infection, type = "pcr"), 
  data.frame(means_mod_ab) %>% mutate(status = antibody_detection, type = "ab")) %>%
  left_join(distinct(plot_dats_mig,species,species)) %>%
  mutate(type_name = ifelse(type == "ab", "Antibodies", "Active infection"),
         species = factor(species, levels = c("Aythya valisineria","Anas acuta","Anser albifrons")),
         species = forcats::fct_relabel(species, function(x) species_labels[x]))


plot_dats_stop <- bind_rows(
  stop_dats %>% filter(!is.na(antibody_detection)) %>% mutate(status = antibody_detection, type = "ab"),
  stop_dats %>% filter(!is.na(active_infection)) %>% mutate(status = active_infection, type = "pcr")) %>%
  mutate(type_name = ifelse(type == "ab", "Antibodies", "Active infection"),
         species = factor(species, levels = c("Aythya valisineria","Anas acuta","Anser albifrons")),
         species = forcats::fct_relabel(species, function(x) species_labels[x]))
mod_dats_stop <- bind_rows(
  as.data.frame(means_stop_pcr) %>% 
    mutate(status = active_infection, type = "pcr",), 
  as.data.frame(means_stop_ab) %>% 
    mutate(status = antibody_detection, type = "ab")
) %>%
  left_join(distinct(plot_dats_stop,species,species))%>%
  mutate(type_name = ifelse(type == "ab", "Antibodies", "Active infection"),
         species = factor(species, levels = c("Aythya valisineria","Anas acuta","Anser albifrons")),
         species = forcats::fct_relabel(species, function(x) species_labels[x]))


labs_mig <- plot_dats_mig %>% bind_rows(mod_dats_mig %>% mutate(doy_start = upper.CL)) %>%
  group_by(species) %>% summarize(doy_start = max(doy_start, na.rm=T)) %>%
  expand_grid(type_name = unique(plot_dats_mig$type_name)) %>%
  arrange(type_name, species) %>%
  filter(!(species == "Gr. white-fronted goose" & type_name == "Active infection")) %>%
  mutate(label = LETTERS[1:5], status = "FALSE")
samp_size_mig <- plot_dats_mig %>% bind_rows(mod_dats_mig %>% mutate(doy_start = upper.CL)) %>%
  group_by(species) %>% summarize(doy_start = min(doy_start, na.rm=T)-11) %>%
  left_join(plot_dats_mig %>% group_by(species,type_name, status) %>% summarize(n = n())) 


plot_mig = ggplot(plot_dats_mig, aes(x = status, y = doy_start, color = status, shape = status)) +
  geom_point(alpha = 0.3, position = position_jitter(width = 0.08, seed = 2)) +
  facet_grid(rows = vars(species), cols = vars(type_name), scales = "free_y", switch = "both", 
             labeller = label_wrap_gen(width=20)) +
  scale_x_discrete("", breaks = c("TRUE","FALSE"), labels = c("Pos.","Neg.")) +
  geom_segment(data = mod_dats_mig, aes(x=status,xend=status,y=lower.CL,yend=upper.CL)) +
  geom_point(data = mod_dats_mig, aes(x=status,y=emmean), size= 2.5) +
  scale_y_continuous("Spring migration start date", 
                     breaks = c(1,32,60,91,121),
                     labels = as.Date(c(1,32,60,91,121), origin = "2012-12-31") %>% format("%d %b") %>% str_remove("0")) +
  theme(legend.position = "none", plot.title = element_text(hjust=0.5), panel.spacing.x = unit(0,"line"),
        strip.background.x = element_blank(), strip.placement.x = "outside",
        strip.text.x = element_text(size = 11, vjust = 1),
        panel.border = element_rect(fill=NA, color = "black")) +
  scale_color_brewer(palette="Dark2") +
  ggtitle("Spring migration start date") +
  geom_text(data = labs_mig, aes(x = status, y = doy_start, label = label), color = "black",
            hjust = 0, vjust = 1, nudge_x = -0.5) +
  geom_text(data = samp_size_mig, aes(x = status, y = doy_start, label = n), color = "black",
            hjust = 0.5, vjust = 0.3, size = 3, fontface = 3) +
  geom_text(data = samp_size_mig %>% filter(!status & species == "Canvasback"),
            aes(x = status, y = doy_start, label = "N="), color = "black",
            hjust = 2.2, vjust = 0.3, size = 3, fontface = 3)



labs_stop <- plot_dats_stop %>% bind_rows(mod_dats_stop %>% mutate(duration = upper.CL)) %>%
  group_by(species) %>% summarize(duration = max(duration, na.rm=T)) %>%
  expand_grid(type_name = unique(plot_dats_mig$type_name)) %>%
  arrange(type_name, species) %>%
  filter(!(species == "Gr. white-fronted goose" & type_name == "Active infection")) %>%
  mutate(label = LETTERS[6:10], status = "FALSE")
samp_size_stop <- plot_dats_stop %>% bind_rows(mod_dats_stop %>% mutate(duration = upper.CL)) %>%
  group_by(species) %>% summarize(duration = min(duration, na.rm=T)-8) %>%
  left_join(plot_dats_stop %>% group_by(species,type_name, status) %>% summarize(n = n())) 

plot_stop = ggplot(plot_dats_stop, aes(x = status, y = duration, color = status, shape = status)) +
  geom_point(alpha = 0.3, position = position_jitter(width = 0.08, seed = 12)) +
  facet_grid(rows = vars(species), cols = vars(type_name), scales = "free_y", switch = "x", 
             labeller = label_wrap_gen(width=20)) +
  scale_y_continuous("Duration of first stopover (days)", breaks = seq(0,60,15)) +
  scale_x_discrete("", breaks = c("TRUE","FALSE"), labels = c("Pos.","Neg.")) +
  geom_segment(data = mod_dats_stop, aes(x=status,xend=status,y=lower.CL,yend=upper.CL)) +
  geom_point(data = mod_dats_stop, aes(x=status,y=emmean), size= 2.5) +
  theme(panel.background = element_rect(size=1), panel.border = element_rect(fill = NA), plot.title = element_text(hjust=0.5))+
  theme(legend.position = "none", plot.title = element_text(hjust=0.5), panel.spacing.x = unit(0,"line"),
        strip.background.x = element_blank(), strip.placement.x = "outside",
        strip.text.x = element_text(size = 11, vjust = 1),
        panel.border = element_rect(fill=NA, color = "black")) +
  scale_color_brewer(palette="Dark2") +
  ggtitle("Duration of first stopover") +
  geom_text(data = labs_stop, aes(x = status, y = duration, label = label), color = "black",
            hjust = 0, vjust = 1, nudge_x = -0.5) +
  geom_text(data = samp_size_stop, aes(x = status, y = duration, label = n), color = "black",
            hjust = 0.5, vjust = 0.3, size = 3, fontface = 3) +
  geom_text(data = samp_size_stop %>% filter(!status & species == "Canvasback"), 
            aes(x = status, y = duration, label = "N="), color = "black",
            hjust = 2.2, vjust = 0.3, size = 3, fontface = 3) 


gridExtra::grid.arrange(plot_mig, plot_stop, ncol=2, widths = c(1,1))


#### migration tracks plots 
#line plots of migration timing separated by species and infection status
full_dat <- bind_rows(mod_dats_pcr_mig %>% mutate(year = as.numeric(as.character(year_start)), status = active_infection, type_name = "Active infection"),
                      mod_dats_ab_mig %>% mutate(year = as.numeric(as.character(year_start)), status = antibody_detection, type_name = "Antibodies") %>% 
                        mutate(year = as.numeric(as.character(year_start)))) %>%
  select(ID, status, type_name, start_used, year_start, species) %>%
  mutate(start_used = str_replace_all(start_used, "2015|2016|2017|2018|2019|2021","2020"),
         start_used= ymd_hms(start_used),
         species = factor(species, levels = c("Aythya valisineria","Anas acuta","Anser albifrons")),
         species = forcats::fct_recode(species, "Gr. white-fronted goose"="Anser albifrons"),
         type_name = factor(type_name, levels = c("Active infection","Antibodies")),
         transparency = ifelse((status & type_name == "Active infection" | !status & type_name == "Antibodies"),
                               1, 0.7))

#get corresponding migration info for each individual
tracks_plot <- left_join(full_dat, mig_tracks %>% select(-species)) %>%
  filter(year(time_local)==year_start & yday(time_local)<190) %>%
  mutate(time_local2 = str_replace_all(time_local, "2015|2016|2017|2018|2019|2021","2020"),
         time_local2 = ymd_hms(time_local2)) 

trackplot_labs <- distinct(tracks_plot, species, type_name) %>% arrange(type_name,species) %>% 
  mutate(lab = LETTERS[1:nrow(.)]) 

track_plot <- tracks_plot %>% 
  filter(time_local2 <= "2020-07-01" & time_local2 >= "2020-01-18") %>% 
  ggplot(aes(x = time_local2, y = net_displacement/1000, color = status)) +
  geom_path(aes(group = ID, alpha = transparency), size=1) +
  lemon::facet_rep_grid(rows = vars(species), cols = vars(type_name), switch = "y") +
  scale_x_datetime(date_labels = "%d %b", date_breaks = "2 months") +
  ylab("Distance from wintering site (km)") +
  scale_color_brewer("Infection status", palette = "Dark2", breaks = c(T,F), labels= c("Pos.","Neg.")) +
  scale_alpha(range = c(0.5,1), guide = "none") +
  geom_rug(data = full_dat, aes(x = start_used + runif(nrow(full_dat),0,1), color = status),
           inherit.aes = F, size = 0.3, length = unit(0.04, "npc")) +
  theme(legend.position = "bottom", plot.title = element_text(hjust=0.5), panel.spacing.x = unit(0,"line"),
        strip.background.x = element_blank(), strip.placement.x = "outside",
        strip.text.x = element_text(size = 11, vjust = 1), 
        panel.border = element_rect(fill=NA, color = "black"), axis.title.x = element_blank()) +
  geom_text(data = trackplot_labs, aes(x = as.POSIXct("2020-01-08"), y = max(tracks_plot$net_displacement/1000), label = lab), color = "black",
            nudge_y = -100)

track_plot





##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##                            LOCAL MOVEMENTS                               ----
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##            ~prep data and calculate movement metrics                     ----
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#### add sunrise times to local movement data
#add sunrise and sunset times, whether fix is during daytime
sunTimes = local_tracks %>%
  mutate(date = as_date(time_local), lat = latitude, lon = longitude) %>%
  select(date,lat,lon) %>%
  suncalc::getSunlightTimes(data = ., tz = "US/Pacific")

mcp_dats <- local_tracks %>% mutate(sunrise = sunTimes$sunrise, sunset = sunTimes$sunset) %>%
  mutate(timeToSunrise = difftime(time_local, sunrise, units = "hours"),
         timeToSunset = difftime(time_local, sunset, units = "hours"),
         rowid = 1:nrow(.))

#add equal-distance projected coordinates
coords_proj = mcp_dats %>% st_as_sf(coords = c("longitude","latitude"), crs=crs("EPSG:4326")) %>%
  st_transform(crs("+proj=aeqd")) %>% st_coordinates()
mcp_dats <- mutate(mcp_dats, y_proj = coords_proj[,2], x_proj = coords_proj[,1])

#nest by ID
#use amt to create tracks
mcp_dats <- mcp_dats %>% 
  nest(track = -c(ID)) %>%
  mutate(track = map(track, ~make_track(.x, .x = x_proj, .y = y_proj, .t = time_local,
                                        rowid = rowid, 
                                        timeToSunrise = timeToSunrise, timeToSunset = timeToSunset)))


#### calculate MCPs
mcp_dats <- mcp_dats %>%
  mutate(track = map(track, ~mutate(.x, date = as_date(t_),
                                    date_sunrise = if_else(timeToSunrise < 0, date-1, date)))) %>% #change date to be split at sunrise 
  mutate(data = map(track, ~nest(.x, data = -c(date_sunrise)))) %>% #within each individual, each day gets a separate df
  unnest(data) %>%
  #filter to days with at least 6 fixes (at least 5 needed to fit an MCP)
  mutate(n_fix = map_int(data, nrow), 
         n_hours = map_int(data, ~length(unique(hour(.x$t_)))),
         timerange_hours = map_dbl(data, ~difftime(max(.x$t_), min(.x$t_), units = "hours"))) %>%
  filter(n_fix >= 6)

#calculate MCPs and distance metrics
mcps <-  mutate(mcp_dats, 
                step_lens = map(data, step_lengths),
                mcp = map(data, ~hr_mcp(.x, levels = 1)),
                distMat = map(data, ~dist(.x[,c("x_","y_")])),
                #get maximum pairwise distance, maximum distance from start point, and MCP area
                total_dist = map_dbl(step_lens, ~sum(.x, na.rm=T)),
                max_dist = map_dbl(distMat, function(x) max(x)),
                max_dist_start = map_dbl(distMat, function(x) max(x[1:(nrow(x)-1)])),
                mcp_area = map_dbl(mcp, ~hr_area(.x)$area)) %>%
  left_join(move_dats_meta) %>%
  mutate(sampling_date = ymd(sampling_date),
         time_since_sampling = difftime(date_sunrise, sampling_date, units = "days"))

mcps_df <- select(mcps, -c(track,data,mcp,distMat)) 



##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##            ~filter data and check correlations                           ----
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#filter to include only the second sampling event per individual if multiple were included
dats_mcp_between <- arrange(mcps_df, desc(sampling_event)) %>%
  distinct(ID, date_sunrise, .keep_all=T) %>%
  filter(time_since_sampling <= 12 & time_since_sampling > 0)  %>%
  mutate(mcp_area = mcp_area/(1000^2))

#where does the relationship between # hours and mcp area asymptote?
dats_mcp_between %>% select(n_hours, mcp_area, total_dist, max_dist) %>%
  pivot_longer(-n_hours, names_to = "metric", values_to = "value") %>%
  ggplot(aes(x = n_hours, y = value)) +
  geom_point(alpha = 0.5) + 
  stat_summary(color = "red", geom = "point", fun = "mean") +
  geom_smooth() + geom_vline(xintercept = 6) +
  scale_x_continuous("Number of hours in fixes", trans = "log10") + 
  scale_y_continuous("",trans = "log10") +
  facet_wrap(~metric, scales = "free", strip.position = "left")


#mcp area vs pcr infection status: replicates per day
dats_mcp_between_pcr <- filter(dats_mcp_between, n_hours >= 6) %>%
  mutate(time_since_sampling = as.integer(time_since_sampling))

sample_sizes <- group_by(dats_mcp_between_pcr, ID) %>% summarize(n = n(), first3 = any(0:3 %in% time_since_sampling))
# kept <- filter(sample_sizes, n>=3 & first3) %>% pull(ID)
kept <- filter(sample_sizes, n>=3) %>% pull(ID)
# kept <- unique(sample_sizes$ID)

dats_mcp_between_pcr <- filter(dats_mcp_between_pcr, ID %in% kept) %>%
  mutate(species = factor(species, levels = c("Aythya valisineria","Anas acuta", "Anas platyrhynchos")))

#check correlations among response variables and with sampling intensity
dats_mcp_between_pcr %>%
  mutate(log_fix = log(n_fix)) %>%
  mutate_at(vars(mcp_area, total_dist, max_dist), log) %>%
  select(mcp_area,total_dist,max_dist,n_fix,log_fix,n_hours) %>%
  cor()



##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##                         ~plots of space use                              ----
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#### plot mcp area vs time since sampling
labs_pts <- bind_cols(species = sort(unique(dats_mcp_between_pcr$species)),
                      y = max(dats_mcp_between_pcr$mcp_area),
                      x = min(dats_mcp_between_pcr$time_since_sampling),
                      label = c("A","B","C"))
labs <- dats_mcp_between_pcr %>% group_by(species, time_since_sampling,active_infection) %>%
  summarize(mn = mean(log(mcp_area)), sd = sd(log(mcp_area),na.rm=T)) %>%
  group_by(species) %>% summarize(x = min(time_since_sampling), y = max(exp(mn+sd),na.rm=T)) %>%
  mutate(label = LETTERS[1:nrow(.)])

axis_bks <- c(.001,.01,.1,1,10,100,1000,10000)
axis_ticks <- bind_cols(start = axis_bks, end = lead(axis_bks)) %>%
  filter(!is.na(end)) %>%
  apply(1, function(x) seq(x[1],x[2],x[1])) %>%
  c() %>% unique()
axis_labs <- ifelse(axis_ticks %in% axis_bks, axis_ticks, "")

mcp_plot_raw <- ggplot(dats_mcp_between_pcr, aes(x = time_since_sampling, y = mcp_area, color = active_infection, shape = active_infection), alpha = 0.7) +
  stat_summary(geom = "line", fun = "mean", position = position_dodge(width = 0.2)) +
  stat_summary(fun.data = "mean_sdl", fun.args = list(mult=1), size = 0.5, position = position_dodge(width = 0.2)) +
  facet_wrap(~species, scales = "free", labeller = labeller(species=species_labels)) + 
  scale_x_continuous("Days since sampling", breaks = seq(2,12,3)) +
  scale_y_continuous("MCP area (sq. km)", trans = "log", 
                     breaks = axis_ticks,
                     labels = axis_labs) + 
  scale_color_manual("Active infection\nstatus (rRT-PCR)", breaks = c(T,F), labels = c("Pos.","Neg."),
                     values = RColorBrewer::brewer.pal(3,"Dark2")[2:1]) +
  scale_shape("Active infection\nstatus (rRT-PCR)", breaks = c(T,F), labels = c("Pos.","Neg.")) +
  theme(axis.ticks.y = element_line(size = 0.4))

mcp_plot_raw +
  geom_text(data = labs, aes(x = x,y = y, label = label), inherit.aes= F)


#### plots for other movement metrics, fig.width=8, fig.height=3.5}
max_dist_plot_raw <- ggplot(dats_mcp_between_pcr, aes(x = time_since_sampling, y = max_dist/1000, color = active_infection, shape = active_infection), alpha = 0.7) +
  stat_summary(geom = "line", fun = "mean", position = position_dodge(width = 0.2)) +
  stat_summary(fun.data = "mean_sdl", fun.args = list(mult=1), size = 0.5, position = position_dodge(width = 0.2)) +
  facet_wrap(~species, scales = "free", labeller = labeller(species=species_labels)) + 
  scale_x_continuous("Days since sampling", breaks = seq(2,12,3)) +
  scale_y_continuous("Maximum daily distance moved (km)", trans = "log", 
                     breaks = axis_ticks,
                     labels = axis_labs) + 
  scale_color_manual("Active infection\nstatus (rRT-PCR)", breaks = c(T,F), labels = c("Pos.","Neg."),
                     values = RColorBrewer::brewer.pal(3,"Dark2")[2:1]) +
  scale_shape("Active infection\nstatus (rRT-PCR)", breaks = c(T,F), labels = c("Pos.","Neg.")) +
  theme(axis.ticks.y = element_line(size = 0.4))
max_dist_plot_raw

total_dist_plot_raw <- ggplot(dats_mcp_between_pcr, aes(x = time_since_sampling, y = total_dist/1000, color = active_infection, shape = active_infection), alpha = 0.7) +
  stat_summary(geom = "line", fun = "mean", position = position_dodge(width = 0.2)) +
  stat_summary(fun.data = "mean_sdl", fun.args = list(mult=1), size = 0.5, position = position_dodge(width = 0.2)) +
  facet_wrap(~species, scales = "free", labeller = labeller(species=species_labels)) + 
  scale_x_continuous("Days since sampling", breaks = seq(2,12,3)) +
  scale_y_continuous("Total daily path length (km)", trans = "log", 
                     breaks = axis_ticks,
                     labels = axis_labs) + 
  scale_color_manual("Active infection\nstatus (rRT-PCR)", breaks = c(T,F), labels = c("Pos.","Neg."),
                     values = RColorBrewer::brewer.pal(3,"Dark2")[2:1]) +
  scale_shape("Active infection\nstatus (rRT-PCR)", breaks = c(T,F), labels = c("Pos.","Neg.")) +
  theme(axis.ticks.y = element_line(size = 0.4))
total_dist_plot_raw



##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##                         ~models of space use                             ----
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#### model mcp area vs time since sampling and infection

mod_mcp_between_pcr <- glmmTMB(log(mcp_area) ~ active_infection + species + log(time_since_sampling) + sex + 
                                 active_infection:species + sex:species +
                                 active_infection:log(time_since_sampling) + species:log(time_since_sampling) + 
                                 active_infection:log(time_since_sampling):species +
                                 log(n_fix) +
                                 (1|ID), data = dats_mcp_between_pcr, na.action = "na.fail")


#model diagnostics
resids <- simulateResiduals(mod_mcp_between_pcr)
par(mfrow = c(2,2))
plotResiduals(resids, log(dats_mcp_between_pcr$time_since_sampling))
plotResiduals(resids, factor(dats_mcp_between_pcr$species))
plotResiduals(resids, factor(dats_mcp_between_pcr$active_infection))

#model summary
summ_mcp_mod <- summary(mod_mcp_between_pcr)
performance::r2_nakagawa(mod_mcp_between_pcr)

#predictions for day 1 
emmeans(mod_mcp_between_pcr, ~ active_infection + species + time_since_sampling, type="response", cov.reduce = min)
#predictions for day 12
emmeans(mod_mcp_between_pcr, ~ active_infection + species + time_since_sampling, type="response", cov.reduce = max)

#get R2 reduction from removing infection status and its interactions from model
mod_mcp_pcr_only <- glmmTMB(log(mcp_area) ~ active_infection +
                              active_infection:species + 
                              active_infection:log(time_since_sampling) + 
                              active_infection:log(time_since_sampling):species +
                              (1|ID), data = dats_mcp_between_pcr, na.action = "na.fail")
performance::r2_nakagawa(mod_mcp_pcr_only)


mod_mcp_no_pcr <- glmmTMB(log(mcp_area) ~ species + log(time_since_sampling) + Sex + 
                            Sex:species + species:log(time_since_sampling) + 
                            log(n_fix) +
                            (1|ID), data = dats_mcp_between_pcr, na.action = "na.fail")
performance::r2_nakagawa(mod_mcp_no_pcr)


#### models for other movement metrics
#model of total distance moved
mod_total_dist_between_pcr <- update(mod_mcp_between_pcr, log(total_dist) ~ .)

#model of maximum daily distance between points
mod_max_dist_between_pcr <- update(mod_mcp_between_pcr, log(max_dist) ~ .)




##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##                             SAMPLE SIZES                                 ----
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
sample_size_by_analysis <- list(mig_ab = mod_dats_ab_mig, mig_pcr = mod_dats_pcr_mig,
                                stopover_ab = stop_dats_ab, stopover_pcr = stop_dats_pcr,
                                mcp_pcr = dats_mcp_between_pcr %>% select(ID,active_infection,species)) %>%
  bind_rows(.id = "model") %>% mutate(model = factor(model, levels = unique(model))) %>%
  select(ID,species,active_infection,antibody_detection,model) %>% distinct() %>%
  group_by(species,model) %>%
  summarize(n = n(),
            pcr_pos = sum(active_infection,na.rm=T), pcr_tot = sum(!is.na(active_infection)),
            pcr_prev = pcr_pos/pcr_tot,
            ab_pos = sum(antibody_detection,na.rm=T), ab_tot = sum(!is.na(antibody_detection)),
            ab_prev = ab_pos/ab_tot) %>%
  mutate(n_pos = ifelse(str_detect(model,"pcr"), pcr_pos, ab_pos), n_neg = n-n_pos,
         species = factor(species),
         species = forcats::fct_relevel(species,"Anser albifrons"),
         pcr_ab = str_detect(model,"pcr"))


labs <- c("Local movements,\nActive infection",
          paste("Migration start",c("Antibodies","Active infection"),sep=",\n"),
          paste("Stopover",c("Antibodies","Active infection"),sep=",\n"))
names(labs) <- c("mcp_pcr","mig_ab","mig_pcr","stopover_ab","stopover_pcr")
labs <- labs[levels(sample_size_by_analysis$model)]

ss_wide <- sample_size_by_analysis %>%
  select(model,pcr_ab,species,n,n_pos,n_neg) %>%
  pivot_longer(cols=c(n_pos,n_neg)) %>%
  mutate(lab = labs[model], lab = factor(lab,levels=labs),
         species = species_labels[as.character(species)],
         species = factor(species, levels = species_labels))
sample_size_by_analysis <- mutate(sample_size_by_analysis, lab = labs[model], lab = factor(lab,levels=labs),
                                  species = species_labels[as.character(species)],
                                  species = factor(species, levels = species_labels))

#dot/pie chart
ggplot(ss_wide, aes(x=n/2,y=value,fill=paste(name,pcr_ab),width=n)) +
  geom_bar(stat="identity",position="fill", color="black") +
  coord_polar("y", clip="off") + 
  geom_text(data=sample_size_by_analysis,
            aes(label = paste(n_pos,n,sep="/"),x=n,y=0), vjust = -0.5, 
            inherit.aes = F, position = position_stack(vjust = -2)) +
  facet_grid(rows=vars(species),cols=vars(lab), switch = "both",
             labeller = label_wrap_gen(width=20)) +
  scale_fill_manual("AIV sample type", breaks = paste(c("n_pos","n_neg","n_pos","n_neg"),c(T,T,F,F)), 
                    values = c("blue","lightgrey","red","lightgrey")) +
  theme_void() + 
  theme(legend.position = "none", 
        strip.text.y.left = element_text(size=11, angle = 0, hjust = 1),
        strip.text.x.bottom = element_text(size=11, angle = 90, hjust = 1))





