# Script information ----
## This includes the main analyses for Taff et al. tree swallow organismal resilience
## Written by Conor Taff: cct663@gmail.com
## Last updated 7/8/2025

# Load packages & settings ----
    # Load packages used in code. These may need to be installed the first time
        pacman::p_load(plyr, lubridate, reshape2, here, tidyverse, lme4, sjPlot, 
                       mgcv, gamm4, gratia, data.table, brms, colorspace, rptR,
                       DHARMa, ggeffects, glmmTMB, ggpubr, Hmisc, viridis)

    # color settings to be used throughout
        fem_color <- "#D55E00"
        mal_color <- "#009E73"
        hb_color <- "#56B4E9"
        fw_color <- "#e1BE6A"
        mg_color <- "#882255"
        
        off_color <- "slateblue"
        on_color <- "coral3"
        cort_color <- "gray50"
        
        robust_color <- "#009292"
        resil_color <- "#B66DFF"
        diff_color <- "#924900"
        
    # function to convert logit to probability
        inv_logit <- function(x){
          exp(x) / (1 + exp(x))
        }

# Load Data ----
    # these are either raw or produced by other scripts and saved
    
    # incubation data from hobos
        all_bouts_ex <- readRDS(here::here("processed_data_output/all_bouts_ex.rds"))
    # provisioning data from rfids
        feeds <- read.delim(here::here("processed_data_output", "feed_per_hour.txt"))
    # capture data from banding database
        captures <- read.delim(here::here("raw_data", "data_by_capture.txt"))
        ith_capture <- captures[captures$location == "Ithaca", ] # subset to Ithaca removing AK/TN/WY
    # nest data from banding database
        nests <- read.delim(here::here("raw_data", "data_by_nest.txt"))
    # list of treatments indicating whether to exclude or not
        exclude_trt <- read.csv(here::here("raw_data", "excluded_treatments.csv"))
        

# Load and wrangle temperature data ----
    # Read in and process temperature data. 
      # This is from the Ithaca Airport Weather Station from 1996-2024
        # provided by Northeast Regional Climate Center (NRCC)
        # load then wrangle columns as needed
              dweath <- read.delim(here::here("raw_data", "ithaca_airport_weather.txt"), 
                                   header = TRUE, stringsAsFactors = FALSE)  
              dweath$date <- as.Date(dweath$date, format = "%m/%d/%y")
              dweath$yday <- yday(dweath$date)
              dweath$temp_C <- (dweath$temp_F - 32) * 5/9
              dweath$year <- year(dweath$date)
              dweath$yr_day_hr <- paste(dweath$year, dweath$yday, dweath$hour, sep = "_")
        
        
        # And this is from the Game Farm Road weather station prior to 1996
            # load then wrangle as needed
              dgfr <- read.delim(here::here("raw_data", "GFR_weather.txt"),
                                 header = TRUE, stringsAsFactors = FALSE)
              dgfr <- subset(dgfr, dgfr$Temp > -50)
              dgfr$date <- as.Date(dgfr$Date, format = "%m/%d/%y")
              dgfr$yday <- yday(dgfr$date)
              dgfr$temp_C <- (dgfr$Temp - 32) * 5/9
              dgfr$year <- year(dgfr$date)
              dgfr$yr_day_hr <- paste(dgfr$year, dgfr$yday, dgfr$Hour, sep = "_")
              dgfr$hour <- dgfr$Hour
        
        # combine the ithaca airport and game farm road data
              w_temp1 <- dweath[, c("year", "yday", "hour", "yr_day_hr", "temp_C")]
              w_temp2 <- dgfr[, c("year", "yday", "hour", "yr_day_hr", "temp_C")]
              w_hour <- rbind(w_temp1, w_temp2)
        
        # cleaned up hourly data for joining with behavior (only field season dates needed)    
            dw2 <- w_hour[, c("yr_day_hr", "temp_C")]
            colnames(dw2)[2] <- "ambient_C"
            dw2 <- dw2[dw2$yday > 110 & dw2$yday < 200, ]
            
        # create object with morning temperature for each day averaged 6-10am
            morning_temp <- w_hour %>%
              filter(hour > 5, hour < 11) %>%
              dplyr::group_by(yday, year) %>%
              dplyr::summarise(m_temp = mean(temp_C, na.rm = TRUE)) %>%
              filter(m_temp > -50, m_temp < 100)
            morning_temp$yr_doy <- paste(morning_temp$year, morning_temp$yday, sep = "_")
        
        # daily data for joining with long term database
            # cleaned up to just save max/min/avg but could be modified
              w_daily <- w_hour %>%
                filter(hour > 5, hour < 21) %>%
                group_by(year, yday) %>%
                summarise(mean_C = mean(temp_C, na.rm = TRUE),
                          high_C = max(temp_C, na.rm = TRUE),
                          low_C = min(temp_C, na.rm = TRUE))
              w_daily$yr_day <- paste(w_daily$year, w_daily$yday, sep = "_")
              
        # more processing of temperature is done in sections below, but this establishes
              # the main data objects needed
              
# Join weather with incubation and provisioning data ----
              
          # about is every bout in one row joined with temperature for that hour
              all_bouts_ex$hour <- round((all_bouts_ex$time_start_full - all_bouts_ex$JDate * 24 * 60 * 60) 
                                         / 60 / 60, 0)
              all_bouts_ex$yr_day_hr <- paste(all_bouts_ex$year, all_bouts_ex$JDate, all_bouts_ex$hour, sep = "_")
              about <- all_bouts_ex %>% 
                filter(exclude_all == 0, is_night == 0, (duration / 60) < 120, duration > 0)
              about <- plyr::join(about, dw2, by = "yr_day_hr", type = "left")
              about$uby <- paste(about$unit, about$nest, about$year, sep = "_")
              
          # each hour feeding with temperature for that hour
              feeds$yr_day_hr <- paste(feeds$year, feeds$yday, feeds$Hour, sep = "_")
              feeds <- plyr::join(feeds, dw2, "yr_day_hr")
              
# Clean up long term capture records ----
  # the long term database includes all kinds of data starting in 1986 and has some errors.
    # this is filtering down to a clean set of records for adults and nestlings
    # by removing impossible values and some treatments etc
              
    # first join capture and nest records
        ith_cap_nest <- plyr::join(ith_capture, nests, "nest_key", "left", "first")   
    
    # ADULTS
        # drop some columns and impossible values
          ith_adults <- ith_cap_nest[, c("nest_key", "encounter_key", "band", "rfid1", "rfid2", "encounter_doy",
                                         "exp_year", "adult_or_nestling", "sex", "individual_experiment", "individual_treatment",
                                         "site", "nest", "capture_num", "age", "mass", "headbill", "flatwing", "bleed1_latency_sec",
                                         "stress_series_type", "bleed1_cort", "bleed2_cort", "bleed3_cort", "release_time",
                                         "nestling_fate", "attempt_num", "clutch_init_doy",
                                         "hatch_doy", "clutch_size", "brood_size_hatching", "nest_fate_doy", "nest_fate", "simple_fate",
                                         "clutch_comp_doy", "num_fledged", "nest_experiment", "nest_treatment")] %>%
            filter(adult_or_nestling == "Adult",
                   mass > 10, mass < 30,
                   encounter_doy > 110, encounter_doy < 210)
        # add column for stage in the breeding cycle
          ith_adults$stage <- ifelse(
            !is.na(ith_adults$hatch_doy),
            ith_adults$encounter_doy - ith_adults$hatch_doy,
            ith_adults$encounter_doy - (ith_adults$clutch_comp_doy + 14)
          )
        # delete impossible values
          ith_adults <- ith_adults %>% filter(stage > -13, stage < 17)
        # remove adults of unknown sex (usually not captured at boxes)
          ith_adults <- subset(ith_adults, ith_adults$sex == "Female" | ith_adults$sex == "Male")
          
      # filtering out manipulative treatments after first capture
          # Step 1: Prepare exclude_trt for joining. This is a list of treatments to exclude
              exclude_trt <- exclude_trt %>% distinct(treatments, .keep_all = TRUE)
              exclude_trt <- exclude_trt %>%
                mutate(treatment = treatments) # Rename to make joining explicit
 
          # Step 2: Join ith_adults with exclude_trt
               ith_adults <- ith_adults %>%
                # Join on individual_treatment and nest_treatment
                left_join(exclude_trt, by = c("individual_treatment" = "treatment")) %>%
                left_join(exclude_trt, by = c("nest_treatment" = "treatment"), suffix = c("_ind", "_nest")) %>%
            
          # Step 3: Apply filtering logic
                filter(
                  # Always include rows where capture_num == 1 because these are pre-treatment
                  capture_num == 1 |
                    # Exclude rows based on exclude_female or exclude_male for capture_num > 1
                    (
                      capture_num > 1 &
                        (
                          (sex == "Female" & exclude_female_ind == "no" & exclude_female_nest == "no") |
                            (sex == "Male" & exclude_male_ind == "no" & exclude_male_nest == "no")
                        )
                    )
                ) %>%
                # Remove unnecessary columns after filtering
                select(-starts_with("exclude_"))
               
      # NESTLINGS
          # filter down to needed columns
               ith_nestlings <- ith_cap_nest[, c("nest_key", "encounter_key", "band", "rfid1", "rfid2", "encounter_doy",
                                              "exp_year", "adult_or_nestling", "sex", "individual_experiment", "individual_treatment",
                                              "site", "nest", "capture_num", "age", "mass", "headbill", "flatwing", "bleed1_latency_sec",
                                              "stress_series_type", "bleed1_cort", "bleed2_cort", "bleed3_cort", "release_time",
                                              "nestling_fate", "attempt_num", "clutch_init_doy",
                                              "hatch_doy", "clutch_size", "brood_size_hatching", "nest_fate_doy", "nest_fate", "simple_fate",
                                              "clutch_comp_doy", "num_fledged", "nest_experiment", "nest_treatment")] %>%
                 filter(adult_or_nestling == "Nestling",
                        mass > 0, mass < 30,
                        encounter_doy > 110, encounter_doy < 210)
               ith_nestlings$stage <- as.numeric(ith_nestlings$age)
               ith_nestlings <- ith_nestlings %>% filter(stage > -1, stage < 18) 
               ith_nestlings$uby <- paste(ith_nestlings$site, ith_nestlings$nest, ith_nestlings$exp_year, sep = "_")
               
          # filtering out manipulative treatments

              # Step 1: Join ith_nestlings with exclude_trt
                    ith_nestlings <- ith_nestlings %>%
              # Join on individual_treatment and nest_treatment
                 left_join(exclude_trt, by = c("individual_treatment" = "treatment")) %>%
                 left_join(exclude_trt, by = c("nest_treatment" = "treatment"), suffix = c("_ind", "_nest")) %>%
                 
              # Step 2: Apply filtering logic
                 filter(exclude_nestling_ind == "no", exclude_nestling_nest == "no")  %>%
              # Remove unnecessary columns after filtering
                 select(-starts_with("exclude_"))
              
              
# Effects of temperature on adult and nestling tree swallow mass/morph ----
    # will be summarized to FIGURE 1 and some supplemental figures
                    
        # add temperature data at different intervals to adult data
                # create a frame with each combination of year and yday captures occured on
                   yr_doy_uni1 <- ith_adults %>%
                     dplyr::group_by(exp_year, encounter_doy) %>%
                     summarise(count = n()) %>%
                     as.data.frame()
                   
                # filter hourly temperature to relevant period and daytime (6am-8pm)
                   wh_day <- w_hour[is.na(w_hour$temp_C) == FALSE & w_hour$hour > 5 & w_hour$hour < 21 &
                                      w_hour$yday > 100 & w_hour$yday < 250, ]
                   
                # for each year/day combination, calculate various lagged temperature averages (note not all used) 
                   for(i in 1:nrow(yr_doy_uni1)){
                     sub1 <- subset(wh_day, wh_day$year == yr_doy_uni1$exp_year[i] &
                                      wh_day$yday > yr_doy_uni1$encounter_doy[i] - 20 &
                                      wh_day$yday < yr_doy_uni1$encounter_doy[i] + 1)
                     sub2 <- sub1[sub1$yday == yr_doy_uni1$encounter_doy[i] & sub1$hour < 11, ]
                     yr_doy_uni1$cap_day_C[i] <- mean(sub2$temp_C)
                     yr_doy_uni1$prior1_C[i] <- mean(c(sub1$temp_C[sub1$yday > (yr_doy_uni1$encounter_doy[i] - 2) &
                                                                    sub1$yday < yr_doy_uni1$encounter_doy[i]], sub2$temp_C))
                     yr_doy_uni1$prior2_C[i] <- mean(c(sub1$temp_C[sub1$yday > (yr_doy_uni1$encounter_doy[i] - 3) &
                                                                    sub1$yday < yr_doy_uni1$encounter_doy[i]], sub2$temp_C))
                     yr_doy_uni1$prior3_C[i] <- mean(c(sub1$temp_C[sub1$yday > (yr_doy_uni1$encounter_doy[i] - 4) &
                                                                    sub1$yday < yr_doy_uni1$encounter_doy[i]], sub2$temp_C))
                     yr_doy_uni1$prior4_C[i] <- mean(c(sub1$temp_C[sub1$yday > (yr_doy_uni1$encounter_doy[i] - 5) &
                                                                    sub1$yday < yr_doy_uni1$encounter_doy[i]], sub2$temp_C))
                     yr_doy_uni1$prior5_C[i] <- mean(c(sub1$temp_C[sub1$yday > (yr_doy_uni1$encounter_doy[i] - 6) &
                                                                    sub1$yday < yr_doy_uni1$encounter_doy[i]], sub2$temp_C))
                     yr_doy_uni1$prior6_C[i] <- mean(c(sub1$temp_C[sub1$yday > (yr_doy_uni1$encounter_doy[i] - 7) &
                                                                    sub1$yday < yr_doy_uni1$encounter_doy[i]], sub2$temp_C))
                     yr_doy_uni1$prior7_C[i] <- mean(c(sub1$temp_C[sub1$yday > (yr_doy_uni1$encounter_doy[i] - 8) &
                                                                    sub1$yday < yr_doy_uni1$encounter_doy[i]], sub2$temp_C))
                   }
                   
                # join the lagged temperature captures to the adult captures
                   yr_doy_uni1$year_capday <- paste(yr_doy_uni1$exp_year, yr_doy_uni1$encounter_doy, sep = "_")
                   ith_adults$year_capday <- paste(ith_adults$exp_year, ith_adults$encounter_doy, sep = "_")
                   ith_adults <- plyr::join(ith_adults, yr_doy_uni1, "year_capday", "left", "first")
              
                # joining may result in some duplicated columns, filter those out        
                   ith_adults <- ith_adults[, !duplicated(names(ith_adults))]
               
      ## ADULT mass loss panel Figure 1A
            # summarize to each day relative to hatching separately for males/females       
                ith_ad2 <- ith_adults %>%
                     dplyr::group_by(sex, stage) %>%
                     dplyr::summarise(mass_mu = mean(mass, na.rm = TRUE),
                                      mass_sd = sd(mass, na.rm = TRUE),
                                      n = n())
                   colnames(ith_ad2)[2] <- "stage"
                   ith_ad2$sem <- ith_ad2$mass_sd / sqrt(ith_ad2$n)
               
            # make the plot for fig 1A
                ad_plot <- ggplot(ith_adults, mapping = aes(x = stage, y = mass, color = sex, fill = sex)) +
                  geom_vline(xintercept = 0, linetype = "dashed") +
                  geom_smooth() +
                  scale_color_manual(values = c(fem_color, mal_color), labels = c("Females", "Males"), name = "") +
                  scale_fill_manual(values = c(fem_color, mal_color), labels = c("Females", "Males"), name = "") +
                  theme_classic() +
                  scale_x_continuous(breaks = seq(-20, 20, 5)) +
                  theme(panel.grid = element_blank(), axis.text = element_text(size = 12), axis.title = element_text(size = 14)) +
                  theme(legend.position = c(0.85,0.85)) +
                  ylab("Mass (grams)") + xlab("Day relative to hatching") +
                  geom_segment(data = ith_ad2, mapping = aes(x = stage, y = mass_mu - sem, xend = stage, yend = mass_mu + sem)) +
                  geom_point(data = ith_ad2, mapping = aes(x = stage, y = mass_mu), shape = 21, color = "black", alpha = 0.7) +
                  coord_cartesian(xlim = c(-12, 16))
             
            # save the figure for use in illustrator   
                ggsave(here::here("output_plots", "fig1_a.svg"), ad_plot, device = "svg", width = 4, height = 2.5, units = "in")
              
      ## ADULT panel C & D. Mass vs. temperature for males and females   
            # fit models to plot partial effects of stage and temperature for males and females
            # one model for each temperature window
              # females day of capture, prior 2 + morning, prior 4 + morning, average temperature
                # mf_m1 <- gamm4(mass ~ s(stage) + s(cap_day_C), random = ~(1|band), data = ith_adults[ith_adults$sex == "Female", ])
                # mf_m3 <- gamm4(mass ~ s(stage) + s(prior2_C), random = ~(1|band), data = ith_adults[ith_adults$sex == "Female", ])
                # mf_m5 <- gamm4(mass ~ s(stage) + s(prior4_C), random = ~(1|band), data = ith_adults[ith_adults$sex == "Female", ])
                
                    # saveRDS(mf_m1, here::here("saved_models", "mf_m1.rds"))
                    # saveRDS(mf_m3, here::here("saved_models", "mf_m3.rds"))
                    # saveRDS(mf_m5, here::here("saved_models", "mf_m5.rds"))
                    mf_m1 <- readRDS(here::here("saved_models", "mf_m1.rds"))
                    mf_m3 <- readRDS(here::here("saved_models", "mf_m3.rds"))
                    mf_m5 <- readRDS(here::here("saved_models", "mf_m5.rds"))
                
              # males day of capture, prior 2 + morning, prior 4 + morning, average temperature
                # mm_m1 <- gamm4(mass ~ s(stage) + s(cap_day_C), random = ~(1|band), data = ith_adults[ith_adults$sex == "Male", ])
                # mm_m3 <- gamm4(mass ~ s(stage) + s(prior2_C), random = ~(1|band), data = ith_adults[ith_adults$sex == "Male", ])
                # mm_m5 <- gamm4(mass ~ s(stage) + s(prior4_C), random = ~(1|band), data = ith_adults[ith_adults$sex == "Male", ])
                
                    # saveRDS(mm_m1, here::here("saved_models", "mm_m1.rds"))
                    # saveRDS(mm_m3, here::here("saved_models", "mm_m3.rds"))
                    # saveRDS(mm_m5, here::here("saved_models", "mm_m5.rds"))
                    mm_m1 <- readRDS(here::here("saved_models", "mm_m1.rds"))
                    mm_m3 <- readRDS(here::here("saved_models", "mm_m3.rds"))
                    mm_m5 <- readRDS(here::here("saved_models", "mm_m5.rds"))
                

            # make plots from these adult models of partial effect of temperature on mass
                
                # FEMALE PLOT (panel C)
                    # extract partial effects from each model
                          d1_fem <- draw(mf_m1, term = "cap_day_C")$data %>%
                            select(cap_day_C = cap_day_C, fit = .estimate, se = .se, cil = .lower_ci, ciu = .upper_ci)
                          d3_fem <- draw(mf_m3, term = "prior2_C")$data %>%
                            select(prior2_C = prior2_C, fit = .estimate, se = .se, cil = .lower_ci, ciu = .upper_ci)
                          d5_fem <- draw(mf_m5, term = "prior4_C")$data %>%
                            select(prior4_C = prior4_C, fit = .estimate, se = .se, cil = .lower_ci, ciu = .upper_ci)
                     
                    # make the plot 
                          fem_part_plot <- ggplot() +
                            geom_abline(slope = 0, intercept = 0, linetype = "61", color = "gray75", size = 0.3) +
                            geom_line(data = d1_fem, aes(x = cap_day_C, y = fit), color = fem_color, linetype = "solid", size = 0.3) +
                            geom_ribbon(data = d1_fem, aes(x = cap_day_C, ymin = cil, ymax = ciu), fill = fem_color, alpha = 0.25) +
                            geom_line(data = d3_fem, aes(x = prior2_C, y = fit), color = fem_color, linetype = "31", size = 0.3) +
                            geom_ribbon(data = d3_fem, aes(x = prior2_C, ymin = cil, ymax = ciu), fill = fem_color, alpha = 0.25) +
                            geom_line(data = d5_fem, aes(x = prior4_C, y = fit), color = fem_color, linetype = "11", size = 0.3) +
                            geom_ribbon(data = d5_fem, mapping = aes(x = prior4_C, ymin = cil, ymax = ciu), fill = fem_color, alpha = 0.25) +
                            labs(x = "Temperature C", y = "Predicted Partial Mass (grams)") +
                            theme_classic() +
                            scale_x_continuous(breaks = seq(0, 30, 10)) +
                          xlim(c(0, 30)) +
                          theme(panel.grid = element_blank(), axis.text = element_text(size = 12), axis.title = element_text(size = 14))
                       
                    # save the plot 
                        ggsave(here::here("output_plots", "fig1_c.svg"), fem_part_plot, device = "svg", width = 2, height = 2, units = "in")
                    
                # MALE PLOT (panel D)
                    # extract partial effects
                        d1_mal <- draw(mm_m1, term = "cap_day_C")$data %>%
                          select(cap_day_C = cap_day_C, fit = .estimate, se = .se, cil = .lower_ci, ciu = .upper_ci)
                        d3_mal <- draw(mm_m3, term = "prior2_C")$data %>%
                          select(prior2_C = prior2_C, fit = .estimate, se = .se, cil = .lower_ci, ciu = .upper_ci)
                        d5_mal <- draw(mm_m5, term = "prior4_C")$data %>%
                          select(prior4_C = prior4_C, fit = .estimate, se = .se, cil = .lower_ci, ciu = .upper_ci)
                      
                    # make the plot  
                        mal_part_plot <- ggplot() +
                          geom_abline(slope = 0, intercept = 0, linetype = "61", color = "gray75", size = .3) +
                          geom_line(data = d1_mal, aes(x = cap_day_C, y = fit), color = mal_color, linetype = "solid", size = .3) +
                          geom_ribbon(data = d1_mal, aes(x = cap_day_C, ymin = cil, ymax = ciu), fill = mal_color, alpha = 0.25) +
                          geom_line(data = d3_mal, aes(x = prior2_C, y = fit), color = mal_color, linetype = "31", size = .3) +
                          geom_ribbon(data = d3_mal, aes(x = prior2_C, ymin = cil, ymax = ciu), fill = mal_color, alpha = 0.25) +
                          geom_line(data = d5_mal, aes(x = prior4_C, y = fit), color = mal_color, linetype = "11", size = .3) +
                          geom_ribbon(data = d5_mal, mapping = aes(x = prior4_C, ymin = cil, ymax = ciu), fill = mal_color, alpha = 0.25) +
                          labs(x = "Temperature C", y = "Predicted Partial Mass (grams)") +
                          theme_classic() +
                          scale_x_continuous(breaks = seq(0, 30, 10)) +
                          scale_y_continuous(breaks = seq(-2, 2, 1)) +
                          xlim(c(0, 30)) +
                          theme(panel.grid = element_blank(), axis.text = element_text(size = 12), axis.title = element_text(size = 14))
                       
                     # save the plot 
                        ggsave(here::here("output_plots", "fig1_d.svg"), mal_part_plot, device = "svg", width = 2, height = 2, units = "in")
                  
        # NESTLING GROWTH PANEL
              # avg adult size to be used as a reference
                    avg_ad_mass <- mean(na.omit(ith_adults$mass[ith_adults$mass < 30 & ith_adults$mass > 10]))
                    avg_ad_hb <- mean(na.omit(ith_adults$headbill[ith_adults$headbill < 33 & ith_adults$headbill > 23]))
                    avg_ad_fw <- mean(na.omit(ith_adults$flatwing[ith_adults$flatwing < 130 & ith_adults$flatwing > 90]))
                    
              # add relative morph columns. these are nestling measures as a % of average adult size
                    ith_nestlings$rel_mass <- ith_nestlings$mass / avg_ad_mass
                    ith_nestlings$rel_hb <- ith_nestlings$headbill / avg_ad_hb
                    ith_nestlings$rel_fw <- ith_nestlings$flatwing / avg_ad_fw
                    
              # Summarize nestling measures by each age in days
                    ith_ne2 <- ith_nestlings %>%
                      dplyr::group_by(stage) %>%
                      dplyr::summarise(mass_rel_mu = mean(rel_mass, na.rm = TRUE),
                                       mass_rel_sd = sd(rel_mass, na.rm = TRUE),
                                       hb_rel_mu = mean(rel_hb, na.rm = TRUE),
                                       hb_rel_sd = mean(rel_hb, na.rm = TRUE),
                                       fw_rel_mu = mean(rel_fw, na.rm = TRUE),
                                       fw_rel_sd = mean(rel_fw, na.rm = TRUE),
                                       n = n())
                    ith_ne2$sem_mass <- ith_ne2$mass_rel_sd / sqrt(ith_ne2$n)
                    ith_ne2$sem_hb <- ith_ne2$hb_rel_sd / sqrt(ith_ne2$n)
                    ith_ne2$sem_fw <- ith_ne2$fw_rel_sd / sqrt(ith_ne2$n)
                  
              # Convert the nestlng measures to a long data frame to make plotting easier
                    ith_ne2long <- data.frame(stage = rep(ith_ne2$stage, 3),
                                              cat = c(rep("Mass", nrow(ith_ne2)), rep("Head", nrow(ith_ne2)), rep("Wing", nrow(ith_ne2))),
                                              rel_meas = c(ith_ne2$mass_rel_mu, ith_ne2$hb_rel_mu, ith_ne2$fw_rel_mu),
                                              rel_sd = c(ith_ne2$mass_rel_sd, ith_ne2$hb_rel_sd, ith_ne2$fw_rel_sd),
                                              sem = c(ith_ne2$sem_mass, ith_ne2$sem_hb, ith_ne2$sem_fw))
                  
              # remove any duplicated columns introduced by joining  
                    ith_nestlings <- ith_nestlings[, !duplicated(names(ith_nestlings))]
                    
              # convert captures to a long format to make plotting easier
                    ith_long <- data.frame(stage = rep(ith_nestlings$stage, 3),
                                           cat = c(rep("Mass", nrow(ith_nestlings)), rep("Head", nrow(ith_nestlings)), rep("Wing", nrow(ith_nestlings))),
                                           measure = c(ith_nestlings$rel_mass*100, ith_nestlings$rel_hb*100, ith_nestlings$rel_fw*100))
                
              # construct the plot    
                    ne_plot <- ggplot(ith_long, mapping = aes(x = stage, y = measure, fill = cat, color = cat)) +
                      geom_smooth() +
                      geom_point(data = ith_ne2long, mapping = aes(x = stage, y = rel_meas*100), shape = 21, color = "black", size = 1.7) +
                      geom_segment(data = ith_ne2long, mapping = aes(x = stage, xend = stage, y = rel_meas*100 + sem*100, yend = rel_meas*100 - sem*100)) +
                      scale_y_continuous(breaks = seq(-25, 100, 25)) +
                      coord_cartesian(ylim = c(0, 108)) +
                      geom_vline(xintercept = 0, linetype = "dashed") +
                      theme_classic() +
                      theme(panel.grid = element_blank(), axis.text = element_text(size = 12), axis.title = element_text(size = 14)) +
                      labs(x = "Day relative to hatching", y = "Percent of Adult Size") +
                      scale_color_manual(values = c("Mass" = mg_color, "Wing" = fw_color, "Head" = hb_color)) +
                      scale_fill_manual(values = c("Mass" = mg_color, "Wing" = fw_color, "Head" = hb_color)) +
                      theme(legend.position = c(0.8, 0.2), legend.title = element_blank())
                
                # save the plot  
                    ggsave(here::here("output_plots", "fig1_b.svg"), ne_plot, device = "svg", width = 4, height = 2.5, units = "in")  
                   
          # NESTLING TEMPERATURE PANELS (E-G) 
                  # add temperature at different intervals to nestling data
                    # get each unique year/day combo
                        yr_doy_uni <- ith_nestlings %>%
                          dplyr::group_by(exp_year, encounter_doy) %>%
                          summarise(count = n()) %>%
                          as.data.frame()
                     
                    # filter down hourly temperature data to daytime (6am-8pm)   
                        wh_day <- w_hour[is.na(w_hour$temp_C) == FALSE & w_hour$hour > 5 & w_hour$hour < 21 &
                                           w_hour$yday > 110 & w_hour$yday < 250, ]
                      
                    # for each year/day combo, calculate lagged temperature. Note slight difference with adults
                        # for day of capture where adults are 6-10am (because captures usually in morning),
                        # but nestlings are 6am-12pm because measures usually around midday
                            for(i in 1:nrow(yr_doy_uni)){
                              sub1 <- subset(wh_day, wh_day$year == yr_doy_uni$exp_year[i] &
                                              wh_day$yday > yr_doy_uni$encounter_doy[i] - 20 &
                                               wh_day$yday < yr_doy_uni$encounter_doy[i] + 1)
                              sub2 <- sub1[sub1$yday == yr_doy_uni$encounter_doy[i] & sub1$hour < 13, ]
                              yr_doy_uni$cap_day_C[i] <- mean(sub2$temp_C)
                              yr_doy_uni$prior1_C[i] <- mean(c(sub1$temp_C[sub1$yday > (yr_doy_uni$encounter_doy[i] - 2) &
                                                                           sub1$yday < yr_doy_uni$encounter_doy[i]], sub2$temp_C))
                              yr_doy_uni$prior2_C[i] <- mean(c(sub1$temp_C[sub1$yday > (yr_doy_uni$encounter_doy[i] - 3) &
                                                                             sub1$yday < yr_doy_uni$encounter_doy[i]], sub2$temp_C))
                              yr_doy_uni$prior3_C[i] <- mean(c(sub1$temp_C[sub1$yday > (yr_doy_uni$encounter_doy[i] - 4) &
                                                                             sub1$yday < yr_doy_uni$encounter_doy[i]], sub2$temp_C))
                              yr_doy_uni$prior4_C[i] <- mean(c(sub1$temp_C[sub1$yday > (yr_doy_uni$encounter_doy[i] - 5) &
                                                                             sub1$yday < yr_doy_uni$encounter_doy[i]], sub2$temp_C))
                              yr_doy_uni$prior5_C[i] <- mean(c(sub1$temp_C[sub1$yday > (yr_doy_uni$encounter_doy[i] - 6) &
                                                                             sub1$yday < yr_doy_uni$encounter_doy[i]], sub2$temp_C))
                              yr_doy_uni$prior6_C[i] <- mean(c(sub1$temp_C[sub1$yday > (yr_doy_uni$encounter_doy[i] - 7) &
                                                                             sub1$yday < yr_doy_uni$encounter_doy[i]], sub2$temp_C))
                              yr_doy_uni$prior7_C[i] <- mean(c(sub1$temp_C[sub1$yday > (yr_doy_uni$encounter_doy[i] - 8) &
                                                                             sub1$yday < yr_doy_uni$encounter_doy[i]], sub2$temp_C))
                            }
                        
                        # join nestlings to temperature data for models
                            yr_doy_uni$year_capday <- paste(yr_doy_uni$exp_year, yr_doy_uni$encounter_doy, sep = "_")
                            ith_nestlings$year_capday <- paste(ith_nestlings$exp_year, ith_nestlings$encounter_doy, sep = "_")
                            ith_nestlings <- plyr::join(ith_nestlings, yr_doy_uni, "year_capday", "left", "first")
                    
                  # NESTLING temperature effect models. three per measure for 1/3/5 day temperature
                        ith_nestlings$uby <- as.factor(ith_nestlings$uby)
                      # mass
                        # nma_mo <- gamm4(mass ~ s(stage) + t2(stage, cap_day_C, bs = c("tp", "cr")), random = ~(1|uby), data = ith_nestlings)
                        # nma_p2 <- gamm4(mass ~ s(stage) + t2(stage, prior2_C, bs = c("tp", "cr")), random = ~(1|uby), data = ith_nestlings)
                        # nma_p4 <- gamm4(mass ~ s(stage) + t2(stage, prior4_C, bs = c("tp", "cr")), random = ~(1|uby), data = ith_nestlings)
                        
                            # saveRDS(nma_mo, here::here("saved_models", "nma_mo.rds"))
                            # saveRDS(nma_p2, here::here("saved_models", "nma_p2.rds"))
                            # saveRDS(nma_p4, here::here("saved_models", "nma_p4.rds"))
                            nma_mo <- readRDS(here::here("saved_models", "nma_mo.rds"))
                            nma_p2 <- readRDS(here::here("saved_models", "nma_p2.rds"))
                            nma_p4 <- readRDS(here::here("saved_models", "nma_p4.rds"))
                      # wing
                        # nfw_mo <- gamm4(flatwing ~ s(stage) + t2(stage, cap_day_C, bs = c("tp", "cr")), random = ~(1|uby), data = ith_nestlings)
                        # nfw_p2 <- gamm4(flatwing ~ s(stage) + t2(stage, prior2_C, bs = c("tp", "cr")), random = ~(1|uby), data = ith_nestlings)
                        # nfw_p4 <- gamm4(flatwing ~ s(stage) + t2(stage, prior4_C, bs = c("tp", "cr")), random = ~(1|uby), data = ith_nestlings)
                        
                            # saveRDS(nfw_mo, here::here("saved_models", "nfw_mo.rds"))
                            # saveRDS(nfw_p2, here::here("saved_models", "nfw_p2.rds"))
                            # saveRDS(nfw_p4, here::here("saved_models", "nfw_p4.rds"))
                            nfw_mo <- readRDS(here::here("saved_models", "nfw_mo.rds"))
                            nfw_p2 <- readRDS(here::here("saved_models", "nfw_p2.rds"))
                            nfw_p4 <- readRDS(here::here("saved_models", "nfw_p4.rds"))
                      # hbill
                        # nhb_mo <- gamm4(headbill ~ s(stage) + t2(stage, cap_day_C, bs = c("tp", "cr")), random = ~(1|uby), data = ith_nestlings)
                        # nhb_p2 <- gamm4(headbill ~ s(stage) + t2(stage, prior2_C, bs = c("tp", "cr")), random = ~(1|uby), data = ith_nestlings)
                        # nhb_p4 <- gamm4(headbill ~ s(stage) + t2(stage, prior4_C, bs = c("tp", "cr")), random = ~(1|uby), data = ith_nestlings)
                        
                            # saveRDS(nhb_mo, here::here("saved_models", "nhb_mo.rds"))
                            # saveRDS(nhb_p2, here::here("saved_models", "nhb_p2.rds"))
                            # saveRDS(nhb_p4, here::here("saved_models", "nhb_p4.rds"))
                            nhb_mo <- readRDS(here::here("saved_models", "nhb_mo.rds"))
                            nhb_p2 <- readRDS(here::here("saved_models", "nhb_p2.rds"))
                            nhb_p4 <- readRDS(here::here("saved_models", "nhb_p4.rds"))
                    
                    
                        # NESTLING MASS by temperature (panel E)
                          # get partial effects
                              pe_nma_mo <- smooth_estimates(nma_mo$gam, smooth = "t2(stage,cap_day_C)", data = data.frame(stage = 8, cap_day_C = seq(8, 29, 0.2)))
                                colnames(pe_nma_mo) <- c("smooth", "type", "by", "estimate", "se", "stage", "cap_day_C")
                                pe_nma_mo$lower_ci = pe_nma_mo$estimate - 1.96 * pe_nma_mo$se
                                pe_nma_mo$upper_ci = pe_nma_mo$estimate + 1.96 * pe_nma_mo$se
                                
                              pe_nma_p2 <- smooth_estimates(nma_p2$gam, smooth = "t2(stage,prior2_C)", data = data.frame(stage = 8, prior2_C = seq(11, 29, 0.2)))
                                colnames(pe_nma_p2) <- c("smooth", "type", "by", "estimate", "se", "stage", "prior2_C")
                                pe_nma_p2$lower_ci = pe_nma_p2$estimate - 1.96 * pe_nma_p2$se
                                pe_nma_p2$upper_ci = pe_nma_p2$estimate + 1.96 * pe_nma_p2$se
                                
                              pe_nma_p4 <- smooth_estimates(nma_p4$gam, smooth = "t2(stage,prior4_C)", data = data.frame(stage = 8, prior4_C = seq(12, 29, 0.2)))
                                colnames(pe_nma_p4) <- c("smooth", "type", "by", "estimate", "se", "stage", "prior4_C")
                                pe_nma_p4$lower_ci = pe_nma_p4$estimate - 1.96 * pe_nma_p4$se
                                pe_nma_p4$upper_ci = pe_nma_p4$estimate + 1.96 * pe_nma_p4$se
                                
                          # make the plot
                              ne_mass_plot <- ggplot(pe_nma_mo, mapping = aes(x = cap_day_C, y = estimate)) +
                                geom_abline(slope = 0, intercept = 0, linetype = "61", color = "gray75", size = 0.3) +
                                geom_line(color = mg_color, size = 0.3) +
                                geom_ribbon(aes(ymin = lower_ci, ymax = upper_ci), alpha = 0.25, fill = mg_color) +
                                geom_line(data = pe_nma_p2, linetype = "31", mapping = aes(x = prior2_C), color = mg_color, size = 0.3) +
                                geom_ribbon(data = pe_nma_p2, mapping = aes(x = prior2_C, ymin = lower_ci, ymax = upper_ci), fill = mg_color, alpha = 0.25) +
                                geom_line(data = pe_nma_p4, linetype = "11", mapping = aes(x = prior4_C), color = mg_color, size = 0.3) +
                                geom_ribbon(data = pe_nma_p4, mapping = aes(x = prior4_C, ymin = lower_ci, ymax = upper_ci), fill = mg_color, alpha = 0.25) +
                                theme_classic() +
                                theme(panel.grid = element_blank(), axis.text = element_text(size = 12), axis.title = element_text(size = 14)) +
                                labs(x = "Temperature (C)", y = "Partial Effect") +
                                scale_x_continuous(breaks = seq(-10,40, 10)) +
                                scale_y_continuous(breaks = seq(-6,6,3))
                              
                              ggsave(here::here("output_plots", "fig1_e.svg"), ne_mass_plot, device = "svg", width = 1.6, height = 2, units = "in")
                              
                        # NESTLING FLATWING by temperature (panel F)
                            # get the partial effects
                              pe_nfw_mo <- smooth_estimates(nfw_mo$gam, smooth = "t2(stage,cap_day_C)", data = data.frame(stage = 8, cap_day_C = seq(10, 29, 0.2)))
                              colnames(pe_nfw_mo) <- c("smooth", "type", "by", "estimate", "se", "stage", "cap_day_C")
                              pe_nfw_mo$lower_ci = pe_nfw_mo$estimate - 1.96 * pe_nfw_mo$se
                              pe_nfw_mo$upper_ci = pe_nfw_mo$estimate + 1.96 * pe_nfw_mo$se
                              
                              pe_nfw_p2 <- smooth_estimates(nfw_p2$gam, smooth = "t2(stage,prior2_C)", data = data.frame(stage = 8, prior2_C = seq(11, 29, 0.2)))
                              colnames(pe_nfw_p2) <- c("smooth", "type", "by", "estimate", "se", "stage", "prior2_C")
                              pe_nfw_p2$lower_ci = pe_nfw_p2$estimate - 1.96 * pe_nfw_p2$se
                              pe_nfw_p2$upper_ci = pe_nfw_p2$estimate + 1.96 * pe_nfw_p2$se
                              
                              pe_nfw_p4 <- smooth_estimates(nfw_p4$gam, smooth = "t2(stage,prior4_C)", data = data.frame(stage = 8, prior4_C = seq(12, 29, 0.2)))
                              colnames(pe_nfw_p4) <- c("smooth", "type", "by", "estimate", "se", "stage", "prior4_C")
                              pe_nfw_p4$lower_ci = pe_nfw_p4$estimate - 1.96 * pe_nfw_p4$se
                              pe_nfw_p4$upper_ci = pe_nfw_p4$estimate + 1.96 * pe_nfw_p4$se
                            
                            # make the plot  
                              ne_wing_plot <- ggplot(pe_nfw_mo, mapping = aes(x = cap_day_C, y = estimate)) +
                                geom_abline(slope = 0, intercept = 0, linetype = "61", color = "gray75", size = 0.3) +
                                geom_line(color = fw_color, size = 0.3) +
                                geom_ribbon(aes(ymin = lower_ci, ymax = upper_ci), alpha = 0.25, fill = fw_color) +
                                geom_line(data = pe_nfw_p2, linetype = "31", mapping = aes(x = prior2_C), color = fw_color, size = 0.3) +
                                geom_ribbon(data = pe_nfw_p2, mapping = aes(x = prior2_C, ymin = lower_ci, ymax = upper_ci), fill = fw_color, alpha = 0.25) +
                                geom_line(data = pe_nfw_p4, linetype = "11", mapping = aes(x = prior4_C), color = fw_color, size = 0.3) +
                                geom_ribbon(data = pe_nfw_p4, mapping = aes(x = prior4_C, ymin = lower_ci, ymax = upper_ci), fill = fw_color, alpha = 0.25) +
                                theme_classic() +
                                theme(panel.grid = element_blank(), axis.text = element_text(size = 12), axis.title = element_text(size = 14)) +
                                labs(x = "Temperature (C)", y = "Partial Effect") +
                                scale_x_continuous(breaks = seq(-10,40, 10)) +
                                coord_cartesian(xlim=c(9,30))
                              
                              ggsave(here::here("output_plots", "fig1_f.svg"), ne_wing_plot, device = "svg", width = 1.6, height = 2, units = "in")
      
                        # NESTLING HEADBILL by temperature
                            # get the partial effects
                              pe_nhb_mo <- smooth_estimates(nhb_mo$gam, smooth = "t2(stage,cap_day_C)", data = data.frame(stage = 8, cap_day_C = seq(8, 29, 0.2)))
                              colnames(pe_nhb_mo) <- c("smooth", "type", "by", "estimate", "se", "stage", "cap_day_C")
                              pe_nhb_mo$lower_ci = pe_nhb_mo$estimate - 1.96 * pe_nhb_mo$se
                              pe_nhb_mo$upper_ci = pe_nhb_mo$estimate + 1.96 * pe_nhb_mo$se
                              
                              pe_nhb_p2 <- smooth_estimates(nhb_p2$gam, smooth = "t2(stage,prior2_C)", data = data.frame(stage = 8, prior2_C = seq(11, 29, 0.2)))
                              colnames(pe_nhb_p2) <- c("smooth", "type", "by", "estimate", "se", "stage", "prior2_C")
                              pe_nhb_p2$lower_ci = pe_nhb_p2$estimate - 1.96 * pe_nhb_p2$se
                              pe_nhb_p2$upper_ci = pe_nhb_p2$estimate + 1.96 * pe_nhb_p2$se
                              
                              pe_nhb_p4 <- smooth_estimates(nhb_p4$gam, smooth = "t2(stage,prior4_C)", data = data.frame(stage = 8, prior4_C = seq(12, 29, 0.2)))
                              colnames(pe_nhb_p4) <- c("smooth", "type", "by", "estimate", "se", "stage", "prior4_C")
                              pe_nhb_p4$lower_ci = pe_nhb_p4$estimate - 1.96 * pe_nhb_p4$se
                              pe_nhb_p4$upper_ci = pe_nhb_p4$estimate + 1.96 * pe_nhb_p4$se
                              
                            # make the plot
                              ne_head_plot <- ggplot(pe_nhb_mo, mapping = aes(x = cap_day_C, y = estimate)) +
                                geom_abline(slope = 0, intercept = 0, linetype = "61", color = "gray75", size = 0.3) +
                                geom_line(color = hb_color, size = 0.3) +
                                geom_ribbon(aes(ymin = lower_ci, ymax = upper_ci), alpha = 0.25, fill = hb_color) +
                                geom_line(data = pe_nhb_p2, linetype = "31", mapping = aes(x = prior2_C), color = hb_color, size = 0.3) +
                                geom_ribbon(data = pe_nhb_p2, mapping = aes(x = prior2_C, ymin = lower_ci, ymax = upper_ci), fill = hb_color, alpha = 0.25) +
                                geom_line(data = pe_nhb_p4, linetype = "11", mapping = aes(x = prior4_C), color = hb_color, size = 0.3) +
                                geom_ribbon(data = pe_nhb_p4, mapping = aes(x = prior4_C, ymin = lower_ci, ymax = upper_ci), fill = hb_color, alpha = 0.25) +
                                theme_classic() +
                                theme(panel.grid = element_blank(), axis.text = element_text(size = 12), axis.title = element_text(size = 14)) +
                                labs(x = "Temperature (C)", y = "Partial Effect") +
                                scale_x_continuous(breaks = seq(-10,40, 10))
                              
                              ggsave(here::here("output_plots", "fig1_g.svg"), ne_head_plot, device = "svg", width = 1.6, height = 2, units = "in")    
                              
                  # Extra plot used to make some elements to use in legends in illustrator
                        line_df <- data.frame(cat = rep(c("Female", "Male", "NMass", "NWing", "NHead"), 3),
                                              lntype = c(rep("solid", 5), rep("31", 5), rep("11", 5)),
                                              pcolor = rep(c(fem_color, mal_color, mg_color, fw_color, hb_color), 3),
                                              xstart = 1, xend = 2, ystart = seq(1, 15, 1), yend = seq(1, 15, 1))
                        lines <- ggplot(line_df, mapping = aes()) +
                          geom_segment(mapping = aes(x = xstart, xend = xend, y = ystart, yend = yend, linetype = lntype), size = 0.6, color = line_df$pcolor) +
                          theme_classic() + theme(panel.grid = element_blank()) +
                          coord_cartesian(xlim = c(0, 3), ylim = c(0, 16)) +
                          guides(fill = "none", color = "none", linetype = "none")
                        ggsave(here::here("output_plots", "lines.svg"), lines, device = "svg", width = 1.8, height = 2, units = "in")
                      
                    
# Effects of temperature on nestling fledging and recruitment ----
    # simple binomial to clutch (or band?) to fledge and fledge to recruit
      # can also cite our papers showing this, but can add model from this dataset
      # what range of days to use for temperature?
                        
      # first for survival from egg to fledge
          surv_nests <- nests
          surv_nests <- plyr::join(surv_nests, ith_capture[ith_capture$adult_or_nestling == "Adult" &
                                                             ith_capture$sex == "Female", ],
                                   "nest_key", "left", "first")      
          surv_nests <- surv_nests[, !duplicated(names(surv_nests))]
          
      # exclude by treatment to get rid of manipulations as above
          exclude_ind <- exclude_trt$treatments[exclude_trt$exclude_female == "yes"]
          surv_nests <- surv_nests[!(surv_nests$individual_treatment %in% exclude_ind), ]
          surv_nests <- surv_nests[!(surv_nests$nest_treatment %in% exclude_ind), ]
          surv_nests <- surv_nests %>%
            filter(location == "Ithaca", species == "TRES", is.na(hatch_doy) == FALSE, is.na(clutch_size) == FALSE,
                   clutch_size > -1, clutch_size < 9, attempt_num == 1) %>%
            select(nest_key, exp_year, clutch_init_doy, clutch_comp_doy, clutch_size, hatch_doy, brood_size_hatching,
                   nest_fate_doy, simple_fate, num_fledged, female_id, male_id) %>%
            filter(simple_fate == "failed" | simple_fate == "fledged") %>%
            filter(is.na(num_fledged) == FALSE, is.na(brood_size_hatching) == FALSE) %>%
            filter(num_fledged > -1, num_fledged < 9, brood_size_hatching > -1, brood_size_hatching < 9) %>%
            filter(num_fledged < brood_size_hatching)
          
      # join with temperature to look at x days prior to day 12. focusing on days 5-12 of growth
          surv_nests$year_capday <- paste(surv_nests$exp_year, surv_nests$hatch_doy + 12, sep = "_")
          surv_nests <- plyr::join(surv_nests, yr_doy_uni1[, 5:12], "year_capday", "left", "first")
          
          ydu1b <- yr_doy_uni1[, 5:12]
          colnames(ydu1b) <- paste("fled", colnames(ydu1b), sep = "_")
          surv_nests$fled_year_capday <- paste(surv_nests$exp_year, surv_nests$hatch_doy + 22, sep = "_")
          surv_nests <- plyr::join(surv_nests, ydu1b, "fled_year_capday", "left", "first")
          
      # add in number recruited from each nest
          all_nestlings <- ith_capture %>%
            filter(location == "Ithaca", adult_or_nestling == "Nestling") %>%
            select(nest_key, band) %>%
            filter(!duplicated(band))
          
          all_adults <- ith_capture %>%
            filter(location == "Ithaca", adult_or_nestling == "Adult") %>%
            select(band) %>%
            filter(!duplicated(band))
          all_adults$recruited <- "yes"
          
          all_nestlings <- plyr::join(all_nestlings, all_adults, "band", "left", "first")
          all_nestlings$recruited <- ifelse(is.na(all_nestlings$recruited) == TRUE, "no", "yes")
          
          for(i in 1:nrow(surv_nests)){
            sub <- subset(all_nestlings, all_nestlings$nest_key == surv_nests$nest_key[i] & all_nestlings$recruited == "yes")
            surv_nests$num_recruit[i] <- nrow(sub)
          }
          
      # make the columns needed for models
          surv_nests$not_fledged <- surv_nests$brood_size_hatching - surv_nests$num_fledged
          surv_nests$not_recruited <- surv_nests$num_fledged - surv_nests$num_recruit
          surv_nests[surv_nests$not_recruited < 0, "not_recruited"] <- 0
          surv_nests$prior7_Cs <- scale(surv_nests$prior7_C)
          
      # fit models
          
          mod_fled <- glmer(cbind(num_fledged, not_fledged) ~ scale(prior7_C) + (1|nest_key), data = surv_nests,
                            family = binomial)
          mod_fled2 <- glmmTMB(cbind(num_fledged, not_fledged) ~ scale(prior7_C) + (1|nest_key), data = surv_nests,
                               family = betabinomial)
          mod_recruit <-glmer(cbind(num_recruit, not_recruited) ~ scale(prior7_C) + (1|nest_key), data = surv_nests,
                              family = binomial)
          mod_recruit2 <- glmmTMB(cbind(num_recruit, not_recruited) ~ scale(fled_prior4_C) + scale(prior4_C) + (1|nest_key), data = surv_nests,
                                  family = binomial, ziformula = ~1)
      
      # generate data needed for plotting    
          min7 <- min(na.omit(surv_nests$prior7_C))
          max7 <- max(na.omit(surv_nests$prior7_C))
          mf2c <- coefficients(summary(mod_fled2))$cond
          fled_est <- data.frame(x = seq(-3, 3, 0.1),
                                 predicted = NA,
                                 conf.low = NA,
                                 conf.high = NA)
          fled_est$predicted <- inv_logit(mf2c[1,1] + fled_est$x * mf2c[2,1])
          fled_est$conf.low <- inv_logit(mf2c[1,1] + fled_est$x * mf2c[2,1] - 1.96 * mf2c[2,2])
          fled_est$conf.high <- inv_logit(mf2c[1,1] + fled_est$x * mf2c[2,1] + 1.96 * mf2c[2,2])
          
           mfled_pred <- as.data.frame(ggeffect(mod_fled, terms = "prior7_C [12:27 by=0.5]", transform = "response"))
        
      # create a plot for fledging (recruiting is not signficant)  
          mfled_plot <- ggplot(mfled_pred, aes(x = x, y= predicted*100)) +
            geom_line() +
            geom_ribbon(aes(ymin = conf.low*100, ymax = conf.high*100), alpha = 0.4) +
            scale_x_continuous(breaks = seq(0, 30, 5)) +
            coord_cartesian(xlim = c(13, 26)) +
            scale_y_continuous(breaks = seq(0, 100, 10), limits = c(5, 70)) +
            theme_classic() +
            theme(panel.grid = element_blank(), axis.text = element_text(size = 12),
                  axis.title = element_text(size = 14)) +
            labs(x = "Daytime Temperature (°C) \n Ages 5-12 Days", y = "Predicted Fledging Percentage")
          
          ggsave(here::here("output_plots", "overall_fled_plot.png"), mfled_plot, device = "png", width = 4, height = 4, units = "in", dpi = 300)
          
          
# Effects of breeding temperature on adult return rate ----
    # make a data frame starting from nests for survival    
        ad_surv <-  nests %>%
            filter(location == "Ithaca") %>%
            filter(site == "Unit_1" | site == "Unit_2" | site == "Unit_4" | site == "Unit_6") %>%
            filter(exp_year < 2024) %>%
            select(nest_key, exp_year, female_id, male_id, clutch_comp_doy, hatch_doy,
                   brood_size_hatching, num_fledged, nest_treatment)
          
        ad_surv <- plyr::join(ad_surv, ith_capture[ith_capture$sex == "Female" & 
                                                     ith_capture$adult_or_nestling == "Adult", 
                                                   c("nest_key", "individual_treatment")],
                              "nest_key", "left", "first")  
        
        ad_surv <- ad_surv[!(ad_surv$individual_treatment %in% exclude_trt$treatments[exclude_trt$exclude_female == "yes"]), ]
        ad_surv <- ad_surv[!(ad_surv$nest_treatment %in% exclude_trt$treatments[exclude_trt$exclude_female == "yes"]), ]
        
        ad_surv2 <- data.frame(band = c(ad_surv$female_id, ad_surv$male_id),
                               exp_year = rep(ad_surv$exp_year, 2),
                               clutch_comp_doy = rep(ad_surv$clutch_comp_doy, 2),
                               hatch_doy = rep(ad_surv$hatch_doy, 2),
                               sex = c(rep("Female", nrow(ad_surv)), rep("Male", nrow(ad_surv))))
        ad_surv2 <- subset(ad_surv2, ad_surv2$band > 1)
        
        ad_surv2$lived <- 0
        for(i in 1:nrow(ad_surv2)){
          sub <- subset(ith_capture, ith_capture$band == ad_surv2$band[i] &
                          ith_capture$exp_year > ad_surv2$exp_year[i])
          if(nrow(sub) > 0){
            ad_surv2$lived[i] <- 1
          }
        }
        
        
        # join with temperature to look at x days prior to day 12. focusing on days 5-12 of growth
          ad_surv2$year_capday <- paste(ad_surv2$exp_year, ad_surv2$hatch_doy + 12, sep = "_")
          ad_surv2 <- plyr::join(ad_surv2, yr_doy_uni1[, 5:12], "year_capday", "left", "first")
          
          ad_surv3 <- na.omit(ad_surv2)
          
          fem_surv <- ad_surv3[ad_surv3$sex == "Female",]
          mal_surv <- ad_surv3[ad_surv3$sex == "Male", ]
          
        # fit models
          m_f_surv <- glm(lived ~ scale(prior7_C), data = fem_surv, family = binomial)
          m_m_surv <- glm(lived ~ scale(prior7_C), data = mal_surv, family = binomial)
          m_b_surv <- glm(lived ~ scale(prior7_C)*sex, data = ad_surv3, family = binomial)
          
          # apparent male survival much lower, but probably because of incomplete surveys, focus on just females
          
          tab_model(m_f_surv, m_m_surv, m_b_surv)
          
          
          fsurv_pred <- as.data.frame(ggeffect(m_f_surv, terms = "prior7_C [12:27 by=0.5]", transform = "response"))
          
          fsurv_plot <- ggplot(fsurv_pred, aes(x = x, y= predicted*100)) +
            geom_line(color = fem_color) +
            geom_ribbon(aes(ymin = conf.low*100, ymax = conf.high*100), fill = fem_color, alpha = 0.4) +
            scale_x_continuous(breaks = seq(0, 30, 5)) +
            coord_cartesian(xlim = c(12.8, 26)) +
            scale_y_continuous(breaks = seq(5, 95, 10), limits = c(34, 58)) +
            theme_classic() +
            theme(panel.grid = element_blank(), axis.text = element_text(size = 12),
                  axis.title = element_text(size = 14)) +
            labs(x = "Daytime Temperature (°C) \n Nestling Ages 5-12 Year n", y = "Percent Returning Year n+1")
          
          ggsave(here::here("output_plots", "female_survive_plot.png"), fsurv_plot, device = "png", width = 4, height = 4, units = "in", dpi = 300)
          
          
          
# Variation in behavioral sensitivity to challenges: population incubation ----
                        
        # join incubation bout data to nest information and filter to useable nests
          # filter nests down to the needed columns
              nests_inc <- nests[nests$site == "Unit_1" | nests$site == "Unit_2" | nests$site == "Unit_4", ]
              nests_inc$unit <- gsub("Unit_", "", nests_inc$site)
              nests_inc$year <- nests_inc$exp_year
              nests_inc <- nests_inc[, c("nest_key", "year", "unit", "nest", "attempt_num", "clutch_init_doy",
                                         "clutch_comp_doy", "clutch_size", "hatch_doy", "brood_size_hatching", 
                                         "day6_brood_mass", "brood_size_day6", "nest_fate_doy", "simple_fate",
                                         "num_fledged", "nest_experiment", "nest_treatment", "female_id", "male_id")]
        
            # join to bouts
                all_bouts_ex$year <- as.integer(all_bouts_ex$year)
                all_bouts <- all_bouts_ex %>%
                  left_join(nests_inc, by = c("unit", "nest", "year"), relationship = "many-to-many")
                
            # need to filter out matches to correct attempt number for nests with two attempts in same year
                all_bouts$day_diff <- all_bouts$JDate - all_bouts$clutch_comp_doy
                
            # filter out the day of capture when incubation could be disrupted
                cap_sub <- ith_capture %>%
                  filter(adult_or_nestling == "Adult", sex == "Female", exp_year > 2012, location == "Ithaca") %>%
                  dplyr::select(nest_key, encounter_doy) 
                colnames(cap_sub) <- c("nest_key", "JDate")
                cap_sub$captured_this_day <- "yes"
                
                all_bouts <- plyr::join(all_bouts, cap_sub, c("nest_key", "JDate"), "left", "first")
                all_bouts <- all_bouts %>%
                  filter(is.na(captured_this_day) == TRUE)
            
            # only saving bouts that are between 4 to 17 days after clutch completion to avoid brooding 
                # and because hobos were typically installed on day 3 of incubation
                # the other joins are either date errors or matched to the wrong nesting attempt
                all_bouts <- all_bouts %>% filter(day_diff > 3, day_diff < 17)
                
            # now join to captures so I can get individual treatments
                all_bouts <- plyr::join(all_bouts, ith_adults[, c("nest_key", "individual_treatment")], 
                                        "nest_key", "left", "first")
                
            # add current hourly temperature
                all_bouts$yr_day_hr <- paste(all_bouts$year, all_bouts$JDate, all_bouts$hour, sep = "_")
                all_bouts <- plyr::join(all_bouts, w_hour, "yr_day_hr", "left", "first") 
                
            # remove some duplicated columns caused by joining
                all_bouts <- all_bouts[, !duplicated(names(all_bouts))]
                
            # remove nighttime bouts because they are typically incubating all night and these aren't scored correctly
                all_bouts <- all_bouts[all_bouts$is_night == 0, ]
                
            # remove bouts that were scored as being in periods that the thermocouple wasn't picking up well
                all_bouts <- all_bouts[all_bouts$exclude_all == 0, ]
                
            # remove extremely long bouts (80 minutes), these are generally scoring errors or bad thermocouple readings
                all_bouts <- all_bouts[all_bouts$duration / 60 < 80, ]
                all_bouts <- all_bouts[is.na(all_bouts$type) == FALSE, ]
                
            # determine how many bouts and range of temperatures are observed for each nest
                ab_stat <- all_bouts %>%
                  dplyr::group_by(nest_key) %>%
                  summarise(count = n(), min_temp = min(temp_C, na.rm = TRUE), max_temp = max(temp_C, na.rm = TRUE))
                ab_stat$temp_range <- ab_stat$max_temp - ab_stat$min_temp
                
                all_bouts <- plyr::join(all_bouts, ab_stat, "nest_key", "left", "first")
                
            # filter out nests that have fewer than 50 bouts scored and less than a 10 degree C temperature range
                all_bouts <- all_bouts[all_bouts$count > 49 & all_bouts$temp_range > 9, ]
                
                all_bouts <- all_bouts[, !duplicated(names(all_bouts))]
                all_bouts$type <- gsub("off", "Off", all_bouts$type)
                all_bouts$type <- gsub("on", "On", all_bouts$type)
              
            # exclud by treatment
                excludes <- exclude_trt$treatments[exclude_trt$exclude_female == "yes"]
                cort_trts <- c("Cort", "CORT_6d", "High", "high")
                all_bouts_cort <- all_bouts[all_bouts$individual_treatment %in% cort_trts, ]
                all_bouts <- all_bouts[!(all_bouts$nest_treatment %in% excludes), ]
                all_bouts <- all_bouts[!(all_bouts$individual_treatment %in% excludes), ]
                  
            # summarize by temperature
                
                all_bouts$round_temp <- round(all_bouts$temp_C, 0)
                avg_bouts <- all_bouts %>%
                  dplyr::group_by(type, round_temp) %>%
                  dplyr::summarise(avg_duration = mean(duration / 60), sd_duration = sd(duration / 60), count_bouts = n())
                avg_bouts$sem_bout <- avg_bouts$sd_duration / sqrt(avg_bouts$count_bouts)
                avg_bouts <- na.omit(avg_bouts)
             
            # figure showing population level incubation response to temperature   
                pop_inc <- ggplot(all_bouts, mapping = aes(x = temp_C)) + 
                    geom_histogram(color = "gray75", fill = "gray75", aes(y = ..density.. * 130 + 8.5), binwidth = 1, alpha = 0.2) +
                    annotate("rect", xmin = -10, xmax = 50, ymin = -10, ymax = 8.1, fill = "white") +
                    geom_smooth( mapping = aes(y = duration / 60, fill = type, color = type)) +
                    scale_color_manual(values = c("slateblue", "coral3")) +
                    scale_fill_manual(values = c("slateblue", "coral3")) +
                    theme_classic() +
                    xlim(c(4, 33)) +
                    #geom_segment(aes(y = 10.1, yend = 10, x = 2.5, xend = 32.5), color = "gray75", size = .8) +
                    coord_cartesian(ylim = c(9, 22), xlim = c(5, 33)) +
                    theme(panel.grid = element_blank(), axis.title = element_text(size = 14), axis.text = element_text(size = 12)) +
                    theme(legend.position = c(0.85, 0.85), legend.title = element_blank()) +
                    labs(x = "Ambient Temperature (C)", y = "Bout Duration (minutes)") #+
                    #geom_segment(data = avg_bouts[avg_bouts$round_temp > 4, ],
                    #             mapping = aes(x = round_temp, xend = round_temp, y = avg_duration - sem_bout, yend = avg_duration + sem_bout, color = type)) +
                    #geom_point(data = avg_bouts[avg_bouts$round_temp > 4, ], 
                    #           mapping = aes(x = round_temp, y = avg_duration, fill = type), shape = 21, color = "black")
                
                ggsave(here::here("output_plots", "fig2a.svg"), pop_inc, device = "svg", width = 4.5, height = 3.2, units = "in")

# Variation in behavioral sensitivity to challenges: individual incubation ----
              # Now plot linear response to declining temperature for each nest
                  # set up the categories in the order I want
                    all_bouts$type2 <- factor(all_bouts$type, levels = c("On", "Off"))
                    
              # filter out repeat observations
                  # out of 209 nests, there are only 25 females observed in multiple nests
                  # this isn't enough to model id + nest attempt separately, so I will randomly exclude all but one obs per female
                  # I need to wrangle some new columns to do that. This cuts down to 184 nests
                      # make a new dataframe of unique nests
                        ind_fems <- data.frame(nest_key = unique(all_bouts$nest_key))
                        ind_fems <- plyr::join(ind_fems, all_bouts[, c("nest_key", "female_id")], "nest_key", "left", "first")
                        for(i in 1:nrow(ind_fems)){
                          ind_fems$nest_count[i] <- nrow(ind_fems[ind_fems$female_id == ind_fems$female_id[i], ])
                        }
                      # randomly select one to include for females with more than one nest
                        set.seed(896)
                        ind_fems <- ind_fems %>%
                          dplyr::group_by(female_id) %>%
                          dplyr::mutate(include = ifelse(nest_count == 1, "yes", "no")) %>%
                          dplyr::mutate(include = ifelse(nest_count > 1 & row_number() == sample(1:n(), 1), "yes", include)) %>%
                          ungroup()
                      # join the random exclusion back to the main dataframe
                        all_bouts <- plyr::join(all_bouts, ind_fems[, c("nest_key", "include")], "nest_key", "left", "first")
                        
                # make a basic model to partition variance
                    all_bouts <- all_bouts[, !duplicated(names(all_bouts))]
                    ab_fem <- all_bouts %>% filter(include == "yes")
                    
                # for on bouts
                    # filter and scale
                      ab_fon <- ab_fem %>%
                        filter(is.na(temp_C) == FALSE)
                      ab_fon$temp_s <- scale(ab_fon$temp_C) 
                      ab_fon$duration_s <- scale(ab_fon$duration/60)
                      ab_fon$duration_m <- ab_fon$duration/60
                      ab_fon <- ab_fon[ab_fon$duration_m > 3, ]
                    
                    # fit a simple model
                      ab_fon1 <- ab_fon[ab_fon$type == "On" & ab_fon$temp_C < 19, ]
                      mi_on <- lmer(duration_m ~ temp_s + (temp_s|nest_key), data = ab_fon1)
                      mi_on_drop <- lmer(duration_m ~ temp_s + (1|nest_key), data = ab_fon1)
                      temp_mean <- mean(ab_fon$temp_C)
                      temp_sd <- sd(ab_fon$temp_C)
                    
                    # exctract random effects and convert back to measurement scale  
                      re_mion <- ranef(mi_on)$nest_key
                      colnames(re_mion) <- c("r_intercept", "r_slope")
                      re_mion$r_intercept <- fixef(mi_on)[1] + re_mion$r_intercept - (fixef(mi_on)[2] * temp_mean / temp_sd)
                      re_mion$r_slope <- (fixef(mi_on)[2] + re_mion$r_slope) / temp_sd
                      re_mion$x_start <- 3
                      re_mion$x_end <- 19
                      re_mion$y_start <- re_mion$r_slope * 3 + re_mion$r_intercept
                      re_mion$y_end <- re_mion$r_slope * 19 + re_mion$r_intercept 
                      
                # for off bouts
                      
                    # fit a simple model
                      ab_fon2 <- ab_fon[ab_fon$type == "Off" & ab_fon$temp_C < 19, ]
                      mi_off <- lmer(duration_m ~ temp_s + (temp_s|nest_key), data = ab_fon2)
                      mi_off_drop <- lmer(duration_m ~ temp_s + (1|nest_key), data = ab_fon2)
                      
                    # exctract random effects and convert back to measurement scale  
                      re_mioff <- ranef(mi_off)$nest_key
                      colnames(re_mioff) <- c("r_intercept", "r_slope")
                      re_mioff$r_intercept <- fixef(mi_off)[1] + re_mioff$r_intercept - (fixef(mi_off)[2] * temp_mean / temp_sd)
                      re_mioff$r_slope <- (fixef(mi_off)[2] + re_mioff$r_slope) / temp_sd
                      re_mioff$x_start <- 3
                      re_mioff$x_end <- 19
                      re_mioff$y_start <- re_mioff$r_slope * 3 + re_mioff$r_intercept
                      re_mioff$y_end <- re_mioff$r_slope * 19 + re_mioff$r_intercept    
                      
                    fig2c_df <- data.frame(x_start = c(re_mion$x_start, re_mioff$x_start),
                                           xend = c(re_mion$x_end, re_mioff$x_end),
                                           y_start = c(re_mion$y_start, re_mioff$y_start),
                                           y_end = c(re_mion$y_end, re_mioff$y_end),
                                           type = factor(c(rep("On", nrow(re_mion)), rep("Off", nrow(re_mioff))), levels = c("On", "Off")))
                    
                # make a figure
                    mu_inc_effect <- fig2c_df %>%
                      dplyr::group_by(type) %>%
                      dplyr::summarise(mu_x = mean(x_start), mu_xend = mean(xend),
                                       mu_y = mean(y_start), mu_yend = mean(y_end))
                    
                    ind_var_inc <- ggplot(fig2c_df) +
                      geom_segment(mapping = aes(x = x_start, xend = xend, y = y_start, yend = y_end, color = type),
                                   size = 0.2, alpha = 0.4) +
                      facet_wrap(~ type) +
                      theme_classic() +
                      theme(panel.grid = element_blank()) +
                      scale_color_manual(values = c("coral3", "slateblue")) +
                      guides(color = "none") +
                      geom_vline(xintercept = 19, size = 0.8, linetype = "dashed", color = "gray40") +
                      theme(strip.text = element_blank(), axis.text = element_text(size = 12), axis.title = element_text(size = 14)) +
                      coord_cartesian(xlim = c(3, 21), ylim = c(0, 34)) +
                      scale_y_continuous(breaks = seq(-10, 50, 10)) +
                      scale_x_continuous(breaks = seq(-10, 40, 10)) +
                      geom_segment(data = mu_inc_effect, aes(x = mu_x, xend = mu_xend, y = mu_y, yend = mu_yend), size = 1.1, color = "black") +
                      labs(x = "temperature", y = "bout duration")
                    
                  
                ggsave(here::here("output_plots", "fig2c.svg"), ind_var_inc, device = "svg", width = 3, height = 2.2, units = "in") 

                # make plot of random effect covariation
                      r1i <- ranef(mi_on)$nest_key
                      colnames(r1i) <- c("Intercept", "Slope")
                      r1i$type <- "On"
                      r1i$nest_key <- rownames(r1i)
                      
                      r2i <- ranef(mi_off)$nest_key
                      colnames(r2i) <- c("Intercept", "Slope")
                      r2i$type <- "Off"
                      r2i$nest_key <- rownames(r2i)
                      
                      inc_rf <- data.frame(rbind(r1i, r2i))
                      
                      inc_covar <- ggplot(inc_rf, mapping = aes(x = Intercept, y = Slope, color = type, fill = type)) +
                        geom_point(size = 0.3, alpha = 0.6) +
                        geom_smooth(method = "lm", alpha = 0.4, size = 0.5) +
                        theme_classic() +
                        theme(panel.grid = element_blank(), axis.text = element_text(size = 12)) +
                        scale_color_manual(values = c("slateblue", "coral3")) +
                        scale_fill_manual(values = c("slateblue", "coral3")) +
                        guides(fill = "none", color = "none") +
                        scale_x_continuous(breaks = seq(-12, 12, 6)) +
                        coord_cartesian(ylim = c(-2.2, 2.2), xlim = c(-7, 10)) +
                        scale_y_continuous(breaks = seq(-6, 6, 2))
                      
                      ggsave(here::here("output_plots", "fig2f.svg"), inc_covar, device = "svg", width = 1.8, height = 1.8, units = "in")
                
                # make barplot of variance explained
                inc_var <- data.frame(type = c(rep("On", 2), rep("Off", 2)),
                                      measure = rep(c("Intercept", "Slope"), 2),
                                      value = c(9.9, 3.2, 7.4, 4.2))
                inc_var$type = factor(inc_var$type, levels = c("Off", "On"))
                inc_var$measure = factor(inc_var$measure, levels = c("Intercept", "Slope"))
                inc_varplot <- ggplot(inc_var, mapping = aes(y = type, x = value, fill = measure)) +
                  geom_bar(stat = "identity", aes(fill = type, color = type), alpha = 0.4, width = 0.6) +
                  scale_fill_manual(values = c("slateblue", "coral3")) +
                  scale_color_manual(values = c("slateblue", "coral3")) +
                  guides(fill = "none", color = "none") +
                  theme_classic() +
                  theme(axis.text.y = element_blank()) +
                  theme(panel.grid = element_blank(), axis.text = element_text(size = 12)) +
                  #theme(axis.line.x = element_blank()) +
                  scale_x_continuous(breaks = seq(0, 40, 10), expand = c(0, 0), limits = c(0, 27))
                ggsave(here::here("output_plots", "fig2e.svg"), inc_varplot, device = "svg", width = 3, height = 1.4, units = "in")
                
          
# Variation in behavioral sensitivity to challenges: population provisioning ----
    # join feeding bout data to nest information and filter to usable nests
          # filter nests down to the needed columns
              nests_fd <- nests[nests$site == "Unit_1" | nests$site == "Unit_2" | nests$site == "Unit_4" | nests$site == "Turkey_Hill", ]
              nests_fd$unit <- gsub("Unit_", "", nests_fd$site)
              nests_fd$unit <- gsub("Turkey_Hill", "T", nests_fd$unit)
              nests_fd$year <- nests_fd$exp_year
              nests_fd <- nests_fd[, c("nest_key", "year", "unit", "nest", "attempt_num", "clutch_init_doy",
                                         "clutch_comp_doy", "clutch_size", "hatch_doy", "brood_size_hatching", 
                                         "day6_brood_mass", "brood_size_day6", "nest_fate_doy", "simple_fate",
                                         "num_fledged", "nest_experiment", "nest_treatment", "female_id", "male_id")]
              
            # make a unit column for feeding
              feeds2 <- feeds %>%
                separate(unitbox, into = c("unit", "nest"), sep = "_")
              feeds2$unit <- gsub("TH", "T", feeds2$unit)
        
            # join to feeding data
                all_feeds <- feeds2 %>%
                  left_join(nests_fd, by = c("unit", "nest", "year"), relationship = "many-to-many")
                
            # need to filter out matches to correct attempt number for nests with two attempts in same year
                all_feeds$day_diff <- all_feeds$yday - all_feeds$hatch_doy
            
            # only saving feeding after hatching day and up to day 14 (because sometimes nestling pits on d15 interfere)
                # the other joins are either date errors or matched to the wrong nesting attempt
                all_feeds <- all_feeds %>% filter(day_diff > 0, day_diff < 15)
                
            # now join to captures so I can get individual treatments
                all_feeds <- plyr::join(all_feeds, ith_adults[, c("nest_key", "individual_treatment")], 
                                        "nest_key", "left", "first")
                
            # add current hourly temperature
                all_feeds$yr_day_hr <- paste(all_feeds$year, all_feeds$yday, all_feeds$Hour, sep = "_")
                all_feeds <- plyr::join(all_feeds, w_hour, "yr_day_hr", "left", "first") 
                
            # remove some duplicated columns caused by joining
                all_feeds <- all_feeds[, !duplicated(names(all_feeds))]
                
            # further filtering
                all_feeds <- all_feeds %>%
                  filter(Hour < 21 & Hour > 5)
                
            # group to get averages for plotting    
                all_feeds$round_temp <- round(all_feeds$temp_C, 0)
                avg_feeds <- all_feeds %>%
                  dplyr::group_by(round_temp) %>%
                  dplyr::summarise(avg_fem = mean(f_feeds, na.rm = TRUE), avg_mal = mean(m_feeds, na.rm = TRUE),
                            sd_fem = sd(f_feeds, na.rm = TRUE), sd_mal = sd(m_feeds, na.rm = TRUE),
                            count = n())
                  avg_feeds$sem_fem <- avg_feeds$sd_fem / sqrt(avg_feeds$count)
                  avg_feeds$sem_mal <- avg_feeds$sd_mal / sqrt(avg_feeds$count)
                  
              # add number alive on d12 banding to potentially filter out nests that failed early
                  d12num <- data.frame(nest_key = unique(nests$nest_key[nests$exp_year > 2012 & nests$location == "Ithaca"]))
                  for(i in 1:nrow(d12num)){
                    temp <- ith_capture %>%
                      filter(adult_or_nestling == "Nestling", age == 12, nest_key == d12num$nest_key[i])
                    d12num$numd12[i] <- nrow(temp)
                  }
                  all_feeds <- plyr::join(all_feeds, d12num, "nest_key", "left", "first")
                  
              # filter out any records after recorded nest failure date
                  all_feeds <- all_feeds %>%
                    filter(yday < nest_fate_doy)
                  
              # add temperature data at different intervals to feeding data
                # create a frame with each combination of year and yday feeding observations occurred on
                   yr_doy_uni2 <- all_feeds %>%
                     dplyr::group_by(year, yday) %>%
                     summarise(count = n()) %>%
                     as.data.frame()
                   
                # for each year/day combination, calculate various lagged temperature averages (note not all used) 
                   for(i in 1:nrow(yr_doy_uni2)){
                     sub1 <- subset(wh_day, wh_day$year == yr_doy_uni2$year[i] &
                                      wh_day$yday > yr_doy_uni2$yday[i] - 20 &
                                      wh_day$yday < yr_doy_uni2$yday[i] + 1)
                     yr_doy_uni2$prior1_C[i] <- mean(c(sub1$temp_C[sub1$yday > (yr_doy_uni2$yday[i] - 2) &
                                                                    sub1$yday < yr_doy_uni2$yday[i]]))
                     yr_doy_uni2$prior2_C[i] <- mean(c(sub1$temp_C[sub1$yday > (yr_doy_uni2$yday[i] - 3) &
                                                                    sub1$yday < yr_doy_uni2$yday[i]]))
                     yr_doy_uni2$prior3_C[i] <- mean(c(sub1$temp_C[sub1$yday > (yr_doy_uni2$yday[i] - 4) &
                                                                    sub1$yday < yr_doy_uni2$yday[i]]))
                     yr_doy_uni2$prior4_C[i] <- mean(c(sub1$temp_C[sub1$yday > (yr_doy_uni2$yday[i] - 5) &
                                                                    sub1$yday < yr_doy_uni2$yday[i]]))
                     yr_doy_uni2$prior5_C[i] <- mean(c(sub1$temp_C[sub1$yday > (yr_doy_uni2$yday[i] - 6) &
                                                                    sub1$yday < yr_doy_uni2$yday[i]]))
                     yr_doy_uni2$prior6_C[i] <- mean(c(sub1$temp_C[sub1$yday > (yr_doy_uni2$yday[i] - 7) &
                                                                    sub1$yday < yr_doy_uni2$yday[i]]))
                     yr_doy_uni2$prior7_C[i] <- mean(c(sub1$temp_C[sub1$yday > (yr_doy_uni2$yday[i] - 8) &
                                                                    sub1$yday < yr_doy_uni2$yday[i]]))
                   }
                   
                # join the lagged temperature captures to the feeding data
                   yr_doy_uni2$year_obsday <- paste(yr_doy_uni2$year, yr_doy_uni2$yday, sep = "_")
                   all_feeds$year_obsday <- paste(all_feeds$year, all_feeds$yday, sep = "_")
                   all_feeds <- plyr::join(all_feeds, yr_doy_uni2, "year_obsday", "left", "first")
                   all_feeds <- all_feeds[, !duplicated(names(all_feeds))]
              
                # figure out what days individuals were captured and exclude them on that day for models
                   # females
                   fem_check <- ith_capture %>% filter(exp_year > 2012, adult_or_nestling == "Adult", sex == "Female")
                   fem_check <- fem_check[, c("band", "encounter_doy")]
                   colnames(fem_check) <- c("female_id", "yday")
                   fem_check$female_capture_day <- "yes"
                   
                   all_feeds <- plyr::join(all_feeds, fem_check, c("female_id", "yday"), "left", "first")
                   all_feeds[is.na(all_feeds$female_capture_day) == TRUE, "female_capture_day"] <- "no"   
                   
                   # males
                   mal_check <- ith_capture %>% filter(exp_year > 2012, adult_or_nestling == "Adult", sex == "Male")
                   mal_check <- mal_check[, c("band", "encounter_doy")]
                   colnames(mal_check) <- c("male_id", "yday")
                   mal_check$male_capture_day <- "yes"
                   
                   all_feeds <- plyr::join(all_feeds, mal_check, c("male_id", "yday"), "left", "first")
                   all_feeds[is.na(all_feeds$male_capture_day) == TRUE, "male_capture_day"] <- "no"   
                   
              # exclude based on treatment
                   excludes <- exclude_trt$treatments[exclude_trt$exclude_female == "yes"]
                   all_feeds_save <- all_feeds #save version before exclusions
                   all_feeds <- all_feeds[!(all_feeds$nest_treatment %in% excludes), ]
                   all_feeds <- all_feeds[!(all_feeds$individual_treatment %in% excludes), ]
                   
              # make the plot  
                pop_feed <- all_feeds %>%
                  #filter(Offset > 5, Offset < 12, Hour > 7, Hour < 16) %>%
                ggplot(mapping = aes(x = temp_C)) +
                  geom_histogram(color = "gray75", fill = "gray75", aes(y = ..density.. * 110 + 3.5), binwidth = 1, alpha = 0.2) +
                  annotate("rect", xmin = -10, xmax = 50, ymin = -10, ymax = 3.5, fill = "white") +
                  geom_smooth(color = fem_color, fill = fem_color, aes(y = f_feeds), formula = y ~ s(x, k = 7)) +
                  geom_smooth(mapping = aes(y = m_feeds), color = mal_color, fill = mal_color) +
                  xlim(c(7,32)) +
                  coord_cartesian(ylim = c(4,14)) +
                  theme_classic() +
                  theme(panel.grid = element_blank(), axis.text = element_text(size = 12), axis.title = element_text(size = 16)) +
                  geom_segment(data = avg_feeds, mapping = aes(x = round_temp, xend = round_temp, 
                                                               y = avg_fem - sem_fem, yend = avg_fem + sem_fem), color = fem_color) +
                  geom_point(data = avg_feeds, mapping = aes(x = round_temp, y = avg_fem), shape = 21, fill = fem_color) +
                  geom_segment(data = avg_feeds, mapping = aes(x = round_temp, xend = round_temp, 
                                                               y = avg_mal - sem_mal, yend = avg_mal + sem_mal), color = mal_color) +
                  geom_point(data = avg_feeds, mapping = aes(x = round_temp, y = avg_mal), shape = 21, fill = mal_color) +
                  labs(x = "Ambient Temperature (C)", y = "Feeding Trips Per Hour")
                
                ggsave(here::here("output_plots", "fig2b.svg"), pop_feed, device = "svg", width = 4.5, height = 3.2, units = "in")
            
                    
# Variation in behavioral sensitivity to challenges: individual provisioning ----     
        
   
                 
              # filter out repeat observations
                  # out of 376 nests, there are 314 unique females and 206 unique males (lots of nests don't have male data)
                  # this isn't enough to model id + nest attempt separately, so I will randomly exclude all but one obs per individual
                  # this needs to be done separately for males and females unlike for incubation above
                  # I need to wrangle some new columns to do that
                      # make a new dataframe of unique nests
                        ind_ffems <- data.frame(nest_key = unique(all_feeds$nest_key))
                        ind_ffems <- plyr::join(ind_ffems, all_feeds[, c("nest_key", "female_id")], "nest_key", "left", "first")
                        for(i in 1:nrow(ind_ffems)){
                          ind_ffems$nest_count[i] <- nrow(ind_ffems[ind_ffems$female_id == ind_ffems$female_id[i], ])
                        }
                        
                        ind_fmals <- data.frame(nest_key = unique(all_feeds$nest_key[is.na(all_feeds$m_feeds) == FALSE]))
                        ind_fmals <- plyr::join(ind_fmals, all_feeds[, c("nest_key", "male_id")], "nest_key", "left", "first")
                        for(i in 1:nrow(ind_fmals)){
                          ind_fmals$nest_count[i] <- nrow(ind_fmals[ind_fmals$male_id == ind_fmals$male_id[i], ])
                        }
                      # randomly select one to include for females and males with more than one nest
                        set.seed(625)
                        ind_ffems <- ind_ffems %>%
                          dplyr::group_by(female_id) %>%
                          dplyr::mutate(include = ifelse(nest_count == 1, "yes", "no")) %>%
                          dplyr::mutate(include = ifelse(nest_count > 1 & row_number() == sample(1:n(), 1), "yes", include)) %>%
                          ungroup()
                        
                        ind_fmals <- ind_fmals %>%
                          dplyr::group_by(male_id) %>%
                          dplyr::mutate(include = ifelse(nest_count == 1, "yes", "no")) %>%
                          dplyr::mutate(include = ifelse(nest_count > 1 & row_number() == sample(1:n(), 1), "yes", include)) %>%
                          ungroup()
                      # join the random exclusion back to the main dataframe
                        all_feedsf <- plyr::join(all_feeds, ind_ffems[, c("nest_key", "include")], "nest_key", "left", "first")
                        all_feedsm <- plyr::join(all_feeds, ind_fmals[, c("nest_key", "include")], "nest_key", "left", "first")
                        
                # make a basic model to partition variance
                        all_feedsf <- all_feedsf[, !duplicated(names(all_feedsf))]
                        all_feedsm <- all_feedsm[, !duplicated(names(all_feedsm))]
                        af_fem <- all_feedsf %>% filter(include == "yes", female_capture_day == "no", f_feeds < 37)
                        af_mal <- all_feedsm %>% filter(include == "yes", male_capture_day == "no", m_feeds < 37)
                    
                # for feeding rate
                    # filter and scale
                      #females
                        af_ffeed <- af_fem %>%
                          filter(temp_C < 25, is.na(temp_C) == FALSE)
                        af_ffeed$temp_s <- scale(af_ffeed$temp_C) 
                        af_ffeed$hour_s <- scale(af_ffeed$Hour) 
                        af_ffeed$age_s <- scale(af_ffeed$Offset)
                        af_ffeed$nest_key <- as.factor(af_ffeed$nest_key)
                      #males
                        af_mfeed <- af_mal %>%
                          filter(temp_C < 25, is.na(temp_C) == FALSE)
                        af_mfeed$temp_s <- scale(af_mfeed$temp_C) 
                        af_mfeed$hour_s <- scale(af_mfeed$Hour) 
                        af_mfeed$age_s <- scale(af_mfeed$Offset)
                    
                    # fit a simple model as basis for extracting random effects
                        # more complex models with other variables were tried but the random effect estimates
                        # are very similar and this is easier to understand
                      # females
                        mffeed <- lmer(f_feeds ~ temp_s + (temp_s|nest_key), data = af_ffeed[af_ffeed$female_capture_day == "no", ])
                        mffeed_drop <- lmer(f_feeds ~ temp_s + (1|nest_key), data = af_ffeed[af_ffeed$female_capture_day == "no", ])
                        temp_meanf <- mean(af_ffeed$temp_C)
                        temp_sdf <- sd(af_ffeed$temp_C)
                        
                        # inspect model. with huge sample size, tests are significant but visual residuals look fine
                          # sf <- simulateResiduals(mffeed)
                          # plot(sf)
                          # testDispersion(sf)
                          # hist(residuals(mffeed))
                        
                      # males
                        mmfeed <- lmer(m_feeds ~ temp_s + (temp_s|nest_key), data = af_mfeed)
                        mmfeed_drop <- lmer(m_feeds ~ temp_s + (1|nest_key), data = af_mfeed)
                        temp_meanm <- mean(af_mfeed$temp_C)
                        temp_sdm <- sd(af_mfeed$temp_C)
                        
                        # inspect model. with huge sample size, tests are significant but visual residuals look fine
                          # sf <- simulateResiduals(mmfeed)
                          # plot(sf)
                          # testDispersion(sf)
                          # hist(residuals(mmfeed))
                    
                    # extract random effects and convert back to measurement scale  
                        # females
                          re_ffeed <- ranef(mffeed)$nest_key  
                          colnames(re_ffeed) <- c("r_intercept", "r_slope")
                          re_ffeed$r_intercept <- fixef(mffeed)[1] + re_ffeed$r_intercept - (fixef(mffeed)[2] * temp_meanf / temp_sdf)
                          re_ffeed$r_slope <- (fixef(mffeed)[2] + re_ffeed$r_slope) / temp_sdf
                          re_ffeed$x_start <- 6
                          re_ffeed$x_end <- 24
                          re_ffeed$y_start <- re_ffeed$r_slope * 6 + re_ffeed$r_intercept
                          re_ffeed$y_end <- re_ffeed$r_slope * 24 + re_ffeed$r_intercept
                          
                        # males
                          re_mfeed <- ranef(mmfeed)$nest_key  
                          colnames(re_mfeed) <- c("r_intercept", "r_slope")
                          re_mfeed$r_intercept <- fixef(mmfeed)[1] + re_mfeed$r_intercept - (fixef(mmfeed)[2] * temp_meanm / temp_sdm)
                          re_mfeed$r_slope <- (fixef(mmfeed)[2] + re_mfeed$r_slope) / temp_sdm
                          re_mfeed$x_start <- 6
                          re_mfeed$x_end <- 24
                          re_mfeed$y_start <- re_mfeed$r_slope * 6 + re_mfeed$r_intercept
                          re_mfeed$y_end <- re_mfeed$r_slope * 24 + re_mfeed$r_intercept
                          
                        # make data frame for plotting
                          fig2d_df <- data.frame(x_start = c(re_ffeed$x_start, re_mfeed$x_start),
                                                 xend = c(re_ffeed$x_end, re_mfeed$x_end),
                                                 y_start = c(re_ffeed$y_start, re_mfeed$y_start),
                                                 y_end = c(re_ffeed$y_end, re_mfeed$y_end),
                                                 type = factor(c(rep("Female", nrow(re_ffeed)), rep("Male", nrow(re_mfeed))), levels = c("Female", "Male")))
                # make a figure
                    mu_feed_effect <- fig2d_df %>%
                      dplyr::group_by(type) %>%
                      dplyr::summarise(mu_x = mean(x_start), mu_xend = mean(xend),
                                       mu_y = mean(y_start), mu_yend = mean(y_end))
                    
                    ind_var_feed <- ggplot(fig2d_df) +
                      geom_segment(mapping = aes(x = x_start, xend = xend, y = y_start, yend = y_end, color = type),
                                   size = 0.2, alpha = 0.4) +
                      facet_wrap(~ type) +
                      theme_classic() +
                      theme(panel.grid = element_blank()) +
                      scale_color_manual(values = c(fem_color, mal_color)) +
                      guides(color = "none") +
                      geom_vline(xintercept = 24, size = 0.8, linetype = "dashed", color = "gray40") +
                      theme(strip.text = element_blank(), axis.text = element_text(size = 12), axis.title = element_text(size = 14)) +
                      coord_cartesian(xlim = c(6, 26), ylim = c(0, 40)) +
                      scale_y_continuous(breaks = seq(-10, 50, 10)) +
                      scale_x_continuous(breaks = seq(-10, 40, 10)) +
                      geom_segment(data = mu_feed_effect, aes(x = mu_x, xend = mu_xend, y = mu_y, yend = mu_yend), size = 1.1, color = "black") +
                      labs(x = "temperature", y = "feed trips")                
                
                
                
                
                  ggsave(here::here("output_plots", "fig2D.svg"), ind_var_feed, device = "svg", width = 3, height = 2.2, units = "in")
              
                
            # make plot of random effect covariation
                r1 <- ranef(mffeed)$nest_key
                colnames(r1) <- c("Intercept", "Slope")
                r1$sex <- "Female"
                r1$nest_key <- rownames(r1)
                
                r2 <- ranef(mmfeed)$nest_key
                colnames(r2) <- c("Intercept", "Slope")
                r2$sex <- "Male"
                r2$nest_key <- rownames(r2)
                
                feed_rf <- data.frame(rbind(r1, r2))
                
                feed_covar <- ggplot(feed_rf, mapping = aes(x = Intercept, y = Slope, color = sex, fill = sex)) +
                  geom_point(size = 0.3, alpha = 0.5) +
                  geom_smooth(method = "lm", alpha = 0.4, size = 0.5) +
                  theme_classic() +
                  theme(panel.grid = element_blank(), axis.text = element_text(size = 12)) +
                  scale_color_manual(values = c(fem_color, mal_color)) +
                  scale_fill_manual(values = c(fem_color, mal_color)) +
                  guides(fill = "none", color = "none") +
                  scale_x_continuous(breaks = seq(-12, 12, 6))
                  
                ggsave(here::here("output_plots", "fig2h.svg"), feed_covar, device = "svg", width = 1.8, height = 1.8, units = "in")
                
            # make barplot of variance explained
                feed_var <- data.frame(sex = c(rep("Female", 2), rep("Male", 2)),
                                       measure = rep(c("Intercept", "Slope")),
                                       value = c(18.6, 3.6, 21.3, 3.2))
                feed_var$sex = factor(feed_var$sex, levels = c("Male", "Female"))
                feed_var$measure = factor(feed_var$measure, levels = c("Intercept", "Slope"))
                feed_varplot <- ggplot(feed_var, mapping = aes(y = sex, x = value, fill = measure)) +
                  geom_bar(stat = "identity", aes(fill = sex, color = sex), alpha = 0.4, width = 0.6) +
                  scale_fill_manual(values = c(mal_color, fem_color)) +
                  scale_color_manual(values = c(mal_color, fem_color)) +
                  guides(fill = "none", color = "none") +
                  theme_classic() +
                  theme(axis.text.y = element_blank()) +
                  theme(panel.grid = element_blank(), axis.text = element_text(size = 12)) +
                  #theme(axis.line.x = element_blank()) +
                  scale_x_continuous(breaks = c(0, 10, 20, 30), expand = c(0, 0), limits = c(0, 27))
                ggsave(here::here("output_plots", "fig2g.svg"), feed_varplot, device = "svg", width = 3, height = 1.4, units = "in")
              

# Variation to lagged temperature: incubation + provisioning ----
  # This is repeating a similar model to those above but using prior day temperature
      # response to prior day temperature (accounting for current temperature) is interpreted
      # as variation in resilience. because prior day temperature is the same for all hourly
      # recordings on the focal day, these models also include a random effect for day of year
       
      # first for incubation bouts
                # join prior temp
                    ab_p <- ab_fon
                    ab_p$year_capday <- paste(ab_p$yr_doy)
                    ab_p <- plyr::join(ab_p, yr_doy_uni1[, c("year_capday", "prior1_C")], "year_capday", "left", "first")
                
                # on bout lagged
                      # fit a simple lagged model
                            ab_fon_lag <- ab_p[ab_p$type == "On" & ab_p$prior1_C < 19, ]
                            ab_fon_lag$p1C_s <- scale(ab_fon_lag$prior1_C)
                            mi_on_lag <- lmer(duration_m ~ temp_s + p1C_s + (p1C_s|nest_key) + (1|year_capday), data = ab_fon_lag)
                            mi_on_lag_drop <- lmer(duration_m ~ temp_s + p1C_s + (1|nest_key) + (1|year_capday), data = ab_fon_lag)
                            temp_mean_l <- mean(na.omit(ab_p$prior1_C))
                            temp_sd_l <- sd(na.omit(ab_p$prior1_C))
                            
                      # exctract random effects and convert back to measurement scale  
                            re_mion_lag <- ranef(mi_on_lag)$nest_key
                            colnames(re_mion_lag) <- c("r_intercept", "r_slope")
                            re_mion_lag$r_intercept <- fixef(mi_on_lag)[1] + re_mion_lag$r_intercept - (fixef(mi_on_lag)[2] * temp_mean / temp_sd)
                            re_mion_lag$r_slope <- (fixef(mi_on_lag)[2] + re_mion_lag$r_slope) / temp_sd
                            re_mion_lag$x_start <- 3
                            re_mion_lag$x_end <- 19
                            re_mion_lag$y_start <- re_mion_lag$r_slope * 3 + re_mion_lag$r_intercept
                            re_mion_lag$y_end <- re_mion_lag$r_slope * 19 + re_mion_lag$r_intercept      
                            
                # off bout lagged
                     # fit a simple lagged model  
                            ab_off_lag <- ab_p[ab_p$type == "Off" & ab_p$prior1_C < 19, ]
                            ab_off_lag$p1C_s <- scale(ab_off_lag$prior1_C)
                            mi_off_lag <- lmer(duration_m ~ temp_s + p1C_s + (p1C_s|nest_key) + (1|year_capday), data = ab_off_lag)
                            mi_off_lag_drop <- lmer(duration_m ~ temp_s + p1C_s + (1|nest_key) + (1|year_capday), data = ab_off_lag)
                            
                    # extract random effects and convert back to measurement scale
                            re_moff_lag <- ranef(mi_off_lag)$nest_key
                            colnames(re_moff_lag) <- c("r_intercept", "r_slope")
                            re_moff_lag$r_intercept <- fixef(mi_off_lag)[1] + re_moff_lag$r_intercept - (fixef(mi_off_lag)[2] * temp_mean / temp_sd)
                            re_moff_lag$r_slope <- (fixef(mi_off_lag)[2] + re_moff_lag$r_slope) / temp_sd
                            re_moff_lag$x_start <- 3
                            re_moff_lag$x_end <- 19
                            re_moff_lag$y_start <- re_moff_lag$r_slope * 3 + re_moff_lag$r_intercept
                            re_moff_lag$y_end <- re_moff_lag$r_slope * 19 + re_moff_lag$r_intercept 
                
      # next for provisioning behavior
                                      
                # fit the female lagged model
                  af_lag <- af_ffeed[af_ffeed$prior1_C < 24, ] 
                  af_lag$p1C_s <- scale(af_lag$prior1_C)
                  mffeed_lag <- lmer(f_feeds ~ temp_s + p1C_s + (p1C_s|nest_key) + (1|year_obsday), data = af_lag)
                  mffeed_lag_drop <- lmer(f_feeds ~ temp_s + p1C_s + (1|nest_key) + (1|year_obsday), data = af_lag)
                  
                # male lagged
                  am_lag <- af_mfeed[af_mfeed$prior1_C < 24, ]
                  am_lag$p1C_s <- scale(am_lag$prior1_C)
                  mmfeed_lag <- lmer(m_feeds ~ temp_s + p1C_s + (p1C_s|nest_key) + (1|year_obsday), data = am_lag,
                                     control = lmerControl(optimizer = "bobyqa"))
                  mmfeed_lag_drop <- lmer(m_feeds ~ temp_s + p1C_s + (1|nest_key) + (1|year_obsday), data = am_lag,
                                     control = lmerControl(optimizer = "bobyqa"))
                  
                  
                  # extract random effects and convert back to measurement scale  
                    # females
                      re_flag <- ranef(mffeed_lag)$nest_key  
                      colnames(re_flag) <- c("r_intercept", "r_slope")
                      re_flag$r_intercept <- fixef(mffeed_lag)[1] + re_flag$r_intercept - (fixef(mffeed_lag)[3] * temp_meanf / temp_sdf)
                      re_flag$r_slope <- (fixef(mffeed_lag)[3] + re_flag$r_slope) / temp_sdf
                      re_flag$x_start <- 6
                      re_flag$x_end <- 24
                      re_flag$y_start <- re_flag$r_slope * 6 + re_flag$r_intercept
                      re_flag$y_end <- re_flag$r_slope * 24 + re_flag$r_intercept
                      
                      
                      mu_feed_lag <- re_flag %>%
                        #dplyr::group_by(type) %>%
                        dplyr::summarise(mu_x = mean(x_start), mu_xend = mean(x_end),
                                         mu_y = mean(y_start), mu_yend = mean(y_end))
                      
                      ggplot(re_flag,aes()) + 
                        geom_segment(aes(x = x_start, xend = x_end, y = y_start, yend = y_end), alpha = 0.3, color = fem_color) +
                        theme_classic() + 
                        labs(x = "Prior Day Avg Temperature", y = "Current Feeding Rate") +
                        geom_segment(data = mu_feed_lag, aes(x = mu_x, xend = mu_xend, y = mu_y, yend = mu_yend), 
                                     color = "black", linewidth = 1.5)
                        
                  
                    # males
                      re_mlag <- ranef(mmfeed_lag)$nest_key  
                      colnames(re_mlag) <- c("r_intercept", "r_slope")
                      re_mlag$r_intercept <- fixef(mmfeed_lag)[1] + re_mlag$r_intercept - (fixef(mmfeed_lag)[3] * temp_meanm / temp_sdm)
                      re_mlag$r_slope <- (fixef(mmfeed_lag)[3] + re_mlag$r_slope) / temp_sdm
                      re_mlag$x_start <- 6
                      re_mlag$x_end <- 24
                      re_mlag$y_start <- re_mlag$r_slope * 6 + re_mlag$r_intercept
                      re_mlag$y_end <- re_mlag$r_slope * 24 + re_mlag$r_intercept
                      
                      
                      mu_feed_lag_m <- re_mlag %>%
                        #dplyr::group_by(type) %>%
                        dplyr::summarise(mu_x = mean(x_start), mu_xend = mean(x_end),
                                         mu_y = mean(y_start), mu_yend = mean(y_end))
                      
                      ggplot(re_mlag, aes()) +
                        geom_segment(aes(x = x_start, xend = x_end, y = y_start, yend = y_end), alpha = 0.3, color = mal_color) +
                        theme_classic() +
                        labs(x = "Prior Day Avg Temperature", y = "Current Feeding Rate") +
                        geom_segment(data = mu_feed_lag_m, aes(x = mu_x, xend = mu_xend, y = mu_y, yend = mu_yend),
                                     color = "black", linewidth = 1.5)
                      
                  
                  # make data frame for plotting
                  fig2d_df <- data.frame(x_start = c(re_ffeed$x_start, re_mfeed$x_start),
                                         xend = c(re_ffeed$x_end, re_mfeed$x_end),
                                         y_start = c(re_ffeed$y_start, re_mfeed$y_start),
                                         y_end = c(re_ffeed$y_end, re_mfeed$y_end),
                                         type = factor(c(rep("Female", nrow(re_ffeed)), rep("Male", nrow(re_mfeed))), levels = c("Female", "Male")))  
                  
                  ggplot()
      
# Model bootstrapping function to sample BLUPS from lmer ----
                  
        # using random effects directly from lmer moddels is anti-conservative because it doesn't
            # propagate uncertainty. I will sample 
              # lm_data should be a two column data frame with 'nest_key' and one predictor    
                # Or 3 columns with nest key and two predictors for binomial models (success/failure)
                # predictors should be standardized before loading if desired. random intercepts & slopes are
                # automatically standardized with scale()
              
                boot_blup <- function(blup_mod_name, lm_data, 
                                      which_re = c("intercept", "slope", "both_single", "both_inter",
                                                   "binom_sin", "binom_x", "binom_tmb_sin", "binom_tmb_x"), 
                                      boots = 10){
                            if(ncol(lm_data) == 2){
                              colnames(lm_data) <- c("nest_key", "predictor")
                            }
                            if(ncol(lm_data) == 3){
                              colnames(lm_data) <- c("nest_key", "success", "failure")
                            }
                              
                        
                              boot_re <- bootMer(blup_mod_name, 
                                                 FUN = function(fit){
                                                   randoms <- ranef(fit)$nest_key
                                                   as.numeric(unlist(randoms))
                                                 },
                                                 nsim = boots,
                                                 re.form = ~(temp_s|nest_key)
                              )
                             
                              nest_keys <- rownames(ranef(blup_mod_name)$nest_key)
                              n_keys <- length(nest_keys) 
                              
                              boot_randoms <- lapply(1:nrow(boot_re$t), function(i){
                                sample <- boot_re$t[i, ]
                                
                                data.frame(
                                  nest_key = nest_keys,
                                  random_intercept = sample[1:n_keys],
                                  random_slope = sample[(n_keys + 1):(2*n_keys)]
                                )
                              })
                              
                            if(which_re == "intercept"){
                              save_out <- data.frame(INT = NA, re_int = NA, cil = NA, ciu = NA)
                              for(i in 1:boots){
                                dat <- plyr::join(boot_randoms[[i]], lm_data, "nest_key", "left", "first")
                                tem_mod <- lm(scale(random_intercept) ~ scale(predictor), data = dat)
                                save_out[i, ] <- c(coef(tem_mod), confint(tem_mod)[2, 1:2])
                                cint <- coefficients(summary(tem_mod))
                                if(nrow(cint) > 1){
                                save_out[i, 3:4] <- c(cint[2, 1] - 1.96 * cint[2, 2], cint[2, 1] + 1.96 * cint[2, 2])
                                }
                              }
                            }  
                              
                            if(which_re == "slope"){
                              save_out <- data.frame(INT = NA, re_slp = NA, cil = NA, ciu = NA)
                              for(i in 1:boots){
                                dat <- plyr::join(boot_randoms[[i]], lm_data, "nest_key", "left", "first")
                                tem_mod <- lm(scale(random_slope) ~ scale(predictor), data = dat)
                                save_out[i, ] <- c(coef(tem_mod), confint(tem_mod)[2, 1:2])
                                cint <- coefficients(summary(tem_mod))
                                if(nrow(cint) > 1){
                                save_out[i, 3:4] <- c(cint[2, 1] - 1.96 * cint[2, 2], cint[2, 1] + 1.96 * cint[2, 2])
                                }
                              }
                            }  
                              
                            if(which_re == "both_inter"){
                              save_out <- data.frame(INT = NA, re_intercept = NA, re_slp = NA, interaction = NA,
                                                     in_cil = NA, in_ciu = NA, slp_cil = NA, slp_ciu = NA,
                                                     x_cil = NA, x_ciu = NA)
                              for(i in 1:boots){
                                dat <- plyr::join(boot_randoms[[i]], lm_data, "nest_key", "left", "first")
                                tem_mod <- lm(scale(predictor) ~ scale(random_intercept) * scale(random_slope), data = dat)
                                save_out[i, 1:4] <- coef(tem_mod)
                                cint <- coefficients(summary(tem_mod))
                                if(nrow(cint) > 2){
                                save_out[i, 5:6] <- c(cint[2, 1] - 1.96 * cint[2, 2], cint[2, 1] + 1.96 *cint[2, 2])
                                save_out[i, 7:8] <- c(cint[3, 1] - 1.96 * cint[3, 2], cint[3, 1] + 1.96 *cint[3, 2])
                                save_out[i, 9:10] <- c(cint[4, 1] - 1.96 * cint[4, 2], cint[4, 1] + 1.96 *cint[4, 2])
                                }
                              }
                            }  
                              
                            if(which_re == "both_single"){
                              save_out <- data.frame(INT = NA, re_intercept = NA, re_slp = NA,
                                                     in_cil = NA, in_ciu = NA, slp_cil = NA, slp_ciu = NA)
                              for(i in 1:boots){
                                dat <- plyr::join(boot_randoms[[i]], lm_data, "nest_key", "left", "first")
                                tem_mod <- lm(scale(predictor) ~ scale(random_intercept) + scale(random_slope), data = dat)
                                save_out[i, 1:3] <- coef(tem_mod)
                                cint <- coefficients(summary(tem_mod))
                                if(nrow(cint) > 2){
                                save_out[i, 4:5] <- c(cint[2, 1] - 1.96 * cint[2, 2], cint[2, 1] + 1.96 *cint[2, 2])
                                save_out[i, 6:7] <- c(cint[3, 1] - 1.96 * cint[3, 2], cint[3, 1] + 1.96 *cint[3, 2])
                                }
                              }
                            }  
                              
                            if(which_re == "binom_sin"){
                              save_out <- data.frame(INT = NA, re_intercept = NA, re_slp = NA,
                                                     in_cil = NA, in_ciu = NA, slp_cil = NA, slp_ciu = NA)
                              for(i in 1:boots){
                                dat <- plyr::join(boot_randoms[[i]], lm_data, "nest_key", "left", "first")
                                tem_mod <- glmer(cbind(success, failure) ~ scale(random_intercept) + scale(random_slope) + (1|nest_key), data = dat,
                                                 family = binomial)
                                cint <- coefficients(summary(tem_mod))
                                if(nrow(cint) > 2){
                                  save_out[i, 1:3] <- cint[1:3, 1]
                                  save_out[i, 4:5] <- c(cint[2, 1] - 1.96 * cint[2, 2], cint[2, 1] + 1.96 * cint[2, 2])
                                  save_out[i, 6:7] <- c(cint[3, 1] - 1.96 * cint[3, 2], cint[3, 1] + 1.96 * cint[3, 2])
                                }
                              }
                            }  
                              
                            if(which_re == "binom_x"){
                              save_out <- data.frame(INT = NA, re_intercept = NA, re_slp = NA, interaction = NA,
                                                     in_cil = NA, in_ciu = NA, slp_cil = NA, slp_ciu = NA,
                                                     x_cil = NA, x_ciu = NA)
                              for(i in 1:boots){
                                dat <- na.omit(plyr::join(boot_randoms[[i]], lm_data, "nest_key", "left", "first"))
                                tem_mod <- glmer(cbind(success, failure) ~ scale(random_intercept) * scale(random_slope) + (1|nest_key), data = dat,
                                                 family = binomial)
                                cint <- coefficients(summary(tem_mod))
                                if(nrow(cint) > 2){
                                save_out[i, 1:4] <- cint[1:4, 1]
                                save_out[i, 5:6] <- c(cint[2, 1] - 1.96 * cint[2, 2], cint[2, 1] + 1.96 * cint[2, 2])
                                save_out[i, 7:8] <- c(cint[3, 1] - 1.96 * cint[3, 2], cint[3, 1] + 1.96 * cint[3, 2])
                                save_out[i, 9:10] <- c(cint[4, 1] - 1.96 * cint[4, 2], cint[4, 1] + 1.96 * cint[4, 2])
                                }
                              }
                            }  
                              
                            if(which_re == "binom_tmb_sin"){
                              save_out <- data.frame(INT = NA, re_intercept = NA, re_slp = NA,
                                                     in_cil = NA, in_ciu = NA, slp_cil = NA, slp_ciu = NA)
                              for(i in 1:boots){
                                dat <- plyr::join(boot_randoms[[i]], lm_data, "nest_key", "left", "first")
                                tem_mod <- glmmTMB(cbind(success, failure) ~ scale(random_intercept) + scale(random_slope) + (1|nest_key), 
                                                   data = dat, family = betabinomial(link = "logit"))
                                cint <- coefficients(summary(tem_mod))$cond
                                if(nrow(cint) > 2){
                                save_out[i, 1:3] <- fixef(tem_mod)[[1]][1:3]
                                save_out[i, 4:5] <- c(cint[2, 1] - 1.96 * cint[2, 2], cint[2, 1] + 1.96 * cint[2, 2])
                                save_out[i, 6:7] <- c(cint[3, 1] - 1.96 * cint[3, 2], cint[3, 1] + 1.96 * cint[3, 2])
                                }
                              }
                            }  
                            
                            if(which_re == "binom_tmb_x"){
                              save_out <- data.frame(INT = NA, re_intercept = NA, re_slp = NA, interaction = NA,
                                                     in_cil = NA, in_ciu = NA, slp_cil = NA, slp_ciu = NA,
                                                     x_cil = NA, x_ciu = NA)
                              for(i in 1:boots){
                                dat <- plyr::join(boot_randoms[[i]], lm_data, "nest_key", "left", "first")
                                tem_mod <- glmmTMB(cbind(success, failure) ~ scale(random_intercept) * scale(random_slope) + (1|nest_key),
                                                   data = dat, family = betabinomial(link = "logit"))
                                cint <- coefficients(summary(tem_mod))$cond
                                if(nrow(cint) > 2){
                                save_out[i, 1:4] <- fixef(tem_mod)[[1]][1:4]
                                save_out[i, 5:6] <- c(cint[2, 1] - 1.96 * cint[2, 2], cint[2, 1] + 1.96 * cint[2, 2])
                                save_out[i, 7:8] <- c(cint[3, 1] - 1.96 * cint[3, 2], cint[3, 1] + 1.96 * cint[3, 2])
                                save_out[i, 9:10] <- c(cint[4, 1] - 1.96 * cint[4, 2], cint[4, 1] + 1.96 * cint[4, 2])
                                }
                              }
                            }  
                              
                              return(save_out)
                }
            
# Effect of natural cort variation on incubation and feeding rate ----

  # remove any from above that don't have treatment info                
    cc <- ith_trts[is.na(ith_trts$simple_treatment) == FALSE,]    # this is cort measures in provisioning
    ccm <- ith_trts[is.na(ith_trts_m$simple_treatment) == FALSE,]    # this is cort measures in provisioning for males
    ci <- ith_trts_pre[is.na(ith_trts_pre$simple_treatment) == FALSE, ] # same thing during incubation
                  
  # subset to the columns needed
    cc2 <-cc [,c("nest_key", "bleed1_cort", "bleed2_cort", "bleed3_cort")] # provisioning
    ccm2 <-ccm [,c("nest_key", "bleed1_cort", "bleed2_cort", "bleed3_cort")] # provisioning for males
    ci2 <-ci [,c("nest_key", "bleed1_cort", "bleed2_cort", "bleed3_cort", "simple_treatment")] # incubation
    
  # join cort data to feeding rate data and separately to incubation data
    ef <- plyr::join(exp_feed, cc2, "nest_key", "left", "first")  # provisioning
    efm <- plyr::join(exp_feed, ccm2, "nest_key", "left", "first")  # provisioning for males
    ei <- plyr::join(ab_fon, ci2, "nest_key", "left", "first") # incubation
    
  # save only control nests to look at natural cort variation
    ef_con <- ef[ef$simple_trt == "Control", ] # provisioning
    ef_conm <- efm[efm$simple_trt == "Control", ] # provisioning for males
    ei_con <- ei[ei$simple_treatment == "Control", ] # incubation
  
  # Fit models for provisioning
        # first for base cort
            ef_con_b_dat <- ef_con[ef_con$bleed1_cort < 50 & ef_con$yday - ef_con$j_hatch < 7, ]
            ef_con_b_dat <- ef_con_b_dat %>%
              select(f_feeds, temp_C, bleed1_cort, nest_key)
            ef_con_b_dat <- na.omit(ef_con_b_dat)
            b_q <- quantile(na.omit(ef_con_b_dat$bleed1_cort), c(0.05, 0.95)) # used for plotting at 5th & 95th percentile
            ef_con_b_dat$bleed1_cort_s <- scale(ef_con_b_dat$bleed1_cort)
            mod_nat_b <- gamm(f_feeds ~ t2(temp_C, bleed1_cort), 
                               random = list(nest_key = ~ 1), data = ef_con_b_dat)
            mod_nat_bs <- gamm(f_feeds ~ s(temp_C) + s(bleed1_cort), 
                              random = list(nest_key = ~ 1), data = ef_con_b_dat)
            mod_nat_bs2 <- gamm(f_feeds ~ s(temp_C), 
                              random = list(nest_key = ~ 1), data = ef_con_b_dat)
            AIC(mod_nat_b$lme, mod_nat_bs$lme, mod_nat_bs2$lme) 
                    # delta AIC to temp*cort interaction model = 25.56
                    # n = 8399 hours feeding, 148 nests
          
        # next for stress cort  
            ef_con_s_dat <- ef_con[ef_con$bleed2_cort < 120 & ef_con$bleed1_cort < 50 &
                                   ef_con$yday - ef_con$j_hatch < 7, ]
            ef_con_s_dat <- ef_con_s_dat %>% select(f_feeds, temp_C, bleed1_cort, bleed2_cort, nest_key)
            ef_con_s_dat <- na.omit(ef_con_s_dat)
            ef_con_s_dat$bleed1_cort_s <- as.numeric(scale(ef_con_s_dat$bleed1_cort))
            b_s <- quantile(na.omit(ef_con_b_dat$bleed3_cort), c(0.05, 0.95))
            mu_b <- mean(na.omit(ef_con_b_dat$bleed1_cort))
            mod_nat_s2 <- gamm(f_feeds ~ t2(temp_C, bleed2_cort, bleed1_cort),
                               random = list(nest_key = ~ 1), data = ef_con_s_dat)
            mod_nat_s1 <- gamm(f_feeds ~ t2(temp_C, bleed1_cort),
                               random = list(nest_key = ~ 1), data = ef_con_s_dat)
            mod_nat_s1b <- gamm(f_feeds ~ t2(temp_C, bleed2_cort),
                               random = list(nest_key = ~ 1), data = ef_con_s_dat)
            mod_nat_s <- gamm(f_feeds ~ t2(temp_C, bleed2_cort) + bleed1_cort_s, 
                               random = list(nest_key = ~ 1), data = ef_con_s_dat)
            mod_nat_sd <- gamm(f_feeds ~ s(temp_C) + s(bleed2_cort), 
                               random = list(nest_key = ~1), data = ef_con_s_dat)
            mod_nat_sd2 <- gamm(f_feeds ~ s(temp_C), 
                               random = list(nest_key = ~1), data = ef_con_s_dat)
            AIC(mod_nat_sd$lme, mod_nat_sd2$lme, mod_nat_s1b$lme)
            #vis.gam(mod_nat_s$gam, plot.type = "contour", nlevels = 20, type = "response", too.far = .2)
                # There is no support for the stress cort * temperature interaction with feeding rate
          
        # make predictions at specific values from the models  
            pred_mnb <- ggpredict(mod_nat_b, terms = c("temp_C [all]", "bleed1_cort [.48, 11.2]"))
            #pred_mns <- ggpredict(mod_nat_s1, terms = c("temp_C [all]", "bleed2_cort [8, 50]", "bleed3_cort [5.5, 31]"))
            #pred_mns$group_level <- paste(pred_mns$group, pred_mns$facet, sep = "_")
            
            gp1 <- ggplot(pred_mnb, aes(x = x, y = predicted, color = group, fill = group)) +
              geom_line() +
              geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = 0.6, color = "transparent") +
              #ggtitle("Base cort") +
              theme_classic() +
              theme(panel.grid = element_blank(), axis.text = element_text(size = 12, family = "Helvetica Light"),
                    axis.title = element_text(size = 14, family = "Helvetica Light"), 
                    title = element_text(family = "Helvetica Light")) +
              labs(x = "Temperature", y = "Hourly feeding rate") +
              scale_color_viridis(discrete = TRUE, option = "H", begin = 0.2, end = 0.7) +
              scale_fill_viridis(discrete = TRUE, option = "H", begin = 0.2, end = 0.7) +
              theme(legend.title = element_blank()) +
              scale_y_continuous(breaks = seq(1, 15, 3))
            
            ggsave(here::here("output_plots", "base_cort_feeding.svg"), gp1, device = "svg", width = 5, height = 3.6,
                   units = "in", dpi = 300)
            
            # gp1b <- ggplot(pred_mns, aes(x = x, y = predicted, color = group, fill = group)) +
            #   geom_line() +
            #   geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = 0.5) +
            #   ggtitle("Stress cort")
            # 
            # gp2 <- ggplot(pred_mns, aes(x = x, y = predicted, color = group, fill = group)) +
            #   geom_line() +
            #   geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = 0.5, color = "transparent") +
            #   facet_wrap(~ facet, labeller = labeller(facet = c("5.5" = "Weak negative feedback", "31" = "Strong negative feedback"))) +
            #   theme(legend.title = element_text("Stress cort"))
              
    # same as above for on and off bout durations
          # first for on bouts
            ei_con$stage <- ei_con$yday - ei_con$hatch_doy
            ei_con_dat <- ei_con[ei_con$type == "Off" &
                                   ei_con$stage > -6 & ei_con$stage < 0 &
                                   ei_con$duration_m < 50, ]
            ei_con_dat <- ei_con_dat %>%
              select(duration_m, temp_C, bleed1_cort, nest_key)
            ei_con_dat <- na.omit(ei_con_dat)
            i_q <- quantile(na.omit(ei_con_dat$bleed1_cort), c(0.05, 0.95)) # used for plotting at 5th & 95th percentile
            mod_nat_ion <- gamm(duration_m ~ t2(temp_C, bleed1_cort), 
                               random = list(nest_key = ~ 1), data = ei_con_dat)
            mod_nat_ion2 <- gamm(duration_m ~ s(temp_C) + s(bleed1_cort), 
                              random = list(nest_key = ~ 1), data = ei_con_dat)
            mod_nat_ion3 <- gamm(duration_m ~ s(temp_C), 
                              random = list(nest_key = ~ 1), data = ei_con_dat)
            AIC(mod_nat_ion$lme, mod_nat_ion2$lme, mod_nat_ion3$lme) 
                    # no evidence for base cort effect on onbout duration
            
            ei_con$stage <- ei_con$yday - ei_con$hatch_doy
            ei_con_dat <- ei_con[ei_con$type == "On" &
                                   ei_con$stage > -6 & ei_con$stage < 0 &
                                   ei_con$duration_m < 50, ]
            ei_con_dat2 <- ei_con_dat %>%
              select(duration_m, temp_C, bleed1_cort, bleed2_cort, nest_key)
            ei_con_dat2 <- na.omit(ei_con_dat2)
            i_q2 <- quantile(na.omit(ei_con_dat2$bleed2_cort), c(0.05, 0.95)) # used for plotting at 5th & 95th percentile
            mod_nat_ions <- gamm(duration_m ~ t2(temp_C, bleed2_cort), 
                                random = list(nest_key = ~ 1), data = ei_con_dat2)
            mod_nat_ion2s <- gamm(duration_m ~ s(temp_C) + s(bleed2_cort), 
                                 random = list(nest_key = ~ 1), data = ei_con_dat2)
            mod_nat_ion3s <- gamm(duration_m ~ s(temp_C), 
                                 random = list(nest_key = ~ 1), data = ei_con_dat2)
            # mod_nat_ion4s <- gamm(duration_m ~ t2(temp_C, bleed2_cort, bleed1_cort), 
            #                      random = list(nest_key = ~ 1), data = ei_con_dat2)
            AIC(mod_nat_ions$lme, mod_nat_ion2s$lme, mod_nat_ion3s$lme) 
            # no evidence for base cort effect on onbout duration
            
                      

# Early life effects ----
    # sample size is small in most cases
            
    # also tried looking at mother-daughter regression...small sample size but again nothing
            
            
    # list of all nestlings
        ith_born <- ith_capture %>%
          filter(adult_or_nestling == "Nestling")
        ith_born <- ith_born[, c("band", "nest_key", "exp_year", "age",
                                 "mass", "headbill", "flatwing")]
        ith_born <- plyr::join(ith_born, nests[, c("nest_key", "exp_year", "clutch_comp_doy", "clutch_size", "hatch_doy", "nest_fate_doy",
                                                   "female_id", "male_id")], "nest_key", "left", "first")
        colnames(ith_born)[1] <- "natal_nest_key"
        colnames(ith_born)[3:ncol(ith_born)] <- paste("natal", colnames(ith_born)[3:ncol(ith_born)], sep = "_")

        # moms <- unique(na.omit(ad_feed$natal_female_id))
        # moms2 <- ad_feed[ad_feed$band %in% moms, ]
        # moms2 <- moms2[, c("band", "intercept", "slope")]
        # colnames(moms2) <- c("natal_female_id", "mom_intercept", "mom_slope")
        # ad_ith2 <- plyr::join(ad_ith, moms2, "natal_female_id", "left", "first")

    # adults in behavior data feeding
        ad_feedf <- as.data.frame(ranef(mffeed)$nest_key)
        colnames(ad_feedf) <- c("intercept", "slope")
        ad_feedf$nest_key <- rownames(ad_feedf)
        ad_feedf <- plyr::join(ad_feedf, nests, "nest_key", "left", "first")
        ad_feedf <- ad_feedf[, c("nest_key", "intercept", "slope", "exp_year", "female_id")]
        ad_feedf$band <- ad_feedf$female_id
        ad_feedf <- plyr::join(ad_feedf, ith_born, "band", "left", "first")
        ad_feedf$born_in_ithaca <- "yes"
        for(i in 1:nrow(ad_feedf)){
          if(is.na(ad_feedf$natal_nest_key[i])==TRUE){ad_feedf$born_in_ithaca[i] <- "no"}
        }
        
      # adults in behavior data feeding males
        ad_feedm <- as.data.frame(ranef(mmfeed)$nest_key)
        colnames(ad_feedm) <- c("intercept", "slope")
        ad_feedm$nest_key <- rownames(ad_feedm)
        ad_feedm <- plyr::join(ad_feedm, nests, "nest_key", "left", "first")
        ad_feedm <- ad_feedm[, c("nest_key", "intercept", "slope", "exp_year", "male_id")]
        ad_feedm$band <- ad_feedm$male_id
        ad_feedm <- plyr::join(ad_feedm, ith_born, "band", "left", "first")
        ad_feedm$born_in_ithaca <- "yes"
        for(i in 1:nrow(ad_feedm)){
          if(is.na(ad_feedm$natal_nest_key[i])==TRUE){ad_feedm$born_in_ithaca[i] <- "no"}
        }        
        
        ad_feed <- bind_rows(ad_feedf, ad_feedm)
        ad_feed <- ad_feed[, !duplicated(names(ad_feed))]
        
    # adults in behavior data incubation on bouts
        ad_inc_on <- as.data.frame(ranef(mi_on)$nest_key)
        colnames(ad_inc_on) <- c("intercept", "slope")
        ad_inc_on$nest_key <- rownames(ad_inc_on)
        ad_inc_on <- plyr::join(ad_inc_on, nests, "nest_key", "left", "first")
        ad_inc_on <- ad_inc_on[, c("nest_key", "intercept", "slope", "exp_year", "female_id")]
        ad_inc_on$band <- ad_inc_on$female_id
        ad_inc_on <- plyr::join(ad_inc_on, ith_born, "band", "left", "first")
        ad_inc_on$born_in_ithaca <- "yes"
        for(i in 1:nrow(ad_inc_on)){
          if(is.na(ad_inc_on$natal_nest_key[i])==TRUE){ad_inc_on$born_in_ithaca[i] <- "no"}
        }    
        ad_inc_on <- ad_inc_on[ad_inc_on$born_in_ithaca == "yes", !duplicated(names(ad_inc_on))]
        
    # adults in behavior data incubation off bouts
        ad_inc_off <- as.data.frame(ranef(mi_off)$nest_key)
        colnames(ad_inc_off) <- c("intercept", "slope")
        ad_inc_off$nest_key <- rownames(ad_inc_off)
        ad_inc_off <- plyr::join(ad_inc_off, nests, "nest_key", "left", "first")
        ad_inc_off <- ad_inc_off[, c("nest_key", "intercept", "slope", "exp_year", "female_id")]
        ad_inc_off$band <- ad_inc_off$female_id
        ad_inc_off <- plyr::join(ad_inc_off, ith_born, "band", "left", "first")
        ad_inc_off$born_in_ithaca <- "yes"
        for(i in 1:nrow(ad_inc_off)){
          if(is.na(ad_inc_off$natal_nest_key[i])==TRUE){ad_inc_off$born_in_ithaca[i] <- "no"}
        }    
        ad_inc_off <- ad_inc_off[ad_inc_off$born_in_ithaca == "yes", !duplicated(names(ad_inc_off))]

    # add temperature data
        day_high <-
          w_hour %>%
          dplyr::group_by(year, yday) %>%
          summarise(max = max(na.omit(temp_C)))

        ad_ith <- ad_feed[ad_feed$born_in_ithaca == "yes", ]
        for(i in 1:nrow(ad_ith)){
          sub <- subset(w_hour, w_hour$year == ad_ith$natal_exp_year[i] &
                          w_hour$hour > 5 & w_hour$hour < 21 &
                          w_hour$yday > ad_ith$natal_hatch_doy[i] - 25 &
                          w_hour$yday < ad_ith$natal_hatch_doy[i] + 25)
          ad_ith$early_inc_T[i] <- mean(na.omit(subset(sub$temp_C, sub$yday < (ad_ith$natal_hatch_doy[i] - 7) &
                                                  sub$yday > (ad_ith$natal_hatch_doy[i] - 14))))
          ad_ith$late_inc_T[i] <- mean(na.omit(subset(sub$temp_C, sub$yday < ad_ith$natal_hatch_doy[i] &
                                                 sub$yday > (ad_ith$natal_hatch_doy[i] - 6))))
          ad_ith$early_prov_T[i] <- mean(na.omit(subset(sub$temp_C, sub$yday > (ad_ith$natal_hatch_doy[i]) &
                                                   sub$yday < (ad_ith$natal_hatch_doy[i] + 11))))
          ad_ith$late_prov_T[i] <- mean(na.omit(subset(sub$temp_C, sub$yday > (ad_ith$natal_hatch_doy[i] + 10) &
                                                  sub$yday < (ad_ith$natal_hatch_doy[i] + 20))))

          sub2i <- subset(day_high, day_high$year == ad_ith$natal_exp_year[i] &
                            day_high$yday < ad_ith$natal_hatch_doy[i] &
                            day_high$yday > ad_ith$natal_hatch_doy[i] - 14)
          sub2p <- subset(day_high, day_high$year == ad_ith$natal_exp_year[i] &
                            day_high$yday > ad_ith$natal_hatch_doy[i] - 1 &
                            day_high$yday < ad_ith$natal_hatch_doy[i] + 14)

          ad_ith$nd_18_inc[i] <- nrow(subset(sub2i, sub2i$max < 18.5))
          ad_ith$nd_18_prov[i] <- nrow(subset(sub2p, sub2p$max < 18.5))
        }

        ad_ith$prov_T <- (ad_ith$early_prov_T + ad_ith$late_prov_T) / 2
        ad_ith$inc_T <- (ad_ith$early_inc_T + ad_ith$late_inc_T) / 2

        ad_ith$ein_s <- scale(ad_ith$early_inc_T)
        ad_ith$lin_s <- scale(ad_ith$late_inc_T)
        ad_ith$epr_s <- scale(ad_ith$early_prov_T)
        ad_ith$lpr_s <- scale(ad_ith$late_prov_T)

        # mmm <- lm(intercept ~ ein_s + lin_s + epr_s + lpr_s + scale(slope) , data = ad_ith)
        # mmm2 <- lm(slope ~ ein_s + lin_s + epr_s + lpr_s + scale(intercept), data = ad_ith)
        # 
        # tab_model(mmm, mmm2)
        
        
      ## Process from the data frame built above
        ad_ith$sex <- "female"
        ad_ith$sex[ad_ith$band == ad_ith$male_id] <- "male"
        ad_ith2 <- ad_ith[, c("nest_key", "band", "sex", "ein_s", "lin_s", "epr_s", "lpr_s")]
        
        # join cort data to feeding rate data and separately to incubation data
          early_life <- plyr::join(exp_feed, ad_ith2, "nest_key", "left", "all")   
          early_life <- early_life %>%
            filter(is.na(lin_s) == FALSE)
        
        # save only control nests to look at natural cort variation
          early_life <- early_life[early_life$simple_trt == "Control",]
          
          early_f_dat <- early_life %>% filter(sex == "female")
          early_f_dat$off_s <- scale(early_f_dat$Offset)
          early_f_dat$hour_s <- scale(early_f_dat$hour)
          el_mod <- gamm(f_feeds ~ t2(temp_C, lin_s),
                         random = list(nest_key = ~ 1), data = early_f_dat)
          early_q <- quantile(na.omit(early_f_dat$lin_s), c(0.05, 0.95))
          el_mod2 <- gamm(f_feeds ~ t2(temp_C),
                         random = list(nest_key = ~ 1), data = early_f_dat)
          el_mod3 <- gamm(f_feeds ~ s(temp_C) + s(lin_s),
                          random = list(nest_key = ~1), data = early_f_dat)
          AIC(el_mod$lme, el_mod2$lme, el_mod3$lme)
              # delta aic for temp * inc_temp interaction = 19.7
          
          pred_early <- ggpredict(el_mod, terms = c("temp_C [all]", "lin_s [-1.754, 1.462]" ))
          
          early_plot <- ggplot(pred_early, aes(x = x, y = predicted, color = group, fill = group)) +
            geom_line() +
            geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = 0.5, color = "transparent") +
            #ggtitle("Early life temperature") +
            xlim(c(8, 33)) +
            theme_classic() +
            theme(panel.grid = element_blank(), axis.text = element_text(size = 12, family = "Helvetica Light"),
                  axis.title = element_text(size = 14, family = "Helvetica Light"), 
                  title = element_text(family = "Helvetica Light")) +
            scale_fill_manual(values = c(off_color, on_color)) +
            scale_color_manual(values = c(off_color, on_color)) +
            labs(y = "Feeding Trips (hourly)", x = "Ambient Temperature (°C)")
          early_plot
          
          ggsave(here::here("output_plots", "early_life_feeding.svg"), early_plot, device = "svg", width = 5, height = 3.6,
                 units = "in", dpi = 300)
          
          
        # this was an earlier iteration with some sliding windows and BLUPs. Not included  
            # ## join to incubation
            #     
            #     early_on <- ab_fon %>% filter(type == "On") %>% 
            #       plyr::join(ad_ith2, "nest_key", "left", "all") %>%
            #       filter(sex == "female") %>%
            #       filter(is.na(lin_s) == FALSE)
            #     
            #     early_off <- ab_fon %>% filter(type == "Off") %>% 
            #       plyr::join(ad_ith2, "nest_key", "left", "all") %>%
            #       filter(sex == "female") %>%
            #       filter(is.na(lin_s) == FALSE)
            #   
            #   
            # # bring in other blups
            #   
            #   ad_ith_c <- plyr::join(ad_ith, cs3_comp, "nest_key", "left", "first")
            #   ad_ith_c <- ad_ith_c[, !duplicated(names(ad_ith_c))]
            #   ad_ith <- ad_ith_c
            #   
            # # make a loop for sliding window analysis
            #   sliding_temp <- data.frame(
            #     start = seq(-13, 18, 1),
            #     end = seq(-9, 22, 1),
            #     est_int = NA,
            #     low_int = NA,
            #     hi_int = NA,
            #     est_slp = NA,
            #     low_slp = NA,
            #     low_int = NA
            #   )
            #   
            #   ad_feedf$sex <- "female"
            #   ad_feedm$sex <- "male"
            #   ad_ith <- bind_rows(ad_feedf, ad_feedm)
            #   #ad_ith <- ad_feedf
            #   ad_ith <- subset(ad_ith, ad_ith$born_in_ithaca == "yes")
            #   for(k in 1:32){
            #     for(i in 1:nrow(ad_ith)){
            #       sub <- subset(w_hour, w_hour$year == ad_ith$natal_exp_year[i] &
            #                       w_hour$hour > 5 & w_hour$hour < 21 &
            #                       w_hour$yday > ad_ith$natal_hatch_doy[i] - 25 &
            #                       w_hour$yday < ad_ith$natal_hatch_doy[i] + 25)
            #       ad_ith$natal_T[i] <- mean(na.omit(subset(sub$temp_C, sub$yday > (ad_ith$natal_hatch_doy[i] - 13 + k) &
            #                                                  sub$yday < (ad_ith$natal_hatch_doy[i] - 6 + k))))
            #       ad_ith$natal_Ts <- scale(ad_ith$natal_T)
            #     }
            #     temp_int <- lm(intercept ~ natal_Ts + scale(slope), data = ad_ith)
            #     temp_slp <- lm(slope ~ natal_Ts + scale(intercept), data = ad_ith)
            #     tic <- coefficients(summary(temp_int))
            #     tis <- coefficients(summary(temp_slp))
            #     
            #     sliding_temp$est_int[k] <- tic[2, 1]
            #     sliding_temp$low_int[k] <- tic[2, 1] - 1.96*tic[2, 2]
            #     sliding_temp$hi_int[k] <- tic[2, 1] + 1.96*tic[2, 2]
            #     
            #     sliding_temp$est_slp[k] <- tis[2, 1]
            #     sliding_temp$low_slp[k] <- tis[2, 1] - 1.96*tis[2, 2]
            #     sliding_temp$hi_slp[k] <- tis[2, 1] + 1.96*tis[2, 2]
            #     
            #     print(k)
            #   }
            #   
            #   
            #   sliding_temp$int_sig <- "white"
            #   sliding_temp$slp_sig <- "white"
            #   for(i in 1:nrow(sliding_temp)){
            #     if(sliding_temp$low_int[i] > 0){sliding_temp$int_sig[i] <- "orange"}
            #     if(sliding_temp$hi_int[i] < 0){sliding_temp$int_sig[i] <- "orange"}
            #     if(sliding_temp$low_slp[i] > 0){sliding_temp$slp_sig[i] <- "slateblue"}
            #     if(sliding_temp$hi_slp[i] < 0){sliding_temp$slp_sig[i] <- "slateblue"}
            #   }
            #   
            #   ggplot(sliding_temp) +
            #     geom_hline(yintercept = 0) +
            #     geom_vline(xintercept = 0, linetype = "dashed", color = "gray50") +
            #     geom_segment(aes(x = start + 0.1, xend = start + 0.1, y = low_int, yend = hi_int), color = "orange") +
            #     geom_point(aes(x = start + 0.1, y = est_int), fill = sliding_temp$int_sig, shape = 21, size = 2) +
            #     geom_segment(aes(x = start - 0.1, xend = start - 0.1, y = low_slp, yend = hi_slp), color = "slateblue") +
            #     geom_point(aes(x = start - 0.1, y = est_slp), fill = sliding_temp$slp_sig, shape = 21, size = 2) +
            #     theme_classic() +
            #     #ylim(-1.6, 1.6) +
            #     theme(panel.grid = element_blank(), axis.text = element_text(size = 12), axis.title = element_text(size = 14)) +
            #     labs(x = "Start of 5-Day Window \n Relative to Hatch", y = "Standardized Effect Size")
        
# Individual random effect estimates vs. fitness metrics ----
    # from the simple linear model BLUPs in figure 2
        # female on bout random effects
          reon <- ranef(mi_on)$nest_key
          colnames(reon) <- c("on_intercept", "on_slope")
          reon$nest_key <- rownames(reon)
          
        # female off bout random effects
          reoff <- ranef(mi_off)$nest_key
          colnames(reoff) <- c("off_intercept", "off_slope")
          reoff$nest_key <- rownames(reoff)
        
        # female feeding random effects
          ref <- ranef(mffeed)$nest_key
          colnames(ref) <- c("fem_feed_intercept", "fem_feed_slope")
          ref$nest_key <- rownames(ref)
          
        # male feeding random effects
          rem <- ranef(mmfeed)$nest_key
          colnames(rem) <- c("mal_feed_intercept", "mal_feed_slope")
          rem$nest_key <- rownames(rem)
        
        # trim nest records down to necessary columns
            nest_red <- nests[, c("nest_key", "exp_year", "attempt_num", "clutch_init_doy",
                                  "clutch_comp_doy", "clutch_size", "hatch_doy", "brood_size_hatching", 
                                  "day6_brood_mass", "nest_fate_doy", "simple_fate", "num_fledged")]
        
        # join to nest records
              reon <- plyr::join(reon, nest_red, "nest_key", "left", "first")
              reoff <- plyr::join(reoff, nest_red, "nest_key", "left", "first")
              ref <- plyr::join(ref, nest_red, "nest_key", "left", "first")
              rem <- plyr::join(rem, nest_red, "nest_key", "left", "first")
              
        # add in number of nestlings that fledged and then recruited to population next year
              recruit <- ith_capture %>%
                filter(adult_or_nestling == "Nestling", age == 12, nestling_fate == "Fledged", exp_year > 2012, location == "Ithaca") %>%
                dplyr::select(nest_key, exp_year, band, age, adult_or_nestling, nestling_fate)
              
              recruit2 <- ith_capture %>%
                filter(adult_or_nestling == "Adult", exp_year > 2013, location == "Ithaca") %>%
                dplyr::select(exp_year, band) %>%
                filter(!duplicated(band))
              
              colnames(recruit2)[1] <- "adult_year"
              
              recruit3 <- plyr::join(recruit, recruit2, "band", "left", "first")
              recruit3 <- recruit3[is.na(recruit3$adult_year) == FALSE, ]        
          
        # add in number of nestlings alive at banding (age 12)
              # for on bouts
                  for(i in 1:nrow(reon)){
                    sub <- subset(ith_capture, ith_capture$nest_key == reon$nest_key[i] &
                                    ith_capture$adult_or_nestling == "Nestling" &
                                    ith_capture$age == 12)
                    reon$num_d12[i] <- length(unique(sub$band))
                    reon$d12_av_mass[i] <- mean(na.omit(sub$mass))
                    reon$d12_av_fw[i] <- mean(na.omit(sub$flatwing)) 
                    reon$d12_av_hb[i] <- mean(na.omit(sub$headbill))
                    
                    sub2 <- subset(recruit3, recruit3$nest_key == reon$nest_key[i])
                    reon$recruits[i] <- nrow(sub2)
                  }
              
              # for off bouts
                  for(i in 1:nrow(reoff)){
                    sub <- subset(ith_capture, ith_capture$nest_key == reoff$nest_key[i] &
                                    ith_capture$adult_or_nestling == "Nestling" &
                                    ith_capture$age == 12)
                    reoff$num_d12[i] <- length(unique(sub$band))
                    reoff$d12_av_mass[i] <- mean(na.omit(sub$mass))
                    reoff$d12_av_fw[i] <- mean(na.omit(sub$flatwing)) 
                    reoff$d12_av_hb[i] <- mean(na.omit(sub$headbill))
                    
                    sub2 <- subset(recruit3, recruit3$nest_key == reoff$nest_key[i])
                    reoff$recruits[i] <- nrow(sub2)
                  }
              
              # for female feeding
                  for(i in 1:nrow(ref)){
                    sub <- subset(ith_capture, ith_capture$nest_key == ref$nest_key[i] &
                                    ith_capture$adult_or_nestling == "Nestling" &
                                    ith_capture$age == 12)
                    ref$num_d12[i] <- length(unique(sub$band))
                    ref$d12_av_mass[i] <- mean(na.omit(sub$mass))
                    ref$d12_av_fw[i] <- mean(na.omit(sub$flatwing)) 
                    ref$d12_av_hb[i] <- mean(na.omit(sub$headbill))
                    
                    sub2 <- subset(recruit3, recruit3$nest_key == ref$nest_key[i])
                    ref$recruits[i] <- nrow(sub2)
                  }
              
              # for male feeding    
                  for(i in 1:nrow(rem)){
                    sub <- subset(ith_capture, ith_capture$nest_key == rem$nest_key[i] &
                                    ith_capture$adult_or_nestling == "Nestling" &
                                    ith_capture$age == 12)
                    rem$num_d12[i] <- length(unique(sub$band))
                    rem$d12_av_mass[i] <- mean(na.omit(sub$mass))
                    rem$d12_av_fw[i] <- mean(na.omit(sub$flatwing)) 
                    rem$d12_av_hb[i] <- mean(na.omit(sub$headbill))
                    
                    sub2 <- subset(recruit3, recruit3$nest_key == rem$nest_key[i])
                    rem$recruits[i] <- nrow(sub2)
                  }

        # combine male/female feed and on/off bout length 
              reon$type <- "On"
              reoff$type <- "Off"
              ref$sex <- "Female"
              rem$sex <- "Male"
              
              reon2 <- reon
              colnames(reon2)[2:3] <- c("intercept", "slope")
              reoff2 <- reoff
              colnames(reoff2)[2:3] <- c("intercept", "slope")
              ref2 <- ref 
              colnames(ref2)[2:3] <- c("intercept", "slope")
              rem2 <- rem
              colnames(rem2)[2:3] <- c("intercept", "slope")
              
              re_inc <- rbind(reon2, reoff2)
              re_feed <- rbind(ref2, rem2)
          
              
        # plots and basic models
              
              # on and off bouts
                    inc_slope_p <- re_inc %>%
                      filter(is.na(num_fledged) == FALSE, simple_fate == "failed" | simple_fate == "fledged") %>%
                      ggplot(mapping = aes(x = as.factor(simple_fate), y = scale(slope), fill = type, color = type)) + 
                      geom_boxplot(alpha = 0.4, outlier.shape = NA) +
                      scale_fill_manual(values = c("slateblue", "coral3")) +
                      scale_color_manual(values = c("slateblue", "coral3")) +
                      theme_classic() +
                      guides(fill = "none", color = "none") +
                      theme(panel.grid = element_blank(), axis.text = element_text(size = 12)) +
                      coord_cartesian(ylim = c(-2.4, 2.4)) 
                      #scale_y_continuous(breaks = seq(-10, 10, 5))
                    
                    inc_int_p <- re_inc %>%
                      filter(is.na(num_fledged) == FALSE, simple_fate == "failed" | simple_fate == "fledged") %>%
                      ggplot(mapping = aes(x = as.factor(simple_fate), y = scale(intercept), fill = type, color = type)) + 
                      geom_boxplot(alpha = 0.4, outlier.shape = NA) +
                      scale_fill_manual(values = c("slateblue", "coral3")) +
                      scale_color_manual(values = c("slateblue", "coral3")) +
                      theme_classic() +
                      guides(fill = "none", color = "none") +
                      theme(panel.grid = element_blank(), axis.text = element_text(size = 12)) +
                      coord_cartesian(ylim = c(-2.4, 2.4)) 
                    #scale_y_continuous(breaks = seq(-10, 10, 5))
                    
                    # basic models
                        reon_sf <- reon[reon$simple_fate == "failed" | reon$simple_fate == "fledged", ]
                        reon_sf$sf_bin <- ifelse(reon_sf$simple_fate == "fledged", 1, 0)
                        msimp_on <- glm(sf_bin ~ on_intercept * on_slope, data = reon_sf, family = binomial)
                        msimp_onr <- glm(sf_bin ~ on_intercept + on_slope, data = reon_sf, family = binomial)
                        
                        reoff_sf <- reoff[reoff$simple_fate == "failed" | reoff$simple_fate == "fledged", ]
                        reoff_sf$sf_bin <- ifelse(reoff_sf$simple_fate == "fledged", 1, 0)
                        msimp_off <- glm(sf_bin ~ off_intercept * off_slope, data = reoff_sf, family = binomial)
                        msimp_offr <- glm(sf_bin ~ off_intercept + off_slope, data = reoff_sf, family = binomial)
                    
              
            # male and female feeding
                        
                        feed_int_p <- re_feed %>%
                          filter(is.na(num_fledged) == FALSE,  brood_size_hatching > 0, simple_fate == "failed" | simple_fate == "fledged") %>%
                          ggplot(mapping = aes(x = as.factor(simple_fate), y = scale(intercept), fill = sex, color = sex)) + 
                          geom_boxplot(alpha = 0.4, outlier.shape = NA) +
                          scale_fill_manual(values = c(fem_color, mal_color)) +
                          scale_color_manual(values = c(fem_color, mal_color)) +
                          #coord_cartesian(ylim = c(-7, 7)) +
                          #scale_y_continuous(breaks = seq(-10, 10, 5)) +
                          theme_classic() +
                          theme(panel.grid = element_blank(), axis.text = element_text(size = 12))
                        
                        feed_slp_p <- re_feed %>%
                          filter(is.na(num_fledged) == FALSE,  brood_size_hatching > 0, simple_fate == "failed" | simple_fate == "fledged") %>%
                          ggplot(mapping = aes(x = as.factor(simple_fate), y = scale(slope), fill = sex, color = sex)) + 
                          geom_boxplot(alpha = 0.4, outlier.shape = NA) +
                          scale_fill_manual(values = c(fem_color, mal_color)) +
                          scale_color_manual(values = c(fem_color, mal_color)) +
                          #coord_cartesian(ylim = c(-7, 7)) +
                          #scale_y_continuous(breaks = seq(-10, 10, 5)) +
                          theme_classic() +
                          theme(panel.grid = element_blank(), axis.text = element_text(size = 12))
                        
                      # basic models
                        ref_sf <- ref[ref$simple_fate == "failed" | ref$simple_fate == "fledged", ]
                        #ref_sf <- ref_sf %>% filter(num_d12 > 0) # results hold even if limiting to nests alive on d12
                        ref_sf$sf_bin <- ifelse(ref_sf$simple_fate == "fledged", 1, 0)
                        msimp_ff <- glm(sf_bin ~ fem_feed_intercept * fem_feed_slope, data = ref_sf, family = binomial)
                        msimp_ffr <- glm(sf_bin ~ fem_feed_intercept + fem_feed_slope, data = ref_sf, family = binomial)
                        
                        rem_sf <- rem[rem$simple_fate == "failed" | rem$simple_fate == "fledged", ]
                        rem_sf$sf_bin <- ifelse(rem_sf$simple_fate == "fledged", 1, 0)
                        msimp_mf <- glm(sf_bin ~ mal_feed_intercept * mal_feed_slope, data = rem_sf, family = binomial)
                        msimp_mfr <- glm(sf_bin ~ mal_feed_intercept + mal_feed_slope, data = rem_sf, family = binomial) 
                        
                        
                        
                        
                        
                        
                        
                        
                        
              re_inc %>%
                filter(is.na(num_fledged) == FALSE) %>%
                ggplot(mapping = aes(x = as.factor(brood_size_hatching), y = intercept, fill = type, color = type)) + 
                geom_boxplot(alpha = 0.4, outlier.shape = NA) +
                scale_fill_manual(values = c("slateblue", "coral3")) +
                scale_color_manual(values = c("slateblue", "coral3")) +
                theme_classic() +
                theme(panel.grid = element_blank(), axis.text = element_text(size = 12)) #+
              #coord_cartesian(ylim = c(-7, 7)) +
              #scale_y_continuous(breaks = seq(-10, 10, 5))
              
              re_inc %>%
                filter(is.na(num_fledged) == FALSE, is.na(brood_size_hatching) == FALSE) %>%
                ggplot(mapping = aes(x = as.numeric(d12_av_hb), y = intercept, fill = type, color = type)) +
                geom_smooth(method = "lm")
              
              
              

              
              re_feed %>%
                filter(is.na(num_fledged) == FALSE, num_fledged > 0, brood_size_hatching > 2) %>%
              ggplot(mapping = aes(x = as.factor(num_fledged), y = slope, fill = sex, color = sex)) + 
                geom_boxplot(alpha = 0.4, outlier.shape = NA) +
                scale_fill_manual(values = c(fem_color, mal_color)) +
                scale_color_manual(values = c(fem_color, mal_color)) +
                theme_classic() +
                theme(panel.grid = element_blank(), axis.text = element_text(size = 12)) +
                coord_cartesian(ylim = c(-3, 3)) +
                scale_y_continuous(breaks = seq(-10, 10, 2))
              
              re_feed %>%
                filter(num_d12 > 0, is.na(num_fledged) == FALSE) %>%
                ggplot(mapping = aes(x = d12_av_fw, y = slope, fill = sex, color = sex)) + 
                geom_jitter(width = 0.2, size = 0.4) +
                geom_smooth(method = "lm")
              
              summary(lm(num_fledged ~ slope*intercept, data = re_feed))
              summary(lm(d12_av_mass ~ slope*intercept + num_d12, data = re_feed[re_feed$sex == "Male", ]))
              
              ggplot(ref, mapping = aes(x = fem_feed_intercept, y = d12_av_hb)) +
                #geom_point() +
                geom_smooth(method = "lm") #+
               # facet_wrap(~as.factor(num_d12), scales = "free")
              
              nestling_fem <- ith_capture[ith_capture$nest_key %in% ref$nest_key, ]
              nestling_fem <- nestling_fem %>%
                filter(adult_or_nestling == "Nestling", age == 12)
              
              
        # Fit models asking if intercept and slope from random effects predict fitness outcomes
            # for each model, I'll fit intercept*slope and remove the interaction if not supported
            # for each effect size, I'll need to generate a confidence interval by bootsrapping
            # from the BLUPs and put together a full table with estimates and CIs
            # there will be a total of 22 models labled m1 to m22 (with reduced models by adding r)
              
              boot.set <- 100 # how many bootstrap iterations to do for blups
              boot.set2 <- 10
              
              # combine data frame. If adding any more columns put them at end or it will break code assignment below
                  est_df <- data.frame(model = paste("m", seq(1, 26, 1), sep = ""),
                                   response = c(rep("hatch", 2), rep("survd12", 4), rep("fledge", 4),
                                                rep("avmass", 4), rep("avhead", 4), rep("avwing", 4),
                                                rep("recruit", 4)),
                                   type = c(rep("binom", 2), rep("binomT", 4), rep("binom", 2), rep("binomT", 2), rep("norm", 12),
                                            rep("binom", 4)),
                                   reduced = rep("yes", 26),
                                   int_est = NA, slp_est = NA, x_est = NA,
                                   int_low = NA, int_hi = NA,
                                   slp_low = NA, slp_hi = NA,
                                   x_low = NA, x_hi = NA,
                                   raw_int_est = NA, raw_slp_est = NA, raw_x_est = NA,
                                   raw_int_low = NA, raw_int_hi = NA,
                                   raw_slp_low = NA, raw_slp_hi = NA,
                                   raw_x_low = NA, raw_x_hi = NA,
                                   group = c("on_bout", "off_bout", rep(c("on_bout", "off_bout", "fem_feed", "mal_feed"), 6))) 
              
              ## Number hatched (incubation only since this is before provisioning starts)
                      # On bouts
                              reon$bsn <- as.numeric(reon$brood_size_hatching)
                              reon$nohatch <- reon$clutch_size - reon$bsn
                              reon_m1 <- reon %>% filter(is.na(nohatch) == FALSE, is.na(brood_size_hatching) == FALSE)
                              reon_m1 <- na.omit(reon_m1[, c("nest_key", "bsn", "nohatch", "on_intercept", "on_slope")])
                              reon_m1[reon_m1$nohatch == -1, "nohatch"] <- 0
                  
                              m1 <- glmer(cbind(bsn, nohatch) ~ scale(on_intercept) * scale(on_slope) + (1|nest_key), 
                                          data = reon_m1, family = binomial)
                              m1c <- coefficients(summary(m1))
                              
                                # check residuals. this is repeated for all models below
                                  # xx <- simulateResiduals(m1)
                                  # plot(xx)
                                  # testDispersion(xx)
                              
                              m1_boot <- boot_blup(mi_on, reon_m1[, 1:3], which_re = "binom_x", boots = boot.set)
                              
                              m1r <- glmer(cbind(bsn, nohatch) ~ scale(on_intercept) + scale(on_slope) + (1|nest_key), 
                                           data = reon_m1, family = binomial)
                              m1rc <- coefficients(summary(m1r))
                              
                              m1r_boot <- boot_blup(mi_on, reon_m1[, 1:3], which_re = "binom_sin", boots = boot.set)
                              
                              # add values to table
                                est_df[1, c(5:6, 8:11)] <- c(exp(colMeans(m1r_boot[, 2:3])),
                                                             exp(quantile(m1r_boot[, 4], .05)),
                                                             exp(quantile(m1r_boot[, 5], .95)),
                                                             exp(quantile(m1r_boot[, 6], .05)),
                                                             exp(quantile(m1r_boot[, 7], .95)))
                                est_df[1, c(7, 12:13)] <- c(exp(mean(m1_boot[, 4])),
                                                            exp(quantile(m1_boot[, 9], .05)),
                                                            exp(quantile(m1_boot[, 10], .95)))
                                est_df[1, c(16, 21:22)] <- c(exp(m1c[4, 1]),
                                                             exp(m1c[4, 1] - 1.96 * m1c[4, 2]),
                                                             exp(m1c[4, 1] + 1.96 * m1c[4, 2]))
                                est_df[1, c(14:15, 17:20)] <- c(exp(m1rc[2:3, 1]),
                                                                exp(m1rc[2, 1] - 1.96 * m1rc[2, 2]),
                                                                exp(m1rc[2, 1] + 1.96 * m1rc[2, 2]),
                                                                exp(m1rc[3, 1] - 1.96 * m1rc[3, 2]),
                                                                exp(m1rc[3, 1] + 1.96 * m1rc[3, 2]))
                                
                                
                              
                      # Off bouts
                              reoff$bsn <- as.numeric(reoff$brood_size_hatching)
                              reoff$nohatch <- reoff$clutch_size - reon$bsn
                              reoff_m2 <- reoff %>% filter(is.na(nohatch) == FALSE, is.na(brood_size_hatching) == FALSE)
                              reoff_m2[reoff_m2$nohatch == -1, "nohatch"] <- 0
                              reoff_m2 <- na.omit(reoff_m2[, c("nest_key", "bsn", "nohatch", "off_intercept", "off_slope")])
                              reoff_m2[reoff_m2$nohatch == -1, "nohatch"] <- 0
                              
                              m2 <- glmer(cbind(bsn, nohatch) ~ scale(off_intercept) * scale(off_slope) + (1|nest_key), 
                                          data = reoff_m2, family = binomial)
                              m2c <- coefficients(summary(m2))
                              
                              m2_boot <- boot_blup(mi_off, reoff_m2[, 1:3], which_re = "binom_x", boots = boot.set)
                              
                              m2r <- glmer(cbind(bsn, nohatch) ~ scale(off_intercept) + scale(off_slope) + (1|nest_key), 
                                           data = reoff_m2, family = binomial)
                              m2rc <- coefficients(summary(m2r))
                              
                              m2r_boot <- boot_blup(mi_off, reoff_m2[, 1:3], which_re = "binom_sin", boots = boot.set)
                              
                              # add values to table
                                  est_df[2, c(5:6, 8:11)] <- c(exp(colMeans(m2r_boot[, 2:3])),
                                                               exp(quantile(m2r_boot[, 4], .05)),
                                                               exp(quantile(m2r_boot[, 5], .95)),
                                                               exp(quantile(m2r_boot[, 6], .05)),
                                                               exp(quantile(m2r_boot[, 7], .95)))
                                  est_df[2, c(7, 12:13)] <- c(exp(mean(m2_boot[, 4])),
                                                              exp(quantile(m2_boot[, 9], .05)),
                                                              exp(quantile(m2_boot[, 10], .95)))
                                  est_df[2, c(16, 21:22)] <- c(exp(m2c[4, 1]),
                                                               exp(m2c[4, 1] - 1.96 * m2c[4, 2]),
                                                               exp(m2c[4, 1] + 1.96 * m2c[4, 2]))
                                  est_df[2, c(14:15, 17:20)] <- c(exp(m2rc[2:3, 1]),
                                                                  exp(m2rc[2, 1] - 1.96 * m2rc[2, 2]),
                                                                  exp(m2rc[2, 1] + 1.96 * m2rc[2, 2]),
                                                                  exp(m2rc[3, 1] - 1.96 * m2rc[3, 2]),
                                                                  exp(m2rc[3, 1] + 1.96 * m2rc[3, 2]))
                              
              ## Number alive on day 12 banding
                      # On bouts
                              reonm3 <- reon %>% filter(is.na(bsn) == FALSE, is.na(num_d12) == FALSE)
                              reonm3[reonm3$num_d12 > reonm3$bsn, "num_d12"] <- 4
                              reonm3$died_by12 <- reonm3$bsn - reonm3$num_d12
                              reonm3 <- na.omit(reonm3[, c("nest_key", "num_d12", "died_by12", "on_intercept", "on_slope")])
                              # m3 <- glmer(cbind(num_d12, died_by12) ~ scale(on_intercept) * scale(on_slope) + (1|nest_key),
                              #             data = reonm3, family = binomial)
                            
                            # model is overdispersed. switching to beta-binomial
                              m3bb <- glmmTMB(cbind(num_d12, died_by12) ~ scale(on_intercept) * scale(on_slope) + (1|nest_key),
                                        data = reonm3, family = betabinomial(link = "logit"))
                              m3bbc <- coefficients(summary(m3bb))$cond
                              
                              m3bb_boot <- boot_blup(mi_on, reonm3[, 1:3], which_re = "binom_tmb_x", boots = boot.set2)
                              
                              m3bbr <- glmmTMB(cbind(num_d12, died_by12) ~ scale(on_intercept) + scale(on_slope) + (1|nest_key),
                                              data = reonm3, family = betabinomial(link = "logit"))
                              m3bbrc <- coefficients(summary(m3bbr))$cond
                              
                              m3bbr_boot <- boot_blup(mi_on, reonm3[, 1:3], which_re = "binom_tmb_sin", boots = boot.set)
                              
                              # add values to table
                                est_df[3, c(5:6, 8:11)] <- c(exp(colMeans(m3bbr_boot[, 2:3])),
                                                             exp(quantile(m3bbr_boot[, 4], .05)),
                                                             exp(quantile(m3bbr_boot[, 5], .95)),
                                                             exp(quantile(m3bbr_boot[, 6], .05)),
                                                             exp(quantile(m3bbr_boot[, 7], .95)))
                                est_df[3, c(7, 12:13)] <- c(exp(mean(m3bb_boot[, 4])),
                                                            exp(quantile(m3bb_boot[, 9], .05)),
                                                            exp(quantile(m3bb_boot[, 10], .95)))
                                est_df[3, c(16, 21:22)] <- c(exp(m3bbc[4, 1]),
                                                             exp(m3bbc[4, 1] - 1.96 * m3bbc[4, 2]),
                                                             exp(m3bbc[4, 1] + 1.96 * m3bbc[4, 2]))
                                est_df[3, c(14:15, 17:20)] <- c(exp(m3bbrc[2:3, 1]),
                                                                exp(m3bbrc[2, 1] - 1.96 * m3bbrc[2, 2]),
                                                                exp(m3bbrc[2, 1] + 1.96 * m3bbrc[2, 2]),
                                                                exp(m3bbrc[3, 1] - 1.96 * m3bbrc[3, 2]),
                                                                exp(m3bbrc[3, 1] + 1.96 * m3bbrc[3, 2]))
                              
                      # Off bouts
                              reoffm4 <- reoff %>% filter(is.na(bsn) == FALSE, is.na(num_d12) == FALSE)
                              reoffm4[reoffm4$num_d12 > reoffm4$bsn, "num_d12"] <- 4
                              reoffm4$died_by12 <- reoffm4$bsn - reoffm4$num_d12
                              reoffm4 <- na.omit(reoffm4[, c("nest_key", "num_d12", "died_by12", "off_intercept", "off_slope")])
                              # m4 <- glmer(cbind(num_d12, died_by12) ~ scale(off_intercept) * scale(off_slope) + (1|nest_key),
                              #              data = reoffm4, family = binomial)
                              
                            # model is overdispersed, switching to beta-binomial
                              m4bb <- glmmTMB(cbind(num_d12, died_by12) ~ scale(off_intercept) * scale(off_slope) + (1|nest_key),
                                              data = reoffm4, family = betabinomial(link = "logit"))
                              m4bbc <- coefficients(summary(m4bb))$cond
                              m4bb_boot <- boot_blup(mi_off, reoffm4[, 1:3], which_re = "binom_tmb_x", boots = boot.set2)
                              
                              m4bbr <- glmmTMB(cbind(num_d12, died_by12) ~ scale(off_intercept) + scale(off_slope) + (1|nest_key),
                                              data = reoffm4, family = betabinomial(link = "logit"))
                              m4bbrc <- coefficients(summary(m4bbr))$cond
                              m4bbr_boot <- boot_blup(mi_off, reoffm4[, 1:3], which_re = "binom_tmb_sin", boots = boot.set)
                              
                              # add values to table
                                est_df[4, c(5:6, 8:11)] <- c(exp(colMeans(m4bbr_boot[, 2:3])),
                                                             exp(quantile(m4bbr_boot[, 4], .05)),
                                                             exp(quantile(m4bbr_boot[, 5], .95)),
                                                             exp(quantile(m4bbr_boot[, 6], .05)),
                                                             exp(quantile(m4bbr_boot[, 7], .95)))
                                est_df[4, c(7, 12:13)] <- c(exp(mean(m4bb_boot[, 4])),
                                                            exp(quantile(m4bb_boot[, 9], .05)),
                                                            exp(quantile(m4bb_boot[, 10], .95)))
                                est_df[4, c(16, 21:22)] <- c(exp(m4bbc[4, 1]),
                                                             exp(m4bbc[4, 1] - 1.96 * m4bbc[4, 2]),
                                                             exp(m4bbc[4, 1] + 1.96 * m4bbc[4, 2]))
                                est_df[4, c(14:15, 17:20)] <- c(exp(m4bbrc[2:3, 1]),
                                                                exp(m4bbrc[2, 1] - 1.96 * m4bbrc[2, 2]),
                                                                exp(m4bbrc[2, 1] + 1.96 * m4bbrc[2, 2]),
                                                                exp(m4bbrc[3, 1] - 1.96 * m4bbrc[3, 2]),
                                                                exp(m4bbrc[3, 1] + 1.96 * m4bbrc[3, 2]))
                              
                      # female feeding
                              ref$bsn <- as.numeric(ref$brood_size_hatching)
                              refm5 <- ref
                              refm5 <- refm5 %>% filter(is.na(bsn) == FALSE, is.na(num_d12) == FALSE)
                              refm5[refm5$num_d12 > refm5$bsn, "bsn"] <- 6
                              refm5$died_by12 <- refm5$bsn - refm5$num_d12
                              refm5 <- na.omit(refm5[, c("nest_key", "num_d12", "died_by12", "fem_feed_intercept", "fem_feed_slope")])
                              # m5 <- glmer(cbind(num_d12, died_by12) ~ scale(fem_feed_intercept) * scale(fem_feed_slope) + (1|nest_key),
                              #             data = refm5, family = binomial)
                              #model is overdispersed. switching to beta binomial
                              m5bb <- glmmTMB(cbind(num_d12, died_by12) ~ scale(fem_feed_intercept) * scale(fem_feed_slope) + (1|nest_key),
                                              data = refm5, family = betabinomial(link = "logit"))
                              m5bbc <- coefficients(summary(m5bb))$cond
                              m5bb_boot <- boot_blup(mffeed, refm5[, 1:3], which_re = "binom_tmb_x", boots = boot.set2)
                              
                              m5bbr <- glmmTMB(cbind(num_d12, died_by12) ~ scale(fem_feed_intercept) + scale(fem_feed_slope) + (1|nest_key),
                                              data = refm5, family = betabinomial(link = "logit"))
                              m5bbrc <- coefficients(summary(m5bbr))$cond
                              m5bbr_boot <- boot_blup(mffeed, refm5[, 1:3], which_re = "binom_tmb_sin", boots = boot.set)
                              
                              # add values to table
                                  est_df[5, c(5:6, 8:11)] <- c(exp(colMeans(m5bbr_boot[, 2:3])),
                                                               exp(quantile(m5bbr_boot[, 4], .05)),
                                                               exp(quantile(m5bbr_boot[, 5], .95)),
                                                               exp(quantile(m5bbr_boot[, 6], .05)),
                                                               exp(quantile(m5bbr_boot[, 7], .95)))
                                  est_df[5, c(7, 12:13)] <- c(exp(mean(m5bb_boot[, 4])),
                                                              exp(quantile(m5bb_boot[, 9], .05)),
                                                              exp(quantile(m5bb_boot[, 10], .95)))
                                  est_df[5, c(16, 21:22)] <- c(exp(m5bbc[4, 1]),
                                                               exp(m5bbc[4, 1] - 1.96 * m5bbc[4, 2]),
                                                               exp(m5bbc[4, 1] + 1.96 * m5bbc[4, 2]))
                                  est_df[5, c(14:15, 17:20)] <- c(exp(m5bbrc[2:3, 1]),
                                                                  exp(m5bbrc[2, 1] - 1.96 * m5bbrc[2, 2]),
                                                                  exp(m5bbrc[2, 1] + 1.96 * m5bbrc[2, 2]),
                                                                  exp(m5bbrc[3, 1] - 1.96 * m5bbrc[3, 2]),
                                                                  exp(m5bbrc[3, 1] + 1.96 * m5bbrc[3, 2]))
                              
                      # male feeding
                              rem$bsn <- as.numeric(rem$brood_size_hatching)
                              rem_m6 <- rem
                              rem_m6 <- rem_m6 %>% filter(is.na(bsn) == FALSE, is.na(num_d12) == FALSE)
                              rem_m6$died_by12 <- rem_m6$bsn - rem_m6$num_d12
                              rem_m6 <- na.omit(rem_m6[, c("nest_key", "num_d12", "died_by12", "mal_feed_intercept", "mal_feed_slope")])
                              # m6 <- glmer(cbind(num_d12, died_by12) ~ scale(mal_feed_intercept) * scale(mal_feed_slope) + (1|nest_key),
                              #             data = rem_m6, family = binomial)
                              # model is overdispersed. switching to beta binomial
                              m6bb <- glmmTMB(cbind(num_d12, died_by12) ~ scale(mal_feed_intercept) * scale(mal_feed_slope) + (1|nest_key),
                                              data = rem_m6, family = betabinomial(link = "logit"))
                              m6bbc <- coefficients(summary(m6bb))$cond
                              m6bb_boot <- boot_blup(mmfeed, rem_m6[, 1:3], which_re = "binom_tmb_x", boots = boot.set2)
                              
                              m6bbr <- glmmTMB(cbind(num_d12, died_by12) ~ scale(mal_feed_intercept) + scale(mal_feed_slope) + (1|nest_key),
                                               data = rem_m6, family = betabinomial(link = "logit"))
                              m6bbrc <- coefficients(summary(m6bbr))$cond
                              m6bbr_boot <- boot_blup(mmfeed, rem_m6[, 1:3], which_re = "binom_tmb_sin", boots = boot.set)
                              
                              # add values to table
                                  est_df[6, c(5:6, 8:11)] <- c(exp(colMeans(m6bbr_boot[, 2:3])),
                                                               exp(quantile(m6bbr_boot[, 4], .05)),
                                                               exp(quantile(m6bbr_boot[, 5], .95)),
                                                               exp(quantile(m6bbr_boot[, 6], .05)),
                                                               exp(quantile(m6bbr_boot[, 7], .95)))
                                  est_df[6, c(7, 12:13)] <- c(exp(mean(m6bb_boot[, 4])),
                                                              exp(quantile(m6bb_boot[, 9], .05)),
                                                              exp(quantile(m6bb_boot[, 10], .95)))
                                  est_df[6, c(16, 21:22)] <- c(exp(m6bbc[4, 1]),
                                                               exp(m6bbc[4, 1] - 1.96 * m6bbc[4, 2]),
                                                               exp(m6bbc[4, 1] + 1.96 * m6bbc[4, 2]))
                                  est_df[6, c(14:15, 17:20)] <- c(exp(m6bbrc[2:3, 1]),
                                                                  exp(m6bbrc[2, 1] - 1.96 * m6bbrc[2, 2]),
                                                                  exp(m6bbrc[2, 1] + 1.96 * m6bbrc[2, 2]),
                                                                  exp(m6bbrc[3, 1] - 1.96 * m6bbrc[3, 2]),
                                                                  exp(m6bbrc[3, 1] + 1.96 * m6bbrc[3, 2]))
                              
            ## Number fledged
                  # On bouts
                              reonm7 <- reon %>% filter(is.na(bsn) == FALSE, is.na(num_fledged) == FALSE)
                              reonm7[reonm7$num_fledged > reonm7$bsn, "num_fledged"] <- 4
                              reonm7[reonm7$num_fledged > reonm7$num_d12, "num_fledged"] <- 4
                              reonm7$num_died <- reonm7$bsn - reonm7$num_fledged
                              reonm7[reonm7$num_died < 0, "num_died"] <- 0
                              reonm7 <- na.omit(reonm7[, c("nest_key", "num_fledged", "num_died", "on_intercept", "on_slope")])
                              m7 <- glmer(cbind(num_fledged, num_died) ~ scale(on_intercept) * scale(on_slope) + (1|nest_key),
                                           data = reonm7, family = binomial)
                              m7c <- coefficients(summary(m7))
                              m7_boot <- boot_blup(mi_on, reonm7[, 1:3], which_re = "binom_x", boots = boot.set2)
                              
                              m7r <- glmer(cbind(num_fledged, num_died) ~ scale(on_intercept) + scale(on_slope) + (1|nest_key),
                                          data = reonm7, family = binomial)
                              m7rc <- coefficients(summary(m7r))
                              m7r_boot <- boot_blup(mi_on, reonm7[, 1:3], which_re = "binom_sin", boots = boot.set)
                              
                              # add values to table
                                  est_df[7, c(5:6, 8:11)] <- c(exp(colMeans(m7r_boot[, 2:3])),
                                                               exp(quantile(m7r_boot[, 4], .05)),
                                                               exp(quantile(m7r_boot[, 5], .95)),
                                                               exp(quantile(m7r_boot[, 6], .05)),
                                                               exp(quantile(m7r_boot[, 7], .95)))
                                  est_df[7, c(7, 12:13)] <- c(exp(mean(m7_boot[, 4])),
                                                              exp(quantile(m7_boot[, 9], .05)),
                                                              exp(quantile(m7_boot[, 10], .95)))
                                  est_df[7, c(16, 21:22)] <- c(exp(m7c[4, 1]),
                                                               exp(m7c[4, 1] - 1.96 * m7c[4, 2]),
                                                               exp(m7c[4, 1] + 1.96 * m7c[4, 2]))
                                  est_df[7, c(14:15, 17:20)] <- c(exp(m7rc[2:3, 1]),
                                                                  exp(m7rc[2, 1] - 1.96 * m7rc[2, 2]),
                                                                  exp(m7rc[2, 1] + 1.96 * m7rc[2, 2]),
                                                                  exp(m7rc[3, 1] - 1.96 * m7rc[3, 2]),
                                                                  exp(m7rc[3, 1] + 1.96 * m7rc[3, 2]))
                            
                              
                      # Off bouts
                              reoffm8 <- reoff %>% filter(is.na(bsn) == FALSE, is.na(num_fledged) == FALSE)
                              reoffm8[reoffm8$num_fledged > reoffm8$bsn, "num_fledged"] <- 4
                              reoffm8$num_died <- reoffm8$bsn - reoffm8$num_fledged
                              reoffm8 <- na.omit(reoffm8[, c("nest_key", "num_fledged", "num_died", "off_intercept", "off_slope")])
                              m8 <- glmer(cbind(num_fledged, num_died) ~ scale(off_intercept) * scale(off_slope) + (1|nest_key),
                                          data = reoffm8, family = binomial)
                              m8c <- coefficients(summary(m8))
                              m8_boot <- boot_blup(mi_off, reoffm8[, 1:3], which_re = "binom_x", boots = boot.set2)
                              
                              m8r <- glmer(cbind(num_fledged, num_died) ~ scale(off_intercept) + scale(off_slope) + (1|nest_key),
                                          data = reoffm8, family = binomial)
                              m8rc <- coefficients(summary(m8r))
                              m8r_boot <- boot_blup(mi_off, reoffm8[, 1:3], which_re = "binom_sin", boots = boot.set)
                              
                              # add values to table
                                  est_df[8, c(5:6, 8:11)] <- c(exp(colMeans(na.omit(m8r_boot[, 2:3]))),
                                                               exp(quantile(m8r_boot[, 4], .05, na.rm = TRUE)),
                                                               exp(quantile(m8r_boot[, 5], .95, na.rm = TRUE)),
                                                               exp(quantile(m8r_boot[, 6], .05, na.rm = TRUE)),
                                                               exp(quantile(m8r_boot[, 7], .95, na.rm = TRUE)))
                                  est_df[8, c(7, 12:13)] <- c(exp(mean(na.omit(m8_boot[, 4]))),
                                                              exp(quantile(m8_boot[, 9], .05, na.rm = TRUE)),
                                                              exp(quantile(m8_boot[, 10], .95, na.rm = TRUE)))
                                  est_df[8, c(16, 21:22)] <- c(exp(na.omit(m8c[4, 1])),
                                                               exp(na.omit(m8c[4, 1] - 1.96 * m8c[4, 2])),
                                                               exp(na.omit(m8c[4, 1] + 1.96 * m8c[4, 2])))
                                  est_df[8, c(14:15, 17:20)] <- c(exp(na.omit(m8rc[2:3, 1])),
                                                                  exp(na.omit(m8rc[2, 1] - 1.96 * m8rc[2, 2])),
                                                                  exp(na.omit(m8rc[2, 1] + 1.96 * m8rc[2, 2])),
                                                                  exp(na.omit(m8rc[3, 1] - 1.96 * m8rc[3, 2])),
                                                                  exp(na.omit(m8rc[3, 1] + 1.96 * m8rc[3, 2])))
                              
                      # female feeding
                              ref$bsn <- as.numeric(ref$brood_size_hatching)
                              refm9 <- ref
                              refm9 <- refm9 %>% filter(is.na(bsn) == FALSE, is.na(num_fledged) == FALSE)
                              refm9[refm9$num_fledged > refm9$bsn & refm9$nest_key == "Unit_2_21_2023_133", "num_fledged"] <- 4
                              refm9[refm9$nest_key == "Unit_2_97_2022_129", "bsn"] <- 6
                              refm9$num_died <- refm9$bsn - refm9$num_fledged
                              refm9 <- na.omit(refm9[, c("nest_key", "num_fledged", "num_died", "fem_feed_intercept", "fem_feed_slope")])
                              
                              # m9 <- glmer(cbind(num_fledged, num_died) ~ scale(fem_feed_intercept) * scale(fem_feed_slope) + (1|nest_key),
                              #             data = refm9, family = binomial,
                              #             control = glmerControl(optimizer = "nloptwrap"))
                              #model is overdispersed. switching to beta binomial
                              m9bb <- glmmTMB(cbind(num_fledged, num_died) ~ scale(fem_feed_intercept) * scale(fem_feed_slope) + (1|nest_key),
                                              data = refm9, family = betabinomial(link = "logit"))
                              m9bbc <- coefficients(summary(m9bb))$cond
                              m9bb_boot <- boot_blup(mffeed, refm9[, 1:3], which_re = "binom_tmb_x", boots = boot.set2)
                              
                              m9bbr <- glmmTMB(cbind(num_fledged, num_died) ~ scale(fem_feed_intercept) + scale(fem_feed_slope) + (1|nest_key),
                                              data = refm9, family = betabinomial(link = "logit"))
                              m9bbrc <- coefficients(summary(m9bbr))$cond
                              m9bbr_boot <- boot_blup(mffeed, refm9[, 1:3], which_re = "binom_tmb_sin", boots = boot.set)
                              
                              # add values to table
                                  est_df[9, c(5:6, 8:11)] <- c(exp(colMeans(m9bbr_boot[, 2:3])),
                                                               exp(quantile(m9bbr_boot[, 4], .05, na.rm = TRUE)),
                                                               exp(quantile(m9bbr_boot[, 5], .95, na.rm = TRUE)),
                                                               exp(quantile(m9bbr_boot[, 6], .05, na.rm = TRUE)),
                                                               exp(quantile(m9bbr_boot[, 7], .95, na.rm = TRUE)))
                                  est_df[9, c(7, 12:13)] <- c(exp(mean(m9bb_boot[, 4])),
                                                              exp(quantile(m9bb_boot[, 9], .05, na.rm = TRUE)),
                                                              exp(quantile(m9bb_boot[, 10], .95, na.rm = TRUE)))
                                  est_df[9, c(16, 21:22)] <- c(exp(m9bbc[4, 1]),
                                                               exp(m9bbc[4, 1] - 1.96 * m9bbc[4, 2]),
                                                               exp(m9bbc[4, 1] + 1.96 * m9bbc[4, 2]))
                                  est_df[9, c(14:15, 17:20)] <- c(exp(m9bbrc[2:3, 1]),
                                                                  exp(m9bbrc[2, 1] - 1.96 * m9bbrc[2, 2]),
                                                                  exp(m9bbrc[2, 1] + 1.96 * m9bbrc[2, 2]),
                                                                  exp(m9bbrc[3, 1] - 1.96 * m9bbrc[3, 2]),
                                                                  exp(m9bbrc[3, 1] + 1.96 * m9bbrc[3, 2]))
                              
                      # male feeding
                              rem$bsn <- as.numeric(rem$brood_size_hatching)
                              rem_m10 <- rem
                              rem_m10$success <- rem_m10$num_fledged
                              rem_m10$failure <- rem_m10$bsn - rem_m10$num_fledged
                              rem_m10 <- na.omit(rem_m10[, c("nest_key", "success", "failure", "mal_feed_intercept", "mal_feed_slope")])
                              # m10 <- glmer(cbind(success, failure) ~ scale(mal_feed_intercept) * scale(mal_feed_slope) + (1|nest_key),
                              #             data = rem_m10, family = binomial)
                              # model is overdispersed. switching to beta binomial
                              m10bb <- glmmTMB(cbind(success, failure) ~ scale(mal_feed_intercept) * scale(mal_feed_slope) + (1|nest_key),
                                              data = rem_m10, family = betabinomial(link = "logit"))
                              m10bbc <- coefficients(summary(m10bb))$cond
                              m10bb_boot <- boot_blup(mmfeed, rem_m10[, 1:3], which_re = "binom_tmb_x", boots = boot.set2)
                              
                              m10bbr <- glmmTMB(cbind(success, failure) ~ scale(mal_feed_intercept) + scale(mal_feed_slope) + (1|nest_key),
                                               data = rem_m10, family = betabinomial(link = "logit"))
                              m10bbrc <- coefficients(summary(m10bbr))$cond
                              m10bbr_boot <- boot_blup(mmfeed, rem_m10[, 1:3], which_re = "binom_tmb_sin", boots = boot.set)
                              rem_m10_s <- rem_m10[, c("nest_key", "success", "failure")]
                              
                              # add values to table
                                  est_df[10, c(5:6, 8:11)] <- c(exp(colMeans(m10bbr_boot[, 2:3])),
                                                               exp(quantile(m10bbr_boot[, 4], .05, na.rm = TRUE)),
                                                               exp(quantile(m10bbr_boot[, 5], .95, na.rm = TRUE)),
                                                               exp(quantile(m10bbr_boot[, 6], .05, na.rm = TRUE)),
                                                               exp(quantile(m10bbr_boot[, 7], .95, na.rm = TRUE)))
                                  est_df[10, c(7, 12:13)] <- c(exp(mean(m10bb_boot[, 4])),
                                                              exp(quantile(m10bb_boot[, 9], .05, na.rm = TRUE)),
                                                              exp(quantile(m10bb_boot[, 10], .95, na.rm = TRUE)))
                                  est_df[10, c(16, 21:22)] <- c(exp(m10bbc[4, 1]),
                                                               exp(m10bbc[4, 1] - 1.96 * m10bbc[4, 2]),
                                                               exp(m10bbc[4, 1] + 1.96 * m10bbc[4, 2]))
                                  est_df[10, c(14:15, 17:20)] <- c(exp(m10bbrc[2:3, 1]),
                                                                  exp(m10bbrc[2, 1] - 1.96 * m10bbrc[2, 2]),
                                                                  exp(m10bbrc[2, 1] + 1.96 * m10bbrc[2, 2]),
                                                                  exp(m10bbrc[3, 1] - 1.96 * m10bbrc[3, 2]),
                                                                  exp(m10bbrc[3, 1] + 1.96 * m10bbrc[3, 2]))
                                  
            ## recruits next year (added later so numbering skips)
                      # on bout length
                              reon$bsn <- as.numeric(reon$brood_size_hatching)
                              reon_m23 <- reon
                              reon_m23 <- reon_m23 %>% filter(is.na(bsn) == FALSE, is.na(num_fledged) == FALSE)
                              reon_m23$num_no_recruit <- reon_m23$num_fledged - reon_m23$recruits
                              reon_m23 <- na.omit(reon_m23[, c("nest_key", "recruits", "num_no_recruit", "on_intercept", "on_slope")])
                              
                              m23 <- glmer(cbind(recruits, num_no_recruit) ~ scale(on_intercept) * scale(on_slope) + (1|nest_key),
                                           data = reon_m23, family = binomial)
                              m23c <- coefficients(summary(m23))
                              m23_boot <- boot_blup(mi_on, reon_m23[, 1:3], which_re = "binom_x", boots = boot.set2)
                              
                              m23r <- glmer(cbind(recruits, num_no_recruit) ~ scale(on_intercept) + scale(on_slope) + (1|nest_key),
                                           data = reon_m23, family = binomial)
                              m23rc <- coefficients(summary(m23r))
                              m23r_boot <- boot_blup(mi_on, reon_m23[, 1:3], which_re = "binom_sin", boots = boot.set)
                              
                              # add values to table
                                  est_df[23, c(5:6, 8:11)] <- c(exp(colMeans(m23r_boot[, 2:3])),
                                                                exp(quantile(m23r_boot[, 4], .05)),
                                                                exp(quantile(m23r_boot[, 5], .95)),
                                                                exp(quantile(m23r_boot[, 6], .05)),
                                                                exp(quantile(m23r_boot[, 7], .95)))
                                  est_df[23, c(7, 12:13)] <- c(exp(mean(m23_boot[, 4])),
                                                               exp(quantile(m23_boot[, 9], .05)),
                                                               exp(quantile(m23_boot[, 10], .95)))
                                  est_df[23, c(16, 21:22)] <- c(exp(m23c[4, 1]),
                                                                exp(m23c[4, 1] - 1.96 * m23c[4, 2]),
                                                                exp(m23c[4, 1] + 1.96 * m23c[4, 2]))
                                  est_df[23, c(14:15, 17:20)] <- c(exp(m23rc[2:3, 1]),
                                                                   exp(m23rc[2, 1] - 1.96 * m23rc[2, 2]),
                                                                   exp(m23rc[2, 1] + 1.96 * m23rc[2, 2]),
                                                                   exp(m23rc[3, 1] - 1.96 * m23rc[3, 2]),
                                                                   exp(m23rc[3, 1] + 1.96 * m23rc[3, 2]))
                                  
                      # off bout length
                              reoff$bsn <- as.numeric(reoff$brood_size_hatching)
                              reoff_m24 <- reoff
                              reoff_m24 <- reoff_m24 %>% filter(is.na(bsn) == FALSE, is.na(num_fledged) == FALSE)
                              reoff_m24$num_no_recruit <- reoff_m24$num_fledged - reoff_m24$recruits
                              reoff_m24 <- na.omit(reoff_m24[, c("nest_key", "recruits", "num_no_recruit", "off_intercept", "off_slope")])
                              
                              m24 <- glmer(cbind(recruits, num_no_recruit) ~ scale(off_intercept) * scale(off_slope) + (1|nest_key),
                                           data = reoff_m24, family = binomial)
                              m24c <- coefficients(summary(m24))
                              m24_boot <- boot_blup(mi_off, reoff_m24[, 1:3], which_re = "binom_x", boots = boot.set2)
                              
                              m24r <- glmer(cbind(recruits, num_no_recruit) ~ scale(off_intercept) + scale(off_slope) + (1|nest_key),
                                           data = reoff_m24, family = binomial)
                              m24rc <- coefficients(summary(m24r))
                              m24r_boot <- boot_blup(mi_off, reoff_m24[, 1:3], which_re = "binom_sin", boots = boot.set)
                              
                              # add values to table
                                  est_df[24, c(5:6, 8:11)] <- c(exp(colMeans(m24r_boot[, 2:3])),
                                                                exp(quantile(m24r_boot[, 4], .05)),
                                                                exp(quantile(m24r_boot[, 5], .95)),
                                                                exp(quantile(m24r_boot[, 6], .05)),
                                                                exp(quantile(m24r_boot[, 7], .95)))
                                  est_df[24, c(7, 12:13)] <- c(exp(mean(m24_boot[, 4])),
                                                               exp(quantile(m24_boot[, 9], .05)),
                                                               exp(quantile(m24_boot[, 10], .95)))
                                  est_df[24, c(16, 21:22)] <- c(exp(m24c[4, 1]),
                                                                exp(m24c[4, 1] - 1.96 * m24c[4, 2]),
                                                                exp(m24c[4, 1] + 1.96 * m24c[4, 2]))
                                  est_df[24, c(14:15, 17:20)] <- c(exp(m24rc[2:3, 1]),
                                                                   exp(m24rc[2, 1] - 1.96 * m24rc[2, 2]),
                                                                   exp(m24rc[2, 1] + 1.96 * m24rc[2, 2]),
                                                                   exp(m24rc[3, 1] - 1.96 * m24rc[3, 2]),
                                                                   exp(m24rc[3, 1] + 1.96 * m24rc[3, 2]))    
                              
                                  
                      # female feeding rate  
                              ref$bsn <- as.numeric(ref$brood_size_hatching)
                              refm25 <- ref
                              refm25 <- refm25 %>% filter(is.na(bsn) == FALSE, is.na(num_fledged) == FALSE)
                              refm25[refm25$num_fledged > refm25$bsn & refm25$nest_key == "Unit_2_21_2023_133", "num_fledged"] <- 4
                              refm25[refm25$nest_key == "Unit_2_97_2022_129", "bsn"] <- 6
                              refm25$num_no_recruit <- refm25$num_fledged - refm25$recruits
                              refm25 <- na.omit(refm25[, c("nest_key", "recruits", "num_no_recruit", "fem_feed_intercept", "fem_feed_slope")])
                              
                              m25 <- glmer(cbind(recruits, num_no_recruit) ~ scale(fem_feed_intercept) * scale(fem_feed_slope) + (1|nest_key),
                                          data = refm25, family = binomial)
                              m25c <- coefficients(summary(m25))
                              m25_boot <- boot_blup(mffeed, refm25[, 1:3], which_re = "binom_x", boots = boot.set2)
                              
                              m25r <- glmer(cbind(recruits, num_no_recruit) ~ scale(fem_feed_intercept) + scale(fem_feed_slope) + (1|nest_key),
                                            data = refm25, family = binomial)
                              m25rc <- coefficients(summary(m25r))
                              m25r_boot <- boot_blup(mffeed, refm25[, 1:3], which_re = "binom_sin", boots = boot.set)
                              
                              # add values to table
                                  est_df[25, c(5:6, 8:11)] <- c(exp(colMeans(m25r_boot[, 2:3])),
                                                                exp(quantile(m25r_boot[, 4], .05)),
                                                                exp(quantile(m25r_boot[, 5], .95)),
                                                                exp(quantile(m25r_boot[, 6], .05)),
                                                                exp(quantile(m25r_boot[, 7], .95)))
                                  est_df[25, c(7, 12:13)] <- c(exp(mean(m25_boot[, 4])),
                                                               exp(quantile(m25_boot[, 9], .05)),
                                                               exp(quantile(m25_boot[, 10], .95)))
                                  est_df[25, c(16, 21:22)] <- c(exp(m25c[4, 1]),
                                                                exp(m25c[4, 1] - 1.96 * m25c[4, 2]),
                                                                exp(m25c[4, 1] + 1.96 * m25c[4, 2]))
                                  est_df[25, c(14:15, 17:20)] <- c(exp(m25rc[2:3, 1]),
                                                                   exp(m25rc[2, 1] - 1.96 * m25rc[2, 2]),
                                                                   exp(m25rc[2, 1] + 1.96 * m25rc[2, 2]),
                                                                   exp(m25rc[3, 1] - 1.96 * m25rc[3, 2]),
                                                                   exp(m25rc[3, 1] + 1.96 * m25rc[3, 2]))
             
                      # male feeding rate
                                rem$bsn <- as.numeric(rem$brood_size_hatching)
                                  remm26 <- rem
                                  remm26 <- remm26 %>% filter(is.na(bsn) == FALSE, is.na(num_fledged) == FALSE)
                                  remm26$num_no_recruit <- remm26$num_fledged - remm26$recruits
                                  remm26 <- na.omit(remm26[, c("nest_key", "recruits", "num_no_recruit", "mal_feed_intercept", "mal_feed_slope")])
                                  
                                  m26 <- glmer(cbind(recruits, num_no_recruit) ~ scale(mal_feed_intercept) * scale(mal_feed_slope) + (1|nest_key),
                                              data = remm26, family = binomial)
                                  m26c <- coefficients(summary(m26))
                                  m26_boot <- boot_blup(mmfeed, remm26[, 1:3], which_re = "binom_x", boots = boot.set2)
                                  
                                  m26r <- glmer(cbind(recruits, num_no_recruit) ~ scale(mal_feed_intercept) + scale(mal_feed_slope) + (1|nest_key),
                                                data = remm26, family = binomial)
                                  m26rc <- coefficients(summary(m26r))
                                  m26r_boot <- boot_blup(mmfeed, remm26[, 1:3], which_re = "binom_sin", boots = boot.set)
                              
                              # add values to table
                                  est_df[26, c(5:6, 8:11)] <- c(exp(colMeans(m26r_boot[, 2:3])),
                                                                exp(quantile(m26r_boot[, 4], .05)),
                                                                exp(quantile(m26r_boot[, 5], .95)),
                                                                exp(quantile(m26r_boot[, 6], .05)),
                                                                exp(quantile(m26r_boot[, 7], .95)))
                                  est_df[26, c(7, 12:13)] <- c(exp(mean(m26_boot[, 4])),
                                                               exp(quantile(m26_boot[, 9], .05)),
                                                               exp(quantile(m26_boot[, 10], .95)))
                                  est_df[26, c(16, 21:22)] <- c(exp(m26c[4, 1]),
                                                                exp(m26c[4, 1] - 1.96 * m26c[4, 2]),
                                                                exp(m26c[4, 1] + 1.96 * m26c[4, 2]))
                                  est_df[26, c(14:15, 17:20)] <- c(exp(m26rc[2:3, 1]),
                                                                   exp(m26rc[2, 1] - 1.96 * m26rc[2, 2]),
                                                                   exp(m26rc[2, 1] + 1.96 * m26rc[2, 2]),
                                                                   exp(m26rc[3, 1] - 1.96 * m26rc[3, 2]),
                                                                   exp(m26rc[3, 1] + 1.96 * m26rc[3, 2]))
                                  
                              
            ## Average day 12 mass 
                      # on bouts
                              m11 <- lm(scale(d12_av_mass) ~ scale(on_intercept) * scale(on_slope), data = reon)
                              m11c <- coefficients(summary(m11))
                              m11_boot <- boot_blup(mi_on, na.omit(reon[, c("nest_key", "d12_av_mass")]), which_re = "both_inter", boots = boot.set2) 
                              m11r <- lm(scale(d12_av_mass) ~ scale(on_intercept) + scale(on_slope), data = reon)
                              m11rc <- coefficients(summary(m11r))
                              m11r_boot <- boot_blup(mi_on, na.omit(reon[, c("nest_key", "d12_av_mass")]), which_re = "both_single", boots = boot.set)
                              
                              # add values to table
                                  est_df[11, c(5:6, 8:11)] <- c(colMeans(m11r_boot[, 2:3]),
                                                                quantile(m11r_boot[, 4], .05),
                                                                quantile(m11r_boot[, 5], .95),
                                                                quantile(m11r_boot[, 6], .05),
                                                                quantile(m11r_boot[, 7], .95))
                                  est_df[11, c(7, 12:13)] <- c(mean(m11_boot[, 4]),
                                                               quantile(m11_boot[, 9], .05),
                                                               quantile(m11_boot[, 10], .95))
                                  est_df[11, c(16, 21:22)] <- c(m11c[4, 1],
                                                                m11c[4, 1] - 1.96 * m11c[4, 2],
                                                                m11c[4, 1] + 1.96 * m11c[4, 2])
                                  est_df[11, c(14:15, 17:20)] <- c(m11rc[2:3, 1],
                                                                   m11rc[2, 1] - 1.96 * m11rc[2, 2],
                                                                   m11rc[2, 1] + 1.96 * m11rc[2, 2],
                                                                   m11rc[3, 1] - 1.96 * m11rc[3, 2],
                                                                   m11rc[3, 1] + 1.96 * m11rc[3, 2])
                      # off bouts
                              m12 <- lm(scale(d12_av_mass) ~ scale(off_intercept) * scale(off_slope), data = reoff)
                              m12c <- coefficients(summary(m12))
                              m12_boot <- boot_blup(mi_off, na.omit(reoff[, c("nest_key", "d12_av_mass")]), which_re = "both_inter", boots = boot.set2) 
                              m12r <- lm(scale(d12_av_mass) ~ scale(off_intercept) + scale(off_slope), data = reoff)
                              m12rc <- coefficients(summary(m12r))
                              m12r_boot <- boot_blup(mi_off, na.omit(reoff[, c("nest_key", "d12_av_mass")]), which_re = "both_single", boots = boot.set) 
                              
                              # add values to table
                                  est_df[12, c(5:6, 8:11)] <- c(colMeans(m12r_boot[, 2:3]),
                                                                quantile(m12r_boot[, 4], .05),
                                                                quantile(m12r_boot[, 5], .95),
                                                                quantile(m12r_boot[, 6], .05),
                                                                quantile(m12r_boot[, 7], .95))
                                  est_df[12, c(7, 12:13)] <- c(mean(m12_boot[, 4]),
                                                               quantile(m12_boot[, 9], .05),
                                                               quantile(m12_boot[, 10], .95))
                                  est_df[12, c(16, 21:22)] <- c(m12c[4, 1],
                                                                m12c[4, 1] - 1.96 * m12c[4, 2],
                                                                m12c[4, 1] + 1.96 * m12c[4, 2])
                                  est_df[12, c(14:15, 17:20)] <- c(m12rc[2:3, 1],
                                                                   m12rc[2, 1] - 1.96 * m12rc[2, 2],
                                                                   m12rc[2, 1] + 1.96 * m12rc[2, 2],
                                                                   m12rc[3, 1] - 1.96 * m12rc[3, 2],
                                                                   m12rc[3, 1] + 1.96 * m12rc[3, 2])
                      # female feeding
                              m13 <- lm(scale(d12_av_mass) ~ scale(fem_feed_intercept) * scale(fem_feed_slope), data = ref)
                              m13c <- coefficients(summary(m13))
                              m13_boot <- boot_blup(mffeed, na.omit(ref[, c("nest_key", "d12_av_mass")]), which_re = "both_inter", boots = boot.set2)
                              m13r <- lm(scale(d12_av_mass) ~ scale(fem_feed_intercept) + scale(fem_feed_slope), data = ref)
                              m13rc <- coefficients(summary(m13r))
                              m13r_boot <- boot_blup(mffeed, na.omit(ref[, c("nest_key", "d12_av_mass")]), which_re = "both_single", boots = boot.set)
                              
                              # add values to table
                                  est_df[13, c(5:6, 8:11)] <- c(colMeans(m13r_boot[, 2:3]),
                                                                quantile(m13r_boot[, 4], .05),
                                                                quantile(m13r_boot[, 5], .95),
                                                                quantile(m13r_boot[, 6], .05),
                                                                quantile(m13r_boot[, 7], .95))
                                  est_df[13, c(7, 12:13)] <- c(mean(m13_boot[, 4]),
                                                               quantile(m13_boot[, 9], .05),
                                                               quantile(m13_boot[, 10], .95))
                                  est_df[13, c(16, 21:22)] <- c(m13c[4, 1],
                                                                m13c[4, 1] - 1.96 * m13c[4, 2],
                                                                m13c[4, 1] + 1.96 * m13c[4, 2])
                                  est_df[13, c(14:15, 17:20)] <- c(m13rc[2:3, 1],
                                                                   m13rc[2, 1] - 1.96 * m13rc[2, 2],
                                                                   m13rc[2, 1] + 1.96 * m13rc[2, 2],
                                                                   m13rc[3, 1] - 1.96 * m13rc[3, 2],
                                                                   m13rc[3, 1] + 1.96 * m13rc[3, 2])
                      # male feeding
                              m14 <- lm(scale(d12_av_mass) ~ scale(mal_feed_intercept) * scale(mal_feed_slope), data = rem)
                              m14c <- coefficients(summary(m14))
                              m14_boot <- boot_blup(mmfeed, na.omit(ref[, c("nest_key", "d12_av_mass")]), which_re = "both_inter", boots = boot.set2)
                              m14r <- lm(scale(d12_av_mass) ~ scale(mal_feed_intercept) + scale(mal_feed_slope), data = rem)
                              m14rc <- coefficients(summary(m14r))
                              m14r_boot <- boot_blup(mmfeed, na.omit(ref[, c("nest_key", "d12_av_mass")]), which_re = "both_single", boots = boot.set)
                              
                              # add values to table
                                  est_df[14, c(5:6, 8:11)] <- c(colMeans(m14r_boot[, 2:3]),
                                                                quantile(m14r_boot[, 4], .05),
                                                                quantile(m14r_boot[, 5], .95),
                                                                quantile(m14r_boot[, 6], .05),
                                                                quantile(m14r_boot[, 7], .95))
                                  est_df[14, c(7, 12:13)] <- c(mean(m14_boot[, 4]),
                                                               quantile(m14_boot[, 9], .05),
                                                               quantile(m14_boot[, 10], .95))
                                  est_df[14, c(16, 21:22)] <- c(m14c[4, 1],
                                                                m14c[4, 1] - 1.96 * m14c[4, 2],
                                                                m14c[4, 1] + 1.96 * m14c[4, 2])
                                  est_df[14, c(14:15, 17:20)] <- c(m14rc[2:3, 1],
                                                                   m14rc[2, 1] - 1.96 * m14rc[2, 2],
                                                                   m14rc[2, 1] + 1.96 * m14rc[2, 2],
                                                                   m14rc[3, 1] - 1.96 * m14rc[3, 2],
                                                                   m14rc[3, 1] + 1.96 * m14rc[3, 2])
              
            ## Average day 12 headbill
                      # on bouts
                              m15 <- lm(scale(d12_av_hb) ~ scale(on_intercept) * scale(on_slope), data = reon)
                              m15c <- coefficients(summary(m15))
                              m15_boot <- boot_blup(mi_on, na.omit(reon[, c("nest_key", "d12_av_hb")]), which_re = "both_inter", boots = boot.set2) 
                              m15r <- lm(scale(d12_av_hb) ~ scale(on_intercept) + scale(on_slope), data = reon)
                              m15rc <- coefficients(summary(m15r))
                              m15r_boot <- boot_blup(mi_on, na.omit(reon[, c("nest_key", "d12_av_hb")]), which_re = "both_single", boots = boot.set)
                              
                              # add values to table
                                  est_df[15, c(5:6, 8:11)] <- c(colMeans(m15r_boot[, 2:3]),
                                                                quantile(m15r_boot[, 4], .05),
                                                                quantile(m15r_boot[, 5], .95),
                                                                quantile(m15r_boot[, 6], .05),
                                                                quantile(m15r_boot[, 7], .95))
                                  est_df[15, c(7, 12:13)] <- c(mean(m15_boot[, 4]),
                                                               quantile(m15_boot[, 9], .05),
                                                               quantile(m15_boot[, 10], .95))
                                  est_df[15, c(16, 21:22)] <- c(m15c[4, 1],
                                                                m15c[4, 1] - 1.96 * m15c[4, 2],
                                                                m15c[4, 1] + 1.96 * m15c[4, 2])
                                  est_df[15, c(14:15, 17:20)] <- c(m15rc[2:3, 1],
                                                                   m15rc[2, 1] - 1.96 * m15rc[2, 2],
                                                                   m15rc[2, 1] + 1.96 * m15rc[2, 2],
                                                                   m15rc[3, 1] - 1.96 * m15rc[3, 2],
                                                                   m15rc[3, 1] + 1.96 * m15rc[3, 2])
                      # off bouts
                              m16 <- lm(scale(d12_av_hb) ~ scale(off_intercept) * scale(off_slope), data = reoff)
                              m16c <- coefficients(summary(m16))
                              m16_boot <- boot_blup(mi_off, na.omit(reoff[, c("nest_key", "d12_av_hb")]), which_re = "both_inter", boots = boot.set2)
                              m16r <- lm(scale(d12_av_hb) ~ scale(off_intercept) + scale(off_slope), data = reoff)
                              m16rc <- coefficients(summary(m16r))
                              m16r_boot <- boot_blup(mi_off, na.omit(reoff[, c("nest_key", "d12_av_hb")]), which_re = "both_single", boots = boot.set)
                              
                              # add values to table
                                  est_df[16, c(5:6, 8:11)] <- c(colMeans(m16r_boot[, 2:3]),
                                                                quantile(m16r_boot[, 4], .05),
                                                                quantile(m16r_boot[, 5], .95),
                                                                quantile(m16r_boot[, 6], .05),
                                                                quantile(m16r_boot[, 7], .95))
                                  est_df[16, c(7, 12:13)] <- c(mean(m16_boot[, 4]),
                                                               quantile(m16_boot[, 9], .05),
                                                               quantile(m16_boot[, 10], .95))
                                  est_df[16, c(16, 21:22)] <- c(m16c[4, 1],
                                                                m16c[4, 1] - 1.96 * m16c[4, 2],
                                                                m16c[4, 1] + 1.96 * m16c[4, 2])
                                  est_df[16, c(14:15, 17:20)] <- c(m16rc[2:3, 1],
                                                                   m16rc[2, 1] - 1.96 * m16rc[2, 2],
                                                                   m16rc[2, 1] + 1.96 * m16rc[2, 2],
                                                                   m16rc[3, 1] - 1.96 * m16rc[3, 2],
                                                                   m16rc[3, 1] + 1.96 * m16rc[3, 2])
                      # female feeding
                              m17 <- lm(scale(d12_av_hb) ~ scale(fem_feed_intercept) * scale(fem_feed_slope), data = ref)
                              m17c <- coefficients(summary(m17))
                              m17_boot <- boot_blup(mffeed, na.omit(ref[, c("nest_key", "d12_av_hb")]), which_re = "both_inter", boots = boot.set2)
                              m17r <- lm(scale(d12_av_hb) ~ scale(fem_feed_intercept) + scale(fem_feed_slope), data = ref)
                              m17rc <- coefficients(summary(m17r))
                              m17r_boot <- boot_blup(mffeed, na.omit(ref[, c("nest_key", "d12_av_hb")]), which_re = "both_single", boots = boot.set)
                              
                              # add values to table
                                  est_df[17, c(5:6, 8:11)] <- c(colMeans(m17r_boot[, 2:3]),
                                                                quantile(m17r_boot[, 4], .05),
                                                                quantile(m17r_boot[, 5], .95),
                                                                quantile(m17r_boot[, 6], .05),
                                                                quantile(m17r_boot[, 7], .95))
                                  est_df[17, c(7, 12:13)] <- c(mean(m17_boot[, 4]),
                                                               quantile(m17_boot[, 9], .05),
                                                               quantile(m17_boot[, 10], .95))
                                  est_df[17, c(16, 21:22)] <- c(m17c[4, 1],
                                                                m17c[4, 1] - 1.96 * m17c[4, 2],
                                                                m17c[4, 1] + 1.96 * m17c[4, 2])
                                  est_df[17, c(14:15, 17:20)] <- c(m17rc[2:3, 1],
                                                                   m17rc[2, 1] - 1.96 * m17rc[2, 2],
                                                                   m17rc[2, 1] + 1.96 * m17rc[2, 2],
                                                                   m17rc[3, 1] - 1.96 * m17rc[3, 2],
                                                                   m17rc[3, 1] + 1.96 * m17rc[3, 2])
                      # male feeding
                              m18 <- lm(scale(d12_av_hb) ~ scale(mal_feed_intercept) * scale(mal_feed_slope), data = rem)
                              m18c <- coefficients(summary(m18))
                              m18_boot <- boot_blup(mmfeed, na.omit(ref[, c("nest_key", "d12_av_hb")]), which_re = "both_inter", boots = boot.set2)
                              m18r <- lm(scale(d12_av_hb) ~ scale(mal_feed_intercept) + scale(mal_feed_slope), data = rem)     
                              m18rc <- coefficients(summary(m18r))
                              m18r_boot <- boot_blup(mmfeed, na.omit(ref[, c("nest_key", "d12_av_hb")]), which_re = "both_single", boots = boot.set)
                              
                              # add values to table
                                  est_df[18, c(5:6, 8:11)] <- c(colMeans(m18r_boot[, 2:3]),
                                                                quantile(m18r_boot[, 4], .05),
                                                                quantile(m18r_boot[, 5], .95),
                                                                quantile(m18r_boot[, 6], .05),
                                                                quantile(m18r_boot[, 7], .95))
                                  est_df[18, c(7, 12:13)] <- c(mean(m18_boot[, 4]),
                                                               quantile(m18_boot[, 9], .05),
                                                               quantile(m18_boot[, 10], .95))
                                  est_df[18, c(16, 21:22)] <- c(m18c[4, 1],
                                                                m18c[4, 1] - 1.96 * m18c[4, 2],
                                                                m18c[4, 1] + 1.96 * m18c[4, 2])
                                  est_df[18, c(14:15, 17:20)] <- c(m18rc[2:3, 1],
                                                                   m18rc[2, 1] - 1.96 * m18rc[2, 2],
                                                                   m18rc[2, 1] + 1.96 * m18rc[2, 2],
                                                                   m18rc[3, 1] - 1.96 * m18rc[3, 2],
                                                                   m18rc[3, 1] + 1.96 * m18rc[3, 2])

            ## Average day 12 flatwing
                      # on bouts
                              m19 <- lm(scale(d12_av_fw) ~ scale(on_intercept) * scale(on_slope), data = reon)
                              m19c <- coefficients(summary(m19))
                              m19_boot <- boot_blup(mi_on, na.omit(reon[, c("nest_key", "d12_av_fw")]), which_re = "both_inter", boots = boot.set2)
                              m19r <- lm(scale(d12_av_fw) ~ scale(on_intercept) + scale(on_slope), data = reon)
                              m19rc <- coefficients(summary(m19r))
                              m19r_boot <- boot_blup(mi_on, na.omit(reon[, c("nest_key", "d12_av_fw")]), which_re = "both_single", boots = boot.set)
                              
                              # add values to table
                                  est_df[19, c(5:6, 8:11)] <- c(colMeans(m19r_boot[, 2:3]),
                                                                quantile(m19r_boot[, 4], .05),
                                                                quantile(m19r_boot[, 5], .95),
                                                                quantile(m19r_boot[, 6], .05),
                                                                quantile(m19r_boot[, 7], .95))
                                  est_df[19, c(7, 12:13)] <- c(mean(m19_boot[, 4]),
                                                               quantile(m19_boot[, 9], .05),
                                                               quantile(m19_boot[, 10], .95))
                                  est_df[19, c(16, 21:22)] <- c(m19c[4, 1],
                                                                m19c[4, 1] - 1.96 * m19c[4, 2],
                                                                m19c[4, 1] + 1.96 * m19c[4, 2])
                                  est_df[19, c(14:15, 17:20)] <- c(m19rc[2:3, 1],
                                                                   m19rc[2, 1] - 1.96 * m19rc[2, 2],
                                                                   m19rc[2, 1] + 1.96 * m19rc[2, 2],
                                                                   m19rc[3, 1] - 1.96 * m19rc[3, 2],
                                                                   m19rc[3, 1] + 1.96 * m19rc[3, 2])
                      # off bouts
                              m20 <- lm(scale(d12_av_fw) ~ scale(off_intercept) * scale(off_slope), data = reoff)
                              m20c <- coefficients(summary(m20))
                              m20_boot <- boot_blup(mi_off, na.omit(reoff[, c("nest_key", "d12_av_fw")]), which_re = "both_inter", boots = boot.set2)
                              m20r <- lm(scale(d12_av_fw) ~ scale(off_intercept) + scale(off_slope), data = reoff)
                              m20rc <- coefficients(summary(m20r))
                              m20r_boot <- boot_blup(mi_off, na.omit(reoff[, c("nest_key", "d12_av_fw")]), which_re = "both_single", boots = boot.set)
                              
                              # add values to table
                                  est_df[20, c(5:6, 8:11)] <- c(colMeans(m20r_boot[, 2:3]),
                                                                quantile(m20r_boot[, 4], .05),
                                                                quantile(m20r_boot[, 5], .95),
                                                                quantile(m20r_boot[, 6], .05),
                                                                quantile(m20r_boot[, 7], .95))
                                  est_df[20, c(7, 12:13)] <- c(mean(m20_boot[, 4]),
                                                               quantile(m20_boot[, 9], .05),
                                                               quantile(m20_boot[, 10], .95))
                                  est_df[20, c(16, 21:22)] <- c(m20c[4, 1],
                                                                m20c[4, 1] - 1.96 * m20c[4, 2],
                                                                m20c[4, 1] + 1.96 * m20c[4, 2])
                                  est_df[20, c(14:15, 17:20)] <- c(m20rc[2:3, 1],
                                                                   m20rc[2, 1] - 1.96 * m20rc[2, 2],
                                                                   m20rc[2, 1] + 1.96 * m20rc[2, 2],
                                                                   m20rc[3, 1] - 1.96 * m20rc[3, 2],
                                                                   m20rc[3, 1] + 1.96 * m20rc[3, 2])
                      # female feeding
                              m21 <- lm(scale(d12_av_fw) ~ scale(fem_feed_intercept) * scale(fem_feed_slope), data = ref)
                              m21c <- coefficients(summary(m21))
                              m21_boot <- boot_blup(mffeed, na.omit(ref[, c("nest_key", "d12_av_fw")]), which_re = "both_inter", boots = boot.set2)
                              ref21 <- na.omit(ref[, c("nest_key", "d12_av_fw", "fem_feed_intercept", "fem_feed_slope")])
                              m21r <- lm(scale(d12_av_fw) ~ scale(fem_feed_intercept) + scale(fem_feed_slope), data = ref21)
                              m21rc <- coefficients(summary(m21r))
                              m21r_boot <- boot_blup(mffeed, na.omit(ref[, c("nest_key", "d12_av_fw")]), which_re = "both_single", boots = boot.set)
                              
                              # add values to table
                                  est_df[21, c(5:6, 8:11)] <- c(colMeans(m21r_boot[, 2:3]),
                                                                quantile(m21r_boot[, 4], .05),
                                                                quantile(m21r_boot[, 5], .95),
                                                                quantile(m21r_boot[, 6], .05),
                                                                quantile(m21r_boot[, 7], .95))
                                  est_df[21, c(7, 12:13)] <- c(mean(m21_boot[, 4]),
                                                               quantile(m21_boot[, 9], .05),
                                                               quantile(m21_boot[, 10], .95))
                                  est_df[21, c(16, 21:22)] <- c(m21c[4, 1],
                                                                m21c[4, 1] - 1.96 * m21c[4, 2],
                                                                m21c[4, 1] + 1.96 * m21c[4, 2])
                                  est_df[21, c(14:15, 17:20)] <- c(m21rc[2:3, 1],
                                                                   m21rc[2, 1] - 1.96 * m21rc[2, 2],
                                                                   m21rc[2, 1] + 1.96 * m21rc[2, 2],
                                                                   m21rc[3, 1] - 1.96 * m21rc[3, 2],
                                                                   m21rc[3, 1] + 1.96 * m21rc[3, 2])
                      # female feeding
                              m22 <- lm(scale(d12_av_fw) ~ scale(mal_feed_intercept) * scale(mal_feed_slope), data = rem)
                              m22c <- coefficients(summary(m22))
                              m22_boot <- boot_blup(mmfeed, na.omit(ref[, c("nest_key", "d12_av_fw")]), which_re = "both_inter", boots = boot.set2)
                              m22r <- lm(scale(d12_av_fw) ~ scale(mal_feed_intercept) + scale(mal_feed_slope), data = rem) 
                              m22rc <- coefficients(summary(m22r))
                              m22r_boot <- boot_blup(mmfeed, na.omit(ref[, c("nest_key", "d12_av_fw")]), which_re = "both_single", boots = boot.set)
                              
                              # add values to table
                                  est_df[22, c(5:6, 8:11)] <- c(colMeans(m22r_boot[, 2:3]),
                                                                quantile(m22r_boot[, 4], .05),
                                                                quantile(m22r_boot[, 5], .95),
                                                                quantile(m22r_boot[, 6], .05),
                                                                quantile(m22r_boot[, 7], .95))
                                  est_df[22, c(7, 12:13)] <- c(mean(m22_boot[, 4]),
                                                               quantile(m22_boot[, 9], .05),
                                                               quantile(m22_boot[, 10], .95))
                                  est_df[22, c(16, 21:22)] <- c(m22c[4, 1],
                                                                m22c[4, 1] - 1.96 * m22c[4, 2],
                                                                m22c[4, 1] + 1.96 * m22c[4, 2])
                                  est_df[22, c(14:15, 17:20)] <- c(m22rc[2:3, 1],
                                                                   m22rc[2, 1] - 1.96 * m22rc[2, 2],
                                                                   m22rc[2, 1] + 1.96 * m22rc[2, 2],
                                                                   m22rc[3, 1] - 1.96 * m22rc[3, 2],
                                                                   m22rc[3, 1] + 1.96 * m22rc[3, 2])
                              
                              
          ## Combine data frame for plotting
                est_df2 <- data.frame(pred = c(rep("intercept", 26), rep("slope", 26)),
                                      est = c(est_df$int_est, est_df$slp_est),
                                      low = c(est_df$int_low, est_df$slp_low),
                                      hi = c(est_df$int_hi, est_df$slp_hi),
                                      est_raw = c(est_df$raw_int_est, est_df$raw_slp_est),
                                      low_raw = c(est_df$raw_int_low, est_df$raw_slp_low),
                                      hi_raw = c(est_df$raw_int_hi, est_df$raw_slp_hi),
                                      resp = rep(est_df$response, 2),
                                      model = rep(est_df$model, 2),
                                      type = rep(est_df$type, 2),
                                      group = rep(est_df$group, 2),
                                      y_pos = c(3.0-.3, 2.95-.3, 
                                                2.1, 2.05, 1.95, 1.9, 
                                                1.1+.3, 1.05+.3, .95+.3, .9+.3, 
                                                3.1-.3, 3.05-.3, 2.95-.3, 2.9-.3, 
                                                2.1, 2.05, 1.95, 1.9, 
                                                1.1+.3, 1.05+.3, .95+.3, .9+.3,
                                                .7, .65, .55, .5),
                                      rs_type = rep(c(rep("survive", 10), rep("morph", 12), rep("survive", 4)), 2)
                )
                est_df2$sig_raw <- "no"
                est_df2$sig_boot <- "no"
                for(i in 1:nrow(est_df2)){
                  if(est_df2$rs_type[i] == "survive"){
                    if(est_df2$low[i] > 1){est_df2$sig_boot[i] <- "yes"}
                    if(est_df2$hi[i] < 1){est_df2$sig_boot[i] <- "yes"}
                    if(est_df2$low_raw[i] > 1){est_df2$sig_raw[i] <- "yes"}
                    if(est_df2$hi_raw[i] < 1){est_df2$sig_raw[i] <- "yes"}
                  }
                  if(est_df2$rs_type[i] == "morph"){
                    if(est_df2$low[i] > 0){est_df2$sig_boot[i] <- "yes"}
                    if(est_df2$hi[i] < 0){est_df2$sig_boot[i] <- "yes"}
                    if(est_df2$low_raw[i] > 0){est_df2$sig_raw[i] <- "yes"}
                    if(est_df2$hi_raw[i] < 0){est_df2$sig_raw[i] <- "yes"}
                  }
                }
                
                est_df2$fill2 <- "non_sig"
                for(i in 1:nrow(est_df2)){
                  if(est_df2$sig_boot[i] == "yes"){est_df2$fill2[i] <- est_df2$group[i]}
                }
                
                survive_p <- est_df2 %>%
                  filter(rs_type == "survive") %>%
                  ggplot(mapping = aes(x = log(est_raw), y = y_pos, fill = fill2, color = group, shape = group)) +
                  geom_vline(xintercept = 0, linetype = "dashed", color = "gray50") +
                  geom_segment(aes(x = log(low), xend = log(hi), y = y_pos, yend = y_pos), linewidth = .6, color = "gray50") +
                  geom_segment(aes(x = log(low_raw), xend = log(hi_raw), y = y_pos, yend = y_pos), linewidth = 1.2) +
                  geom_point(size = 2.5, color = "black") +
                  scale_shape_manual(values = c(22, 22, 21, 21)) +
                  facet_wrap(~ pred, scale = "free") +
                  scale_fill_manual(values = c(lighten(fem_color, .1), mal_color, "white", darken("coral3", .1))) +
                  scale_color_manual(values = c(lighten(fem_color, .1), mal_color, "slateblue", darken("coral3", .1))) +
                  scale_linetype_manual(values = c("dotted", "solid")) +
                  theme_classic() +
                  #scale_y_continuous(labels = c("Recruit", "Fledge", "Day12", "Hatch"), breaks = c(0.6, 1.3, 2, 2.7)) +
                  coord_cartesian(xlim = c(-1.8, 1.8), ylim = c(.4, 3)) +
                  theme(strip.text = element_blank(), panel.grid = element_blank(), axis.text = element_text(size = 14),
                       axis.title = element_text(size = 16), axis.ticks.y = element_blank(),
                       axis.line.x.top = element_blank(), axis.line.y = element_blank(), axis.text.y = element_blank()) +
                  guides(fill = "none", color = "none", shape = "none") +
                  labs(x = "Log Odds Ratio", y = "")
                
                ggsave(here::here("output_plots", "survive_plot.svg"), survive_p, device = "svg", width = 4, height = 3.3, units = "in")
                
                grow_p <- est_df2 %>%
                  filter(rs_type == "morph") %>%
                  ggplot(mapping = aes(x = est_raw, y = y_pos, fill = fill2, color = group, shape = group)) +
                  geom_vline(xintercept = 0, linetype = "dashed", color = "gray50") +
                  geom_segment(aes(x = low, xend = hi, y = y_pos, yend = y_pos), linewidth = .6, color = "gray50") +
                  geom_segment(aes(x = low_raw, xend = hi_raw, y = y_pos, yend = y_pos), linewidth = 1.2) +
                  geom_point(size = 2.5, color = "black") +
                  scale_shape_manual(values = c(22, 22, 21, 21)) +
                  facet_wrap(~ pred, scale = "free") +
                  scale_fill_manual(values = c(lighten(fem_color, .2), "white", mal_color, "slateblue", darken("coral3", .2))) +
                  scale_color_manual(values = c(lighten(fem_color, .2), mal_color, "slateblue", darken("coral3", .2))) +
                  theme_classic() +
                  #scale_y_continuous(labels = c("Wing", "Head", "Mass"), breaks = c(1.3, 2, 2.7)) +
                  coord_cartesian(xlim = c(-.75, .75), ylim = c(1, 3)) +
                  scale_x_continuous(breaks = seq(-2, 2, 0.5)) +
                  theme(strip.text = element_blank(), panel.grid = element_blank(), axis.text = element_text(size = 14),
                        axis.title = element_text(size = 16), axis.ticks.y = element_blank(),
                        axis.line.x.top = element_blank(), axis.line.y = element_blank(), axis.text.y = element_blank()) +
                  guides(fill = "none", color = "none", shape = "none") +
                  labs(x = "Standardized Effect Size", y = "")
                
                ggsave(here::here("output_plots", "growth_plot.svg"), grow_p, device = "svg", width = 4, height = 3, units = "in")
                  
                              
                
                
        # Model predicted plots for nestling day 12 morphology (plots with data points in supplement)
                
                pmi <- as.data.frame(ggeffect(m13r, terms = "fem_feed_intercept [-3:3 by=0.1]"))
                pmi$predicted <- pmi$predicted * sd(na.omit(ref$d12_av_mass)) + mean(na.omit(ref$d12_av_mass)) 
                pmi$conf.low <- pmi$conf.low * sd(na.omit(ref$d12_av_mass)) + mean(na.omit(ref$d12_av_mass)) 
                pmi$conf.high <- pmi$conf.high * sd(na.omit(ref$d12_av_mass)) + mean(na.omit(ref$d12_av_mass)) 
                
                pms <- as.data.frame(ggeffect(m13r, terms = "fem_feed_slope [-3:3 by=0.1]"))
                pms$predicted <- pms$predicted * sd(na.omit(ref$d12_av_mass)) + mean(na.omit(ref$d12_av_mass)) 
                pms$conf.low <- pms$conf.low * sd(na.omit(ref$d12_av_mass)) + mean(na.omit(ref$d12_av_mass)) 
                pms$conf.high <- pms$conf.high * sd(na.omit(ref$d12_av_mass)) + mean(na.omit(ref$d12_av_mass)) 
                
                phi <- as.data.frame(ggeffect(m17r, terms = "fem_feed_intercept [-3:3 by=0.1]"))
                phi$predicted <- phi$predicted * sd(na.omit(ref$d12_av_hb)) + mean(na.omit(ref$d12_av_hb)) 
                phi$conf.low <- phi$conf.low * sd(na.omit(ref$d12_av_hb)) + mean(na.omit(ref$d12_av_hb)) 
                phi$conf.high <- phi$conf.high * sd(na.omit(ref$d12_av_hb)) + mean(na.omit(ref$d12_av_hb)) 
                
                phs <- as.data.frame(ggeffect(m17r, terms = "fem_feed_slope [-3:3 by=0.1]"))
                phs$predicted <- phs$predicted * sd(na.omit(ref$d12_av_hb)) + mean(na.omit(ref$d12_av_hb)) 
                phs$conf.low <- phs$conf.low * sd(na.omit(ref$d12_av_hb)) + mean(na.omit(ref$d12_av_hb)) 
                phs$conf.high <- phs$conf.high * sd(na.omit(ref$d12_av_hb)) + mean(na.omit(ref$d12_av_hb))
                
                pfi <- as.data.frame(ggeffect(m21r, terms = "fem_feed_intercept [-3:3 by=0.1]"))
                pfi$predicted <- pfi$predicted * sd(na.omit(ref$d12_av_fw)) + mean(na.omit(ref$d12_av_fw)) 
                pfi$conf.low <- pfi$conf.low * sd(na.omit(ref$d12_av_fw)) + mean(na.omit(ref$d12_av_fw)) 
                pfi$conf.high <- pfi$conf.high * sd(na.omit(ref$d12_av_fw)) + mean(na.omit(ref$d12_av_fw)) 
                
                pfs <- as.data.frame(ggeffect(m21r, terms = "fem_feed_slope [-3:3 by=0.1]"))
                pfs$predicted <- pfs$predicted * sd(na.omit(ref$d12_av_fw)) + mean(na.omit(ref$d12_av_fw)) 
                pfs$conf.low <- pfs$conf.low * sd(na.omit(ref$d12_av_fw)) + mean(na.omit(ref$d12_av_fw)) 
                pfs$conf.high <- pfs$conf.high * sd(na.omit(ref$d12_av_fw)) + mean(na.omit(ref$d12_av_fw))
                
                
                pred_mass <- ggplot(pmi, aes(x = x, y = predicted)) +
                  geom_line(color = mg_color) +
                  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = 0.4, fill = mg_color) +
                  geom_line(data = pms, color = mg_color, linetype = "dashed") +
                  geom_ribbon(data = pms, aes(ymin = conf.low, ymax = conf.high), alpha = 0.2, fill = mg_color) +
                  theme_classic() +
                  scale_x_continuous(breaks = seq(-6, 6, 2)) +
                  coord_cartesian(xlim = c(-2.8, 2.8), ylim = c(15, 23)) +
                  theme(axis.text = element_text(size = 12), axis.title = element_text(size = 14)) +
                  labs(x = "BLUP (sd units)", y = "Mass (grams)")
                
                pred_hb <- ggplot(phi, aes(x = x, y = predicted)) +
                  geom_line(color = hb_color) +
                  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = 0.4, fill = hb_color) +
                  geom_line(data = phs, color = hb_color, linetype = "dashed") +
                  geom_ribbon(data = phs, aes(ymin = conf.low, ymax = conf.high), alpha = 0.2, fill = hb_color) +
                  theme_classic() +
                  scale_x_continuous(breaks = seq(-6, 6, 2)) +
                  scale_y_continuous(breaks = seq(19, 29, 1)) +
                  coord_cartesian(xlim = c(-2.8, 2.8), ylim = c(23.2, 27)) +
                  theme(axis.text = element_text(size = 12), axis.title = element_text(size = 14)) +
                  labs(x = "BLUP (sd units)", y = "Head (mm)")
                
                pred_fw <- ggplot(pfi, aes(x = x, y = predicted)) +
                  geom_line(color = fw_color) +
                  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = 0.4, fill = fw_color) +
                  geom_line(data = pfs, color = fw_color, linetype = "dashed") +
                  geom_ribbon(data = pfs, aes(ymin = conf.low, ymax = conf.high), alpha = 0.2, fill = fw_color) +
                  theme_classic() +
                  scale_x_continuous(breaks = seq(-6, 6, 2)) +
                  scale_y_continuous(breaks = seq(10, 70, 10)) +
                  coord_cartesian(xlim = c(-2.8, 2.8), ylim = c(35, 60)) +
                  theme(axis.text = element_text(size = 12), axis.title = element_text(size = 14)) +
                  labs(x = "BLUP (sd units)", y = "Wing (mm)")
                
                morph_plot <- ggpubr::ggarrange(pred_mass, pred_hb, pred_fw, nrow = 1)
                ggsave(here::here("output_plots", "morph_plot.svg"), morph_plot, device = "svg", width = 4, height = 1.9, units = "in")
                
                
                
        # Point estimates plot for survival
                surv_ests <- data.frame(
                  response = c("1hatch", "1hatch", rep("2day12", 6), rep("3fledge", 4), 
                               "1hatch", "1hatch", "2day12", "2day12",
                               rep("4recruit", 4)),
                  z_blup = rep(c(-1, 1), 10),
                  est = NA,
                  low = NA,
                  high = NA,
                  group = c(rep("on", 4), "fem", "fem", "mal", "mal", "fem", "fem", "mal", "mal",
                            rep("off", 4), "fem", "fem", "mal", "mal"),
                  dodge = c(-.2, -.2, 
                            -.6, -.6, .2, .2, .6, .6, 
                            -.2, -.2, .2, .2, 
                            .2, .2, 
                            -.2, -.2, 
                            .2, .2, -.2, -.2)
                )
                
              
                  
                  sd.set <- 2
                  
                  
                  surv_ests[1, 3] <- inv_logit(m1rc[1,1] + -sd.set * m1rc[2,1])
                  surv_ests[1, 4] <- inv_logit(m1rc[1,1] + -sd.set * m1rc[2,1] - 1.96 * m1rc[2,2])
                  surv_ests[1, 5] <- inv_logit(m1rc[1,1] + -sd.set * m1rc[2,1] + 1.96 * m1rc[2,2])
                  surv_ests[2, 3] <- inv_logit(m1rc[1,1] + sd.set * m1rc[2,1])
                  surv_ests[2, 4] <- inv_logit(m1rc[1,1] + sd.set * m1rc[2,1] - 1.96 * m1rc[2,2])
                  surv_ests[2, 5] <- inv_logit(m1rc[1,1] + sd.set * m1rc[2,1] + 1.96 * m1rc[2,2])
                  surv_ests[3, 3] <- inv_logit(m3bbrc[1,1] + -sd.set * m3bbrc[2,1])
                  surv_ests[3, 4] <- inv_logit(m3bbrc[1,1] + -sd.set * m3bbrc[2,1] - 1.96 * m3bbrc[2,2])
                  surv_ests[3, 5] <- inv_logit(m3bbrc[1,1] + -sd.set * m3bbrc[2,1] + 1.96 * m3bbrc[2,2])
                  surv_ests[4, 3] <- inv_logit(m3bbrc[1,1] + sd.set * m3bbrc[2,1])
                  surv_ests[4, 4] <- inv_logit(m3bbrc[1,1] + sd.set * m3bbrc[2,1] - 1.96 * m3bbrc[2,2])
                  surv_ests[4, 5] <- inv_logit(m3bbrc[1,1] + sd.set * m3bbrc[2,1] + 1.96 * m3bbrc[2,2])
                  surv_ests[5, 3] <- inv_logit(m5bbrc[1,1] + -sd.set * m5bbrc[2,1])
                  surv_ests[5, 4] <- inv_logit(m5bbrc[1,1] + -sd.set * m5bbrc[2,1] - 1.96 * m5bbrc[2,2])
                  surv_ests[5, 5] <- inv_logit(m5bbrc[1,1] + -sd.set * m5bbrc[2,1] + 1.96 * m5bbrc[2,2])
                  surv_ests[6, 3] <- inv_logit(m5bbrc[1,1] + sd.set * m5bbrc[2,1])
                  surv_ests[6, 4] <- inv_logit(m5bbrc[1,1] + sd.set * m5bbrc[2,1] - 1.96 * m5bbrc[2,2])
                  surv_ests[6, 5] <- inv_logit(m5bbrc[1,1] + sd.set * m5bbrc[2,1] + 1.96 * m5bbrc[2,2])
                  surv_ests[7, 3] <- inv_logit(m6bbrc[1,1] + -sd.set * m6bbrc[2,1])
                  surv_ests[7, 4] <- inv_logit(m6bbrc[1,1] + -sd.set * m6bbrc[2,1] - 1.96 * m6bbrc[2,2])
                  surv_ests[7, 5] <- inv_logit(m6bbrc[1,1] + -sd.set * m6bbrc[2,1] + 1.96 * m6bbrc[2,2])
                  surv_ests[8, 3] <- inv_logit(m6bbrc[1,1] + sd.set * m6bbrc[2,1])
                  surv_ests[8, 4] <- inv_logit(m6bbrc[1,1] + sd.set * m6bbrc[2,1] - 1.96 * m6bbrc[2,2])
                  surv_ests[8, 5] <- inv_logit(m6bbrc[1,1] + sd.set * m6bbrc[2,1] + 1.96 * m6bbrc[2,2])
                  surv_ests[9, 3] <- inv_logit(m9bbrc[1,1] + -sd.set * m9bbrc[2,1])
                  surv_ests[9, 4] <- inv_logit(m9bbrc[1,1] + -sd.set * m9bbrc[2,1] - 1.96 * m9bbrc[2,2])
                  surv_ests[9, 5] <- inv_logit(m9bbrc[1,1] + -sd.set * m9bbrc[2,1] + 1.96 * m9bbrc[2,2])
                  surv_ests[10, 3] <- inv_logit(m9bbrc[1,1] + sd.set * m9bbrc[2,1])
                  surv_ests[10, 4] <- inv_logit(m9bbrc[1,1] + sd.set * m9bbrc[2,1] - 1.96 * m9bbrc[2,2])
                  surv_ests[10, 5] <- inv_logit(m9bbrc[1,1] + sd.set * m9bbrc[2,1] + 1.96 * m9bbrc[2,2])
                  surv_ests[11, 3] <- inv_logit(m10bbrc[1,1] + -sd.set * m10bbrc[2,1])
                  surv_ests[11, 4] <- inv_logit(m10bbrc[1,1] + -sd.set * m10bbrc[2,1] - 1.96 * m10bbrc[2,2])
                  surv_ests[11, 5] <- inv_logit(m10bbrc[1,1] + -sd.set * m10bbrc[2,1] + 1.96 * m10bbrc[2,2])
                  surv_ests[12, 3] <- inv_logit(m10bbrc[1,1] + sd.set * m10bbrc[2,1])
                  surv_ests[12, 4] <- inv_logit(m10bbrc[1,1] + sd.set * m10bbrc[2,1] - 1.96 * m10bbrc[2,2])
                  surv_ests[12, 5] <- inv_logit(m10bbrc[1,1] + sd.set * m10bbrc[2,1] + 1.96 * m10bbrc[2,2])
                  surv_ests[13, 3] <- inv_logit(m2rc[1,1] + -sd.set * m2rc[2,1])
                  surv_ests[13, 4] <- inv_logit(m2rc[1,1] + -sd.set * m2rc[2,1] - 1.96 * m2rc[2,2])
                  surv_ests[13, 5] <- inv_logit(m2rc[1,1] + -sd.set * m2rc[2,1] + 1.96 * m2rc[2,2])
                  surv_ests[14, 3] <- inv_logit(m2rc[1,1] + sd.set * m2rc[2,1])
                  surv_ests[14, 4] <- inv_logit(m2rc[1,1] + sd.set * m2rc[2,1] - 1.96 * m2rc[2,2])
                  surv_ests[14, 5] <- inv_logit(m2rc[1,1] + sd.set * m2rc[2,1] + 1.96 * m2rc[2,2])
                  surv_ests[15, 3] <- inv_logit(m4bbrc[1,1] + -sd.set * m4bbrc[2,1])
                  surv_ests[15, 4] <- inv_logit(m4bbrc[1,1] + -sd.set * m4bbrc[2,1] - 1.96 * m4bbrc[2,2])
                  surv_ests[15, 5] <- inv_logit(m4bbrc[1,1] + -sd.set * m4bbrc[2,1] + 1.96 * m4bbrc[2,2])
                  surv_ests[16, 3] <- inv_logit(m4bbrc[1,1] + sd.set * m4bbrc[2,1])
                  surv_ests[16, 4] <- inv_logit(m4bbrc[1,1] + sd.set * m4bbrc[2,1] - 1.96 * m4bbrc[2,2])
                  surv_ests[16, 5] <- inv_logit(m4bbrc[1,1] + sd.set * m4bbrc[2,1] + 1.96 * m4bbrc[2,2])
                  surv_ests[17, 3] <- inv_logit(m25rc[1,1] - sd.set * m25rc[2,1])
                  surv_ests[17, 4] <- inv_logit(m25rc[1,1] - sd.set * m25rc[2,1] - 1.96 * m25rc[2,2])
                  surv_ests[17, 5] <- inv_logit(m25rc[1,1] - sd.set * m25rc[2,1] + 1.96 * m25rc[2,2])
                  surv_ests[18, 3] <- inv_logit(m25rc[1,1] + sd.set * m25rc[2,1])
                  surv_ests[18, 4] <- inv_logit(m25rc[1,1] + sd.set * m25rc[2,1] - 1.96 * m25rc[2,2])
                  surv_ests[18, 5] <- inv_logit(m25rc[1,1] + sd.set * m25rc[2,1] + 1.96 * m25rc[2,2])
                  surv_ests[19, 3] <- inv_logit(m26rc[1,1] - sd.set * m26rc[2,1])
                  surv_ests[19, 4] <- inv_logit(m26rc[1,1] - sd.set * m26rc[2,1] - 1.96 * m26rc[2,2])
                  surv_ests[19, 5] <- inv_logit(m26rc[1,1] - sd.set * m26rc[2,1] + 1.96 * m26rc[2,2])
                  surv_ests[20, 3] <- inv_logit(m26rc[1,1] + sd.set * m26rc[2,1])
                  surv_ests[20, 4] <- inv_logit(m26rc[1,1] + sd.set * m26rc[2,1] - 1.96 * m26rc[2,2])
                  surv_ests[20, 5] <- inv_logit(m26rc[1,1] + sd.set * m26rc[2,1] + 1.96 * m26rc[2,2])
                  
            surv_ests$group2 <- surv_ests$group
            surv_ests$group2[19:20] <- "off"
                  
          surv_plot <- ggplot(surv_ests, mapping = aes(x = z_blup + dodge, y = est*100, fill = group2, shape = group)) +
            geom_segment(aes(xend = z_blup + dodge, y = low*100, yend = high*100, color = group), linewidth = 1) +
            geom_point(size = 1.7) +
            scale_shape_manual(values = c(22, 22, 21, 21)) +
            facet_wrap(~ response, scales = "free_y", nrow = 1) +
            scale_fill_manual(values = c(lighten(fem_color, .1), mal_color, "white", darken("coral3", .1))) +
            scale_color_manual(values = c(lighten(fem_color, .1), mal_color, "slateblue", darken("coral3", .1))) +
            theme_classic() +
            scale_x_continuous(breaks = c(-5, -1, 1, 5), limits = c(-2.25, 2.25)) +
            #scale_y_continuous(breaks = c(0, 5, 10, 20, 40, 60, 70, 80, 90, 100), limits = c(0, 14)) +
            theme(panel.grid = element_blank(), axis.text.x = element_blank(),
                  strip.text = element_blank(), axis.text.y = element_text(size = 12), 
                  axis.title = element_blank()) +
            guides(color = "none", shape = "none", fill = "none") +
            labs(x = "BLUP intercept", y = "Predicted Survival (%)")
          
          ggsave(here::here("output_plots", "surv_plot.svg"), surv_plot, device = "svg", width = 4, height = 1.5, units = "in")

# 
          
          ref_t <- ref
          ref_t$band_day <- ref_t$hatch_doy + 12
          ref_t$year_capday <- paste(ref_t$exp_year, ref_t$band_day, sep = "_")
          temp_jn <- yr_doy_uni1[, 4:12]
          ref_t <- plyr::join(ref_t, temp_jn, "year_capday", "left", "first")
          
          mtf1 <- lm(scale(d12_av_hb) ~ scale(fem_feed_intercept) + scale(prior5_C)*scale(fem_feed_slope), data = ref_t)
          
          ref_t$bsn <- as.numeric(ref_t$brood_size_hatching)
          ref_t$died_d12 <- ref_t$bsn - ref_t$num_d12
          ref_t$died_fledging <- ref_t$bsn - ref_t$num_fledged
          ref_t <- ref_t[is.na(ref_t$died_fledging) == FALSE & is.na(ref_t$num_fledged) == FALSE,]
          ref_t[ref_t$died_d12 < 0, "died_d12"] <- 0
          ref_t[ref_t$died_fledging < 0, "died_fledging"] <- 0
          
          mtfs <- glmer(cbind(num_fledged, died_fledging) ~ scale(fem_feed_intercept) + scale(fem_feed_slope) + (1|nest_key), data = ref_t,
                        family = binomial)
          mtfs2 <- glmmTMB(cbind(num_fledged, died_fledging) ~ scale(fem_feed_intercept)* scale(prior5_C) + scale(fem_feed_slope) + (1|nest_key),
                           data = ref_t, family = betabinomial)
          
# Between year repeatability ----
    
      # not a big enough sample size to look at incubation bout repeatability    
    
      # # on and off bouts
      #         # make a basic model to partition variance
      #               all_bouts <- all_bouts[, !duplicated(names(all_bouts))]
      #               ab_fem2 <- all_bouts 
      #               
      #           # for on bouts
      #               # filter and scale
      #                 xab_fon2 <- ab_fem2 %>%
      #                   filter(is.na(temp_C) == FALSE)
      #                 xab_fon2$temp_s <- scale(xab_fon2$temp_C) 
      #                 xab_fon2$duration_s <- scale(xab_fon2$duration/60)
      #                 xab_fon2$duration_m <- xab_fon2$duration/60
      #                 xab_fon2 <- xab_fon2[xab_fon2$duration_m > 3, ]
      #               
      #               # fit a simple model
      #               # 
      #                 xab_fon12 <- xab_fon2[xab_fon2$type == "On" & xab_fon2$temp_C < 19, ]
      #                 mi_on2 <- lmer(duration_m ~ temp_s + (temp_s|nest_key), data = xab_fon12)
      #                 temp_mean2 <- mean(xab_fon2$temp_C)
      #                 temp_sd2 <- sd(xab_fon2$temp_C)
      #               
      #               # exctract random effects and convert back to measurement scale  
      #                 re_mion2 <- ranef(mi_on2)$nest_key
      #                 colnames(re_mion2) <- c("r_intercept", "r_slope")
      #                 re_mion2$r_intercept <- fixef(mi_on2)[1] + re_mion2$r_intercept - (fixef(mi_on2)[2] * temp_mean2 / temp_sd2)
      #                 re_mion2$r_slope <- (fixef(mi_on2)[2] + re_mion2$r_slope) / temp_sd2
      #                 re_mion2$nest_key <- rownames(re_mion2)
      #                 re_mion2 <- plyr::join(re_mion2, nests[, c("nest_key", "female_id", "exp_year")], "nest_key", "left", "first")
      #                 for(i in 1:nrow(re_mion2)){
      #                   re_mion2$num_obs[i] <- nrow(subset(re_mion2, re_mion2$female_id == re_mion2$female_id[i]))
      #                 }
      #                 on_rpt <- re_mion2[re_mion2$num_obs > 1, ]
      #                 
                       ages <- ith_capture[, c("band", "exp_year", "age")]
                       colnames(ages)[1] <- "female_id"
      #                 on_rpt <- plyr::join(on_rpt, ages, c("female_id", "exp_year"), "left", "first")
      #                 
      #                 on_rep1 <- rpt(scale(r_intercept) ~ age + scale(r_slope) + (1|female_id), grname = "female_id", data = on_rpt, datatype = "Gaussian")
      #                 on_rep2 <- rpt(scale(r_slope) ~ age + scale(r_intercept) + (1|female_id), grname = "female_id", data = on_rpt, datatype = "Gaussian")
      #                 
      #                 anova(lm(r_slope ~ as.factor(female_id), data = on_rpt))
      #                  
      #                 
      #           # for off bouts
      #                 
      #               # fit a simple model
      #                 xab_fon2b <- xab_fon2[xab_fon2$type == "Off" & xab_fon2$temp_C < 19, ]
      #                 mi_off2 <- lmer(duration_m ~ temp_s + (temp_s|nest_key), data = xab_fon2b)
      #                 
      #               # exctract random effects and convert back to measurement scale  
      #                 re_mioff2 <- ranef(mi_off2)$nest_key
      #                 colnames(re_mioff2) <- c("r_intercept", "r_slope")
      #                 re_mioff2$r_intercept <- fixef(mi_off2)[1] + re_mioff2$r_intercept - (fixef(mi_off2)[2] * temp_mean2 / temp_sd2)
      #                 re_mioff2$r_slope <- (fixef(mi_off2)[2] + re_mioff2$r_slope) / temp_sd2
      #                 re_mioff2$nest_key <- rownames(re_mioff2)
      #                 re_mioff2 <- plyr::join(re_mioff2, nests[, c("nest_key", "female_id", "exp_year")], "nest_key", "left", "first")
      #                 for(i in 1:nrow(re_mioff2)){
      #                   re_mioff2$num_obs[i] <- nrow(subset(re_mioff2, re_mioff2$female_id == re_mioff2$female_id[i]))
      #                 }
      #                 off_rpt <- re_mioff2[re_mioff2$num_obs > 1, ]
      #                 off_rpt <- plyr::join(off_rpt, ages, c("female_id", "exp_year"), "left", "first")
      #                 
      #                 off_rep1 <- rpt(scale(r_intercept) ~ age + scale(r_slope) + (1|female_id), grname = "female_id", data = off_rpt, datatype = "Gaussian")
      #                 off_rep2 <- rpt(scale(r_slope) ~ age + scale(r_intercept) + (1|female_id), grname = "female_id", data = off_rpt, datatype = "Gaussian")
      #                 
      #                 anova(lm(r_intercept ~ as.factor(female_id), data = off_rpt))
                
                
                
      # feeding rate          
          # starting with data frames from above
                        af_fem2 <- all_feedsf %>% filter(female_capture_day == "no", f_feeds < 37)
                        af_mal2 <- all_feedsm %>% filter(male_capture_day == "no", m_feeds < 37)
                    
                # for feeding rate
                    # filter and scale
                      #females
                        af_ffeed2 <- af_fem2 %>%
                          filter(temp_C < 25, is.na(temp_C) == FALSE)
                        af_ffeed2$temp_s <- scale(af_ffeed2$temp_C) 
                        af_ffeed2$hour_s <- scale(af_ffeed2$Hour) 
                        af_ffeed2$age_s <- scale(af_ffeed2$Offset)
                        af_ffeed2$nest_key <- as.factor(af_ffeed2$nest_key)
                      #males
                        af_mfeed2 <- af_mal2 %>%
                          filter(temp_C < 25, is.na(temp_C) == FALSE)
                        af_mfeed2$temp_s <- scale(af_mfeed2$temp_C) 
                        af_mfeed2$hour_s <- scale(af_mfeed2$Hour) 
                        af_mfeed2$age_s <- scale(af_mfeed2$Offset)
                    
                    # fit a simple model as basis for extracting random effects
                        # more complex models with other variables were tried but the random effect estimates
                        # are very similar and this is easier to understand
                      # females
                        mffeed2 <- lmer(f_feeds ~ temp_s + (temp_s|nest_key), data = af_ffeed2)
                        temp_meanf2 <- mean(af_ffeed2$temp_C)
                        temp_sdf2 <- sd(af_ffeed2$temp_C)
                        
                        # inspect model. with huge sample size, tests are significant but visual residuals look fine
                          # sf <- simulateResiduals(mffeed)
                          # plot(sf)
                          # testDispersion(sf)
                          # hist(residuals(mffeed))
                        
                      # males
                        mmfeed2 <- lmer(m_feeds ~ temp_s + (temp_s|nest_key), data = af_mfeed2, control = lmerControl(optimizer = "bobyqa"))
                        temp_meanm2 <- mean(af_mfeed2$temp_C)
                        temp_sdm2 <- sd(af_mfeed2$temp_C)
                        
                        # inspect model. with huge sample size, tests are significant but visual residuals look fine
                          # sf <- simulateResiduals(mmfeed)
                          # plot(sf)
                          # testDispersion(sf)
                          # hist(residuals(mmfeed))
                    
                    # extract random effects and convert back to measurement scale  
                        # females
                              re_ffeed2 <- ranef(mffeed2)$nest_key  
                              colnames(re_ffeed2) <- c("r_intercept", "r_slope")
                              re_ffeed2$r_intercept <- fixef(mffeed2)[1] + re_ffeed2$r_intercept - (fixef(mffeed2)[2] * temp_meanf2 / temp_sdf2)
                              re_ffeed2$r_slope <- (fixef(mffeed2)[2] + re_ffeed2$r_slope) / temp_sdf2
                              re_ffeed2$nest_key <- rownames(re_ffeed2)
                              re_ffeed2 <- plyr::join(re_ffeed2, nests[, c("nest_key", "female_id", "exp_year")], "nest_key", "left", "first")
                              for(i in 1:nrow(re_ffeed2)){
                                re_ffeed2$num_obs[i] <- nrow(subset(re_ffeed2, re_ffeed2$female_id == re_ffeed2$female_id[i]))
                              }
                              
                               ref_rpt <- re_ffeed2[re_ffeed2$num_obs > 1, ]
                               ref_rpt <- plyr::join(ref_rpt, ages, c("female_id", "exp_year"), "left", "first")
                              # 
                               rep1 <- rpt(scale(r_intercept) ~ age + scale(r_slope) + (1 | female_id), grname = "female_id", data = ref_rpt, datatype = "Gaussian",
                                   nboot = 2000)
                               rep2 <- rpt(scale(r_slope) ~ age + scale(r_intercept) + (1 | female_id), grname = "female_id", data = ref_rpt, datatype = "Gaussian",
                                   nboot = 500)
                               ref_rpt$band <- ref_rpt$female_id
                          
                        # males
                              re_mfeed2 <- ranef(mmfeed2)$nest_key  
                              colnames(re_mfeed2) <- c("r_intercept", "r_slope")
                              re_mfeed2$r_intercept <- fixef(mmfeed2)[1] + re_mfeed2$r_intercept - (fixef(mmfeed2)[2] * temp_meanm2 / temp_sdm2)
                              re_mfeed2$r_slope <- (fixef(mmfeed2)[2] + re_mfeed2$r_slope) / temp_sdm2
                              re_mfeed2$nest_key <- rownames(re_mfeed2)
                              re_mfeed2 <- plyr::join(re_mfeed2, nests[, c("nest_key", "male_id", "exp_year")], "nest_key", "left", "first")
                              for(i in 1:nrow(re_mfeed2)){
                                re_mfeed2$num_obs[i] <- nrow(subset(re_mfeed2, re_mfeed2$male_id == re_mfeed2$male_id[i]))
                              }
                              
                            
                              rem_rpt <- re_mfeed2[re_mfeed2$num_obs > 1, ]
                              
                              repm1 <- rpt(scale(r_intercept) ~ scale(r_slope) + (1 | male_id), grname = "male_id", data = rem_rpt, datatype = "Gaussian",
                                           nboot = 2000)
                              repm2 <- rpt(scale(r_slope) ~ scale(r_intercept) + (1 | male_id), grname = "male_id", data = rem_rpt, datatype = "Gaussian",
                                           nboot = 2000)
                               rem_rpt$band <- rem_rpt$male_id
                              
                              comb_rpt <- bind_rows(ref_rpt, rem_rpt)
                              
                              
                              rpt_comb <- rpt(scale(r_intercept) ~ (1 | band) + (1 | exp_year), grname = "band", data = comb_rpt, datatype = "Gaussian",
                                              nboot = 2000)
                              rpt_comb2 <- rpt(scale(r_slope) ~ scale(r_intercept) + (1 | band) + (1|exp_year), grname = "band", data = comb_rpt, datatype = "Gaussian",
                                               nboot = 2000)
                              
                            # account for mate
                              # ref_m_rpt <- ref_rpt[, c("female_id", "exp_year", "nest_key", "r_intercept", "r_slope")]
                              # colnames(ref_m_rpt) <- c("female_id", "exp_year", "nest_key", "fr_intercept", "fr_slope")
                              # rem_f_rpt <- rem_rpt[, c("male_id", "exp_year", "nest_key", "r_intercept", "r_slope")]
                              # colnames(rem_f_rpt) <- c("male_id", "exp_year", "nest_key", "mr_intercept", "mr_slope")
                              # 
                              # comb2_rpt <- plyr::join(ref_m_rpt, rem_f_rpt, "nest_key", "left", "first")
                              # comb2_rpt <- na.omit(comb2_rpt[, !duplicated(names(comb2_rpt))])
                              # 
                              # rptc2 <- rpt(scale(fr_slope) ~ scale(mr_slope) + scale(mr_intercept) + scale(fr_intercept) + (1|female_id), grname = "female_id", data = comb2_rpt, datatype = "Gaussian")
                              # 
                   
# Change in feed rate and return latency ----
  # in addition to BLUPs from models above, calculate the change in feeding rate before vs
    # after capture in females (resilience) and latency to return to box after capture
        
  # latency comes directly from the processed RFID as next reading of the bird after release time
          # read in return latency
            return_latency <- read.csv(here::here("processed_data_output", "return_latencies.csv"))
            
          # subset to males  
                ret_males <- return_latency %>%
                  filter(sex == "Male", ret_lat_m < 300) %>%
                  select(nest_key, encounter_doy, ret_lat_m)
                colnames(ret_males) <- c("nest_key", "encounter_doy", "mal_retlatm")
                ret_males <- plyr::join(ret_males, nests[, c("nest_key", "hatch_doy")], "nest_key", "left", "first")
                ret_males$Offset <- ret_males$encounter_doy - ret_males$hatch_doy
                ret_males <- ret_males %>% filter(Offset > 2, !duplicated(nest_key)) %>%
                  select(nest_key, mal_retlatm)
                colnames(ret_males) <- c("nest_key", "mal_prov_latency")
            
          # subset to females (in incubation or in provisioning)
                ret_females <- return_latency %>%
                  filter(sex == "Female", ret_lat_m < 300) %>%
                  select(nest_key, encounter_doy, ret_lat_m)
                colnames(ret_females) <- c("nest_key", "encounter_doy", "fem_retlatm")
                ret_females <- plyr::join(ret_females, nests[, c("nest_key", "hatch_doy")], "nest_key", "left", "first")
                ret_females$Offset <- ret_females$encounter_doy - ret_females$hatch_doy
                ret_females$stage <- "exclude"
                ret_females <- ret_females %>% filter(is.na(Offset) == FALSE)
                for(i in 1:nrow(ret_females)){
                  if(ret_females$Offset[i] < -2){ret_females$stage[i] <- "incubation"}
                  if(ret_females$Offset[i] > 4){ret_females$stage[i] <- "provisioning"}
                }
                ret_females <- ret_females %>%
                  dplyr::group_by(nest_key, stage) %>%
                  summarise(fem_retlatm = mean(fem_retlatm, na.rm = TRUE)) %>%
                  pivot_wider(id_cols = c(nest_key), names_from = stage, values_from = fem_retlatm) %>%
                  select(nest_key, incubation, provisioning)
                colnames(ret_females) <- c("nest_key", "fem_inc_latency", "fem_prov_latency")
                
          # make descriptive figure
                comb_latency <- data.frame(lat_min = c(ret_males$mal_prov_latency, ret_females$fem_inc_latency, ret_females$fem_prov_latency) / 60,
                                           Sex = c(rep("Males", nrow(ret_males)), rep("Females", nrow(ret_females)*2)),
                                           stage = c(rep("Provisioning", nrow(ret_males)), rep("Incubation", nrow(ret_females)),
                                                     rep("Provisioning", nrow(ret_females))))
                
                return_plot <- ggplot(comb_latency, mapping = aes(x = lat_min, fill = Sex)) +
                  geom_histogram(breaks = seq(0, 5, .5), alpha = 0.6, color = "gray25") +
                  facet_wrap(~Sex + stage, scales = "free_y") +
                  scale_fill_manual(values = c(fem_color, mal_color)) +
                  theme_classic() +
                  theme(panel.grid = element_blank(), axis.title = element_text(size = 14), axis.text = element_text(size = 12)) +
                  labs(y = "Count", x = "Return latency (hours)")
                
                ggsave(here::here("output_plots", "return_plot.png"), return_plot, device = "png", width = 7.2, height = 2.6, units = "in", dpi = 300)
        
                
    # Change in feeding rate after capture. Recorded only for females during provisioning.
        # We need to do a bunch of wrangling to get the data for this
              
              # start with a data frame that has each capture day for all females    
                    cap_days <- all_feeds %>%
                      filter(female_capture_day == "yes") %>%
                      select(nest_key, female_id, year, yday, Offset)
                    cap_days <- unique(cap_days)
                
              # split this out into dataframes that are the day before and day after capture
                    dbef_capture <- cap_days[, c("nest_key", "year", "yday")]
                    dbef_capture$yday <- dbef_capture$yday - 1
                    dbef_capture$dbef_capture <- "yes"
                    
                    daft_capture <- cap_days[, c("nest_key", "year", "yday")]
                    daft_capture$yday <- daft_capture$yday + 1
                    daft_capture$daft_capture <- "yes"
                
              # now join feeding data to capture day-before/day-of/day-after
                    cap_sens <- all_feeds
                    cap_sens <- plyr::join(cap_sens, dbef_capture, c("nest_key", "year", "yday"))
                    cap_sens <- plyr::join(cap_sens, daft_capture, c("nest_key", "year", "yday"))
                    cap_sens$dbef_capture <- ifelse(is.na(cap_sens$dbef_capture) == TRUE, "no", "yes")
                    cap_sens$daft_capture <- ifelse(is.na(cap_sens$daft_capture) == TRUE, "no", "yes")
                
              # filter out to only day before, of, after and calculate feeding rate from noon to 5pm (at least a few hours after capture)
                    cap_sens2 <- cap_sens %>%
                      filter(dbef_capture == "yes" | daft_capture == "yes" | female_capture_day == "yes") %>%
                      filter(Hour > 12 & Hour < 17, Offset > 2, Offset < 10) %>%
                      dplyr::group_by(nest_key, female_capture_day, dbef_capture, daft_capture) %>%
                      summarise(avg_feeds = mean(f_feeds, na.rm = TRUE))
                
                
              # Filter out rows where multiple conditions are "yes". These have some joining error
                    cap_sens2_filtered <- cap_sens2 %>%
                      filter(
                        (female_capture_day == "yes" & dbef_capture == "no" & daft_capture == "no") |
                          (female_capture_day == "no" & dbef_capture == "yes" & daft_capture == "no") |
                          (female_capture_day == "no" & dbef_capture == "no" & daft_capture == "yes")
                      )
                
                # Pivot wider and wrangle for plotting
                      cap_sens2_wide <- cap_sens2_filtered %>%
                        pivot_wider(
                          names_from = c(female_capture_day, dbef_capture, daft_capture),
                          values_from = avg_feeds,
                          names_glue = "{female_capture_day}_{dbef_capture}_{daft_capture}_avg_feeds"
                        ) %>%
                        rename(
                          after = `no_no_yes_avg_feeds`,  # Female capture = no, DBef = no, DAft = yes
                          before = `no_yes_no_avg_feeds`,  # Female capture = no, DBef = yes, DAft = no
                          during = `yes_no_no_avg_feeds`    # Female capture = yes, DBef = no, DAft = no
                        )
                
                        #cap_sens2_wide$diff <- cap_sens2_wide$after - cap_sens2_wide$before   # this was to day after, not using
                        cap_sens2_wide$diff <- cap_sens2_wide$during - cap_sens2_wide$before   # this is to day of
                        
                # make a descriptive plot showing change in feeding
                        # cs2l <- cap_sens2_wide %>%
                        #   pivot_longer(cols = c(after, during, before), names_to = "type")
                        # type2 <- data.frame(type = c("after", "during", "before"), type2 = c("3after", "2during", "1before"))
                        # cs2l <- plyr::join(cs2l, type2, "type", "left")
                
                
                      # ggplot(cs2l, aes(x = type2, y = value)) + 
                      #   geom_line(aes(group = nest_key), alpha = 0.4, color = "gray75", linewidth = 0.4) +
                      #   geom_boxplot(fill = fem_color, alpha = 0.7, color = "gray25", width = 0.4, outlier.shape = NA) +
                      #   guides(color = "none") +
                      #   theme_classic() +
                      #   theme(panel.grid = element_blank(), axis.title = element_text(size = 14), axis.text = element_text(size = 12)) +
                      #   scale_x_discrete(labels = c("Day before", "Capture day", "Day after")) +
                      #   labs(y = "Hourly feeding rate \n (12pm - 5pm)", x = "")
                
                
        
# Compare all behavioral measures to each other ----
    # test for correlations between each pair of inc/prov blups + sensitivity to capture 
    # collecting all the measures together from each section above
          
        # start with on bout slope intercept
            behav_comp <- ranef(mi_on)$nest_key
            colnames(behav_comp) <- c("on_int", "on_slp")
            behav_comp$nest_key <- rownames(behav_comp)
            behav_comp <- behav_comp[, c("nest_key", "on_int", "on_slp")]
            
        # add off bout slope and intercept
            off_tmp <- ranef(mi_off)$nest_key
            colnames(off_tmp) <- c("off_int", "off_slp")
            off_tmp$nest_key <- rownames(off_tmp)
            
            behav_comp <- plyr::join(behav_comp, off_tmp, "nest_key", type = "full")
            
        # add female feeding slope/intercept
            ff_tmp <- ranef(mffeed)$nest_key
            colnames(ff_tmp) <- c("f_feed_int", "f_feed_slp")
            ff_tmp$nest_key <- rownames(ff_tmp)
            
            behav_comp <- plyr::join(behav_comp, ff_tmp, "nest_key", type = "full")
            
        # add male feeding slope/intercept
            mf_tmp <- ranef(mmfeed)$nest_key
            colnames(mf_tmp) <- c("m_feed_int", "m_feed_slp")
            mf_tmp$nest_key <- rownames(mf_tmp)
            
            behav_comp <- plyr::join(behav_comp, mf_tmp, "nest_key", type = "full")
            
        # add lagged on bout slope
            on_tmp_lag <- ranef(mi_on_lag)$nest_key
            colnames(on_tmp_lag) <- c("on_int_lag", "on_slp_lag")
            on_tmp_lag$nest_key <- rownames(on_tmp_lag)
            on_tmp_lag <- on_tmp_lag[, c("nest_key", "on_slp_lag")]
            
            behav_comp <- plyr::join(behav_comp, on_tmp_lag, "nest_key", type = "full")
            
        # add lagged off bout slope
            off_tmp_lag <- ranef(mi_off_lag)$nest_key
            colnames(off_tmp_lag) <- c("off_int_lag", "off_slp_lag")
            off_tmp_lag$nest_key <- rownames(off_tmp_lag)
            off_tmp_lag <- off_tmp_lag[, c("nest_key", "off_slp_lag")]
            
            behav_comp <- plyr::join(behav_comp, off_tmp_lag, "nest_key", type = "full")
            
        # add lagged feeding rate female
            ffeed_tmp_lag <- ranef(mffeed_lag)$nest_key
            colnames(ffeed_tmp_lag) <- c("f_feed_int_lag", "f_feed_slp_lag")
            ffeed_tmp_lag$nest_key <- rownames(ffeed_tmp_lag)
            ffeed_tmp_lag <- ffeed_tmp_lag[, c("nest_key", "f_feed_slp_lag")]
            
            behav_comp <- plyr::join(behav_comp, ffeed_tmp_lag, "nest_key", type = "full")
         
        # add lagged feeding rate male
            mfeed_tmp_lag <- ranef(mmfeed_lag)$nest_key
            colnames(mfeed_tmp_lag) <- c("m_feed_int_lag", "m_feed_slp_lag")
            mfeed_tmp_lag$nest_key <- rownames(mfeed_tmp_lag)
            mfeed_tmp_lag <- mfeed_tmp_lag[, c("nest_key", "m_feed_slp_lag")]
            
            behav_comp <- plyr::join(behav_comp, mfeed_tmp_lag, "nest_key", type = "full")  
            
        # add feeding difference
            behav_comp <- plyr::join(behav_comp, cap_sens2_wide[, c("nest_key", "diff")], "nest_key", type = "full")
            
        # add return latency
            # females
              behav_comp <- plyr::join(behav_comp, ret_females, "nest_key", type = "full")
            # males
              behav_comp <- plyr::join(behav_comp, ret_males, "nest_key", type = "full")
              
    # standardize and drop nest key
          behav_comp_s <- behav_comp %>% select(-nest_key)
          behav_comp_s <- as.data.frame(apply(behav_comp_s, 2, scale))
          
          behav_comp_s <- behav_comp_s %>% filter(!is.na(on_int) | !is.na(f_feed_int))
          
    # for robustness and resilience columns, sometimes positive values = more resilient and others = more sensitive
          # to facilitate comparisons, I'm converting them so that in all cases + = more/resilient/robust
          
        
          behav_comp_s$on_slp_lag <- behav_comp_s$on_slp_lag * -1
          behav_comp_s$off_slp_lag <- behav_comp_s$off_slp_lag * -1
          behav_comp_s$fem_inc_latency <- behav_comp_s$fem_inc_latency * -1
          behav_comp_s$fem_prov_latency <- behav_comp_s$fem_prov_latency * -1
          behav_comp_s$mal_prov_latency <- behav_comp_s$mal_prov_latency * -1
          behav_comp_s$diff <- behav_comp_s$diff * -1
          behav_comp_s$m_feed_slp <- behav_comp_s$m_feed_slp * -1
          behav_comp_s$m_feed_slp_lag <- behav_comp_s$m_feed_slp_lag * -1
          behav_comp_s$f_feed_slp <- behav_comp_s$f_feed_slp * -1
          behav_comp_s$f_feed_slp_lag <- behav_comp_s$f_feed_slp_lag * -1
        
    ## make a plot
        
      # get correlation  
        cor_results <- Hmisc::rcorr(as.matrix(behav_comp_s), type = "pearson")
        cor_matrix <- cor_results$r
        cor_matrix <- cor_matrix[, c(1, 2, 9, 3, 4, 10, 14, 16, 15, 13, 7, 8, 12, 5, 6, 11)]
        p_matrix <- cor_results$P
        p_matrix <- p_matrix[, c(1, 2, 9, 3, 4, 10, 14, 16, 15, 13, 7, 8, 12, 5, 6, 11)]
        
      # Melt the correlation and p-value matrices into long format for ggplot2
        cor_melt <- reshape2::melt(cor_matrix)
        p_melt <- reshape2::melt(p_matrix)  
        
      # Combine correlation and p-values into one data frame
        cor_data <- data.frame(cor_melt, p_value = p_melt$value)

      # Add significance stars
        cor_data$significance <- ifelse(cor_data$p_value <= 0.001, "***", 
                                      ifelse(cor_data$p_value <= 0.01, "**", 
                                             ifelse(cor_data$p_value <= 0.05, "*", "")))
        
      # Ensure Var1 and Var2 are characters before converting to numeric indices
        cor_data$Var1 <- as.character(cor_data$Var1)
        cor_data$Var2 <- as.character(cor_data$Var2)
        
      # Find numeric indices from column names (the variables in your correlation matrix)
        col_names <- colnames(cor_matrix)
        cor_data$Var1_index <- match(cor_data$Var1, col_names)
        cor_data$Var2_index <- match(cor_data$Var2, col_names)
        
      # Filter for lower triangle: only rows where Var1_index > Var2_index
        cor_data <- cor_data[cor_data$Var1_index > cor_data$Var2_index, ]
        
      # Plot using ggplot2
        blup_cors <- ggplot(cor_data, aes(Var1_index, Var2_index, fill = value)) + 
          geom_tile(color = "gray25") +
          geom_text(aes(label = sprintf("%.2f", value)), color = "black", size = 3) + 
          geom_text(aes(label = significance), color = "black", size = 4, vjust = 1.5) +  # Adjust vjust to move stars down
          scale_fill_gradient2(low = "coral3", high = "slateblue", mid = "white", midpoint = 0) + 
          theme_minimal() +
          scale_x_continuous(breaks = seq(1, length(col_names), by = 1), labels = col_names) +
          scale_y_continuous(breaks = seq(1, length(col_names), by = 1), labels = col_names) +
          theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
          labs(x = "", y = "") +
          theme(panel.grid = element_blank())
        
        ggsave(here::here("output_plots", "blup_cors.svg"), blup_cors, device = "svg",
               width = 8.4, height = 7, units = "in")
        
        
    # Make a network plot showing signifcant correlations between behaviors
        library(igraph)
        library(tidygraph)
        library(ggraph)
        
        # characteristics of nodes
            node_info <- data.frame(
              name = unique(c(cor_data$Var1, cor_data$Var2)),
              color_category = c("robust", "difference", "robust", "difference",
                                 "robust", "difference", "robust", rep("resilience", 8),
                                 "difference"),
              #shape_category = c(rep("incubation", 3), rep("provisioning", 4),
              #                   rep("incubation", 2), rep("provisioning", 3),
              #                   "incubation", rep("provisioning", 2), "incubation")
              shape_category = c(rep("Female", 5), "Male", "Male", rep("Female", 3),
                                 "Male", rep("Female", 3), "Male", "Female")
            )
            
            node_info$color_category <- factor(node_info$color_category)
            node_info$shape_category <- factor(node_info$shape_category)
        
        # filter correlations to only significant correlations
            cor_data_filtered <- subset(cor_data, p_value > 0.05 & p_value < 0.1)
            cor_data_filtered[nrow(cor_data_filtered) + 1, ] <-
              c("mal_prov_latency", "on_int", 0, 0.5, NA, 8, 1)
            cor_data_filtered$value <- as.numeric(cor_data_filtered$value)
            cor_data_filtered$p_value <- as.numeric(cor_data_filtered$p_value)
            
        # create edge data frame
            edges <- cor_data_filtered[, c("Var1", "Var2", "value")]
            colnames(edges) <- c("from", "to", "weight")
            #edges$weight <- abs(edges$weight)
            
        # reorder in the way I want them
            custom_order <- c("on_int", "on_slp", "on_slp_lag", 
                              "off_int", "off_slp", "off_slp_lag",
                              "fem_inc_latency", "mal_prov_latency", "fem_prov_latency",
                              "diff", "m_feed_int", "m_feed_slp",
                              "m_feed_slp_lag", "f_feed_int", "f_feed_slp", "f_feed_slp_lag")
            
            
        # create graph
            graph <- graph_from_data_frame(edges, directed = FALSE)
            
        # convert to tidygraph and add node attributes
            graph_tbl <- as_tbl_graph(graph) %>%
              activate(nodes) %>%
              left_join(node_info, by = c("name"))
            graph_tbl <- permute(graph_tbl, match(V(graph_tbl)$name, custom_order))  
            
            #graph_tbl_filtered <- delete_edges(graph_tbl, E(graph_tbl)[p_value > 0.05])
            
        # define colors and shapes
            color_palette <- c("robust" = robust_color, "difference" = diff_color, "resilience" = resil_color)
            #shape_palette <- c("incubation" = 21, "provisioning" = 22)
            shape_palette <- c("Female" = 21, "Male" = 22)
            

            
        # plot network graph
            set.seed(seed = 189)
            network_plot <- ggraph(graph_tbl, layout = "circle") +
              geom_edge_link(aes(color = weight), edge_width = 1.5) +
              geom_node_point(aes(fill = color_category, shape = shape_category), size = 7) +
              geom_node_text(aes(label = name), repel = TRUE) +
              #scale_edge_width(range = c(1, 2.5)) + 
              scale_fill_manual(values = color_palette) +
              scale_shape_manual(values = shape_palette) +
              scale_edge_color_gradient2(low = "coral3", high = "slateblue", mid = "white", midpoint = 0) +
              theme_void()
            
            ggsave(here::here("output_plots", "behav_network.svg"), network_plot,
                   device = "svg", width = 7, height = 6, units = "in", dpi = 300)
        
        
# Aggregate effects of temperature from incubation ----
  # plots of aggregate influence of bout duration on absolute egg temperature and temperature swing
  # We only report this at aggregate level because exact egg placement can have a big effect on absolute
    # egg temperature so it isn't super reliable as an individual level measure
           #models   
              all_bouts$duration_m <- all_bouts$duration / 60
              maxm <- lmer(max_temp ~ duration_m*temp_C + (1|nest_key), data = all_bouts[all_bouts$type2 == "On", ])
              minm <- lmer(min_temp ~ duration_m*temp_C + (1|nest_key), data = all_bouts[all_bouts$type2 == "Off", ])
              
              mswingon <- lmer(delta ~ duration_m*temp_C + (1|nest_key), data = all_bouts[all_bouts$type2 == "On", ])
              mswingoff <- lmer(delta ~ duration_m*temp_C + (1|nest_key), data = all_bouts[all_bouts$type2 == "Off", ]) 
              
              maxm_s <- lmer(scale(max_temp) ~ scale(duration_m)*scale(temp_C) + (1|nest_key), data = all_bouts[all_bouts$type2 == "On", ])
              minm_s <- lmer(scale(min_temp) ~ scale(duration_m)*scale(temp_C) + (1|nest_key), data = all_bouts[all_bouts$type2 == "Off", ])
              
              mswingon_s <- lmer(scale(delta) ~ scale(duration_m)*scale(temp_C) + (1|nest_key), data = all_bouts[all_bouts$type2 == "On", ])
              mswingoff_s <- lmer(scale(delta) ~ scale(duration_m)*scale(temp_C) + (1|nest_key), data = all_bouts[all_bouts$type2 == "Off", ])
             
            # plots 
              maxp <- as.data.frame(ggeffect(maxm, terms = c("duration_m [0:80 by = 0.1]", "temp_C [15, 20, 25]")))
              
              tp_1 <- ggplot(maxp, aes(x = x, y = predicted, linetype = group)) +
                geom_line(color = on_color) +
                geom_ribbon(aes(x = x, ymin = conf.low, ymax = conf.high), fill = on_color, alpha = 0.3, color = "transparent") +
                theme_classic() +
                theme(panel.grid = element_blank(), axis.text = element_text(size = 12),
                      element_text(size = 14)) +
                coord_cartesian(xlim = c(3, 35), ylim = c(32.6, 36.7)) +
                labs(x = "On Bout Duration (min)", y = "Ending Temperature (C)")
              
              minp <- as.data.frame(ggeffect(minm, terms = c("duration_m [0:80 by = 0.1]", "temp_C [15, 20, 25]")))
              
              tp_2 <- ggplot(minp, aes(x = x, y = predicted, linetype = group)) +
                geom_line(color = off_color) +
                geom_ribbon(aes(x = x, ymin = conf.low, ymax = conf.high), alpha = 0.3, color = "transparent", fill = off_color) +
                theme_classic() +
                theme(panel.grid = element_blank(), axis.text = element_text(size = 12),
                      element_text(size = 14)) +
                coord_cartesian(xlim = c(3, 35), ylim = c(26, 33)) +
                labs(x = "Off Bout Duration (min)", y = "Ending Temperature (C)")
              
              swingonp <- as.data.frame(ggeffect(mswingon, terms = c("duration_m [0:80 by = 0.1]", "temp_C [15, 20, 25]")))
              
              tp_3 <- ggplot(swingonp, aes(x = x, y = predicted, linetype = group)) +
                geom_line(color = on_color) +
                geom_ribbon(aes(x = x, ymin = conf.low, ymax = conf.high), alpha = 0.3, color = "transparent", fill = on_color) +
                theme_classic() +
                theme(panel.grid = element_blank(), axis.text = element_text(size = 12),
                      element_text(size = 14)) +
                coord_cartesian(xlim = c(3, 35), ylim = c(2.1, 5.5)) +
                labs(x = "Off Bout Duration (min)", y = "Temperature Swing (C)")
              
              swingoffp <- as.data.frame(ggeffect(mswingoff, terms = c("duration_m [0:80 by = 0.1]", "temp_C [15, 20, 25]")))
              
              tp_4 <- ggplot(swingoffp, aes(x = x, y = predicted, linetype = group)) +
                geom_line(color = off_color) +
                geom_ribbon(aes(x = x, ymin = conf.low, ymax = conf.high), alpha = 0.3, color = "transparent", fill = off_color) +
                theme_classic() +
                theme(panel.grid = element_blank(), axis.text = element_text(size = 12),
                      element_text(size = 14)) +
                coord_cartesian(xlim = c(3, 35), ylim = c(-6.5, -2)) +
                labs(x = "Off Bout Duration (min)", y = "Temperature Swing (C)")
              
              tp_all <- ggpubr::ggarrange(tp_1, tp_2, tp_3, tp_4)
              ggsave(here::here("output_plots", "bout_temperature.svg"), tp_all, device = "svg",
                     width = 7, height = 7, units = "in")
        
        
        