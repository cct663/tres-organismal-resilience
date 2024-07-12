# Description ----
    # This code is written to process HOBO files to identify on and off
    # bout lengths and characteristics. 
    # Last updated July 2024.
    # Written by Conor Taff: cct663@gmail.com
    
    # INSTRUCTIONS:
    # This code is designed to be run in a folder that has nested 
    # folders inside it. The main folder contains the R-Studio project
    # that is used to process the code. The sub folders required are 
    # listed below along with what is expected to be in each folder.
    #
    # hobo_csvs: contains csv files of HOBO data to process with no 
    #   header and three columns for date, time, temperature. Files should 
    #   be named as site.nestbox.mm.dd.yy.csv
    #
    # reference_files: includes one file used for matching to HOBOS
    #       i: 'exclude_sections.csv' details sections that should not be scored
    #             * initially created by the code, then fill in sections to exclude manually
    #
    # hobo_plots: all output figures, with subfolder 'initial_diagnostic'
    #
    # scripts: scripts used for processing are stored here
    #
    # processed_data_output: other output files are written here
    #

# Load packages ----
    # Load packages used in code. These may need to be installed the first time
    
        pacman::p_load(plyr, rethinking, modeest, here, lubridate, scales, svMisc, suncalc)

# Optional settings ----
    # Modify settings that change how the code runs.
    
    # Set the tolerance threshold. The TOTAL temperature change in degrees has to be larger than this 
    # amount in order to count as a bout transition. This does not require that a particular change
    # occur between two consecutive readings. This has all been tested in Celcius, but in theory
    # it should work regardless of whether readings are in F or C, except that you'll want a different
    # number. Larger numbers will get rid of more false transitions but will also miss some real ones.
        tol <- 1.5
    
    # What is the span to use for wiggliness of loess fit. This number represents the number
    # of consecutive temperature readings used for local regression. Making the number smaller
    # will result in a more wiggly line. About 75-125 seems to work pretty well for records that are taken
    # 10 seconds apart. I think that this number should be reduced to 1/3 for records taken 30 seconds 
    # apart, but I have not tested that.
        SPAN <- 80
    
    # There is a secondary filter built in that attempts to clean up nighttime bouts. The problem is that at night
    # there is often a very long on bout, but the temperature oscillates up and down a bit. The main algorithm tends
    # to score a few on and off bouts in there and can even end up counting the entire thing as an off bout if there
    # is a period at the beginning that goes down a little bit. This attempts to clean that up by counting any
    # reading taken at night that is above a set temperature as an 'on' bout regardless of the slope. There are two
    # settings that control the behavior.
    
    # Padding for sunrise and sunset. The time of sunrise and sunset is used to determine what counts as a night or 
    # day bout, but it is pretty clear that normal on off behavior typically doesn't start until a while after 
    # official sunrise and continues for a while after official sunset (at least for TRES in Ithaca NY). This
    # allows you to set a padding buffer separately for sunrise and sunset that essentially pushes the times back
    # by a set amount before applying the nighttime rules. Times should be specificed in seconds. Set to 0 to turn off.
        rise_pad <- 3e3
        set_pad <- 6e3
    
    # Readings above this value are set as 'on' even if temperature is fluctuating a bit if they occur during the night.
    # This helps to clean up some bad scoring when temperature is fairly stable all night long. To turn off change
    # to a very high number like 100.
        temp_lev <- 31

# Settings for sunrise and sunset ----        
      # This sets up the twilight time data to be used to determine 
      # day and night time bouts later on in the code.        
      
      # Set the earliest and latest day of any file (ok to provide a buffer even earlier and later)
          early.day <- as.Date("2013-01-01") # set earlier than earliest possible date
          late.day <- as.Date("2024-08-01") # set later than latest possible date
          
      # Set location for sunrise and sunset times 
      # currently using a single location. could be adjusted to allow different locations for different files
          lat.sun <- 42.503348
          lon.sun <- -76.453737
          
      # use suncalc package to get sunrise and sunset times for the year
      # be sure to adjust time zone if needed
          sun <- getSunlightTimes(date = seq.Date(early.day, late.day, by = 1), 
                              lat = lat.sun, lon = lon.sun, keep = c("sunrise", "sunset"), tz = "America/New_York")
      
      # change format around
          sun$JDate <- yday(sun$date)
          sun$rise <- hour(sun$sunrise) * 60 * 60 + minute(sun$sunrise) * 60
          sun$set <- hour(sun$sunset) * 60 * 60 + minute(sun$sunset) * 60
          sun$rise2 <- sun$rise + sun$JDate * 24 * 60 * 60
          sun$set2 <- sun$set + sun$JDate * 24 * 60 * 60
          sun$yr_doy <- paste(year(sun$date), sun$JDate, sep = "_")
          sun <- sun[, c("date", "yr_doy", "JDate", "rise", "set", "rise2", "set2")]
          sun <- subset(sun, sun$JDate > 100 & sun$JDate < 200)  # restrict to dates of breeding to speed up

# Data frame setup ----          
      # Make an object that has the names of all files to be processed
      
          mycsv <- dir(path = here::here("0_input_data"), pattern = ".csv")
      
      # Make the data frame that will hold output from each list
      
          all_bouts <- as.data.frame(matrix(nrow = 0, ncol = 23))
          colnames(all_bouts) <- c("type", "delta", "start", "end", "mu_temp", "sd_temp", "JDate", "min_temp", "max_temp",
                                   "time_start", "time_end", "duration", "unit", "nest", "file", "date", "begin.time", 
                                   "end.time", "rise", "set", "rise2", "set2", "is_night")
          
# Loop to process each raw file ----         
      # Beginning of a big for loop that will process each file in the folder  
      
          for(p in 1:length(mycsv)){	
                start2 <- Sys.time()
            
            #Read in the first file
                data <- read.csv(paste(here::here("0_input_data/"), mycsv[p], sep = ""), header = FALSE)
            
            # Change column names
                colnames(data) <- c("date", "time", "temp")
            
            # Add julian date
                data$JDate <- yday(as.POSIXct(data$date, format = "%m/%d/%y"))
                data$yr_doy <- paste(year(as.POSIXct(data$date, format = "%m/%d/%y")), data$JDate, sep = "_")
                
            # Split the time string into hours minutes and seconds
                data$time<-as.character(data$time)
                data$hour <- substr(data$time, 1, nchar(data$time) - 6)
                data$minute <- substr(data$time, nchar(data$time) - 4, nchar(data$time) - 3)
                data$sec <- substr(data$time, nchar(data$time) - 1, nchar(data$time))
                
            # Make a new time that is seconds after midnight
                data$time2 <- as.numeric(data$hour) * 60 * 60 +
                  as.numeric(data$minute) * 60 + as.numeric(data$sec)
                
            # Make a new time that is seconds after January 1st. This is unwieldy but makes all 
            # the plotting and calculations much easier since it combines date and time into
            # a single continuous measure for each year.
                data$time4 <- data$time2 + data$JDate * 24 * 60 * 60
            
            # Fit a loess smoothed regression. Adjustments to the degree of wiggliness are done above
            # by change the 'SPAN' value. This could be done on the raw values, but the loess
            # smooths out places where temp jumps up or down only briefly.
                loessM <- loess(temp ~ time4, data = data, span = SPAN/nrow(data))
                smoothedM <- predict(loessM)
                
            # This looks at the difference in temperature from time n to time n + 1 and identified
            # locations where the sign flips from positive to negative or vice versa. These are
            # inferred to be potential bout transition points.
            # first for the lower points
                upper <- which(diff(sign(diff(smoothedM))) == -2)
            # then for the upper points
                lower <- which(diff(sign(diff(smoothedM))) == 2)
            # Cut to the same length
                if(length(upper) > length(lower)){upper <- upper[1:length(lower)]}
                if(length(lower) > length(upper)){lower <- lower[1:length(upper)]}
                
            # This puts the upper and lower points together into a data frame. The first part represents
            # off bouts (upper to lower) and the second part is on bouts (lower to upper). The offset
            # in brackets is to make them line up right. Column names are changed to 'start' and 'end'
            # and another column is added to identify on vs. off bouts.
            
                ifelse(upper[1] < lower[1],
                       bouts <- as.data.frame(cbind(c(upper, lower), 
                                                    c(lower, upper[2:length(upper)], NA))),
                       bouts <- as.data.frame(cbind(c(lower, upper), 
                                                    c(upper, lower[2:length(lower)], NA))))
                colnames(bouts) <- c("start", "end")
                
            # Calculate the total change in temperature from the start to the end of each bout.
                for(i in 1:nrow(bouts)){
                  bouts$delta[i] <- smoothedM[bouts[i, "end"]] - smoothedM[bouts[i, "start"]]
                }
                mu.delt.1 <- mean(na.omit(bouts$delta[1:length(upper)]))
                mu.delt.2 <- mean(na.omit(bouts$delta[(length(upper) + 1):(length(upper) * 2)]))
                
                
                ifelse(mu.delt.1 > mu.delt.2,
                       bouts$type <- c(rep("on", length(upper)), rep("off", length(lower))),
                       bouts$type <- c(rep("off", length(lower)), rep("on", length(upper)))
                )
            
            # Add an index value that is used to sort the on and off bouts together into sequential 
            # order. This is necessary because I add on and off bouts as separate blocks to the 
            # data frame above and I want them in time sequence.
                bouts$index <- c(seq(1:length(upper)), seq(1:length(lower)))
            
            # Remove any bouts that do not have delta temperature values greater than the tolerance limit 
            # set at the top of the code. On bouts should be > tol and off bouts < -tol.
                bouts2 <- subset(bouts, bouts$type == "on" & bouts$delta > tol |
                               bouts$type == "off" & bouts$delta < -tol)
            
            # Reorder the bouts into temporal sequnce based on the start time.
                bouts4 <- bouts2[order(bouts2$start),]
            
            # This is a bit tricky. What the loop here does is to extend the identified bouts that have had short bouts
            # removed from the middle of them based on the temperature threshold above. For example, if two long off
            # bouts were separated by a short on bout that only went up 1 degres, that would result in the on bout
            # being deleted, but still two separate off bouts. This extends the first of those off bouts all the way
            # through the deleted on bout AND the subsequent off bout (and this could go on for an indefinite number
            # of bouts). This has the effect of cleaning up areas where temperature is fairly stable for a long time
            # and in previous iterations would have been scored as lots of short on/off bouts.
                for(i in 1:(nrow(bouts4)-2)){
                  if(bouts4$end[i] == bouts4$start[i + 1]){
                    bouts4$end.new[i] <- bouts4$start[i + 1]
                  }
                  if(bouts4$end[i] != bouts4$start[i + 1]){
                    if(sum(bouts4$type[i:nrow(bouts4)] != bouts4$type[i]) > 0){
                      bouts4$end.new[i] <- bouts4$start[min(which(bouts4$type[i:nrow(bouts4)] != bouts4$type[i])) + i - 1]
                    }
                  }
            }
            ## NOTE: if a file ends with a string of on or off bouts in a row this will give an error about
            # 'no non-missing' arguments, but this can be ignored. 
            
            # Now that the first bout is extended when necessary, this deletes duplicate bouts, which are identified
            # based on shared end times. Note, these are really nested rather than fully duplicated bouts, but they
            # need to bre cleaned out. It is important that the bouts are sorted in temporal order here because only
            # the first one is saved.
                bouts5 <- subset(bouts4, !duplicated(bouts4$end.new))
                bouts5 <- subset(bouts5, is.na(bouts5$end.new) == FALSE)
                bouts5 <- subset(bouts5, bouts5$end.new > bouts5$start)
                
            # Calculate a new delta temperature value now that the bouts are extended
                for(i in 1:nrow(bouts5)){
                  bouts5$delta2[i] <- data$temp[bouts5[i, "end.new"]] - data$temp[bouts5[i, "start"]]
                }
            
            ## Add on off to data
                data$on_off <- rep(NA)
                for(i in 1:nrow(bouts5)){
                  data[bouts5[i, "start"] : bouts5[i, "end.new"], "on_off"] <- bouts5$type[i]
                }
            
            ## Add in sunrise and sunset times
                data <- join(data, sun, "yr_doy")
            ## Add in the padding buffer for sunrise and sunset
                data$rise3 <- data$rise2 + rise_pad
                data$set3 <- data$set2 + set_pad
            ## Add column for if the record is at night
                data$after_rise <- data$rise3 - data$time4  # positive if before sunrise
                data$after_set <- data$set3 - data$time4    # positive if before sunset
                data$sun_sign <- sign(data$after_rise) + sign(data$after_set)
                joiner <- data.frame(sun_sign = c(-2, 0, 2),
                                     is_night = c(1, 0, 1))
                data <- join(data, joiner, "sun_sign")
            
            ## Change on off to numbers
                data$on_off <- gsub("off", 0, data$on_off)
                data$on_off <- gsub("on", 1, data$on_off)
                data$on_off <- as.numeric(data$on_off)
                
            ## change nighttime bouts over certain temperature to on bouts (see above)
                data.s <- data
                data <- subset(data, is.na(data$is_night) == FALSE)
                data[data$temp > temp_lev & data$is_night == 1, 10] <- 1
                
            # Make new list of on off break points based on cleaned up version
            # first for the lower points
                upper2 <- which(diff(data$on_off) == 1) + 1
            # then for the upper points
                lower2 <- which(diff(data$on_off) == -1) + 1
            # these need to be the same length, check and trim one off
                if(length(upper2) > length(lower2)){lower2 <- lower2[1:length(upper2)]}
                if(length(lower2) > length(upper2)){upper2 <- upper2[1:length(lower2)]}
                
            # This puts the upper and lower points together into a data frame. The first part represents
            # off bouts (upper to lower) and the second part is on bouts (lower to upper). The offset
            # in brackets is to make them line up right. Column names are changed to 'start' and 'end'
            # and another column is added to identify on vs. off bouts.
                ifelse(upper2[1] < lower2[1],
                       xbouts <- as.data.frame(cbind(c(upper2, lower2), 
                                                     c(lower2, upper2[2:length(upper2)], NA))),
                       xbouts <- as.data.frame(cbind(c(lower2, upper2), 
                                                     c(upper2, lower2[2:length(lower2)], NA))))
                colnames(xbouts) <- c("start", "end")
            
            # Calculate the total change in temperature from the start to the end of each bout.
                for(i in 1:nrow(xbouts)){
                  xbouts$delta[i] <- smoothedM[xbouts[i, "end"]] - smoothedM[xbouts[i, "start"]]
                }
                mu.delt.1 <- mean(na.omit(xbouts$delta[1:length(upper2)]))
                mu.delt.2 <- mean(na.omit(xbouts$delta[(length(upper2) + 1):(length(upper2) * 2)]))
                
                
                ifelse(mu.delt.1 > mu.delt.2,
                       xbouts$type <- c(rep("on", length(upper2)), rep("off", length(lower2))),
                       xbouts$type <- c(rep("off", length(lower2)), rep("on", length(upper2)))
                )
            
            # Add an index value that is used to sort the on and off bouts together into sequential 
            # order. This is necessary because I add on and off bouts as separate blocks to the 
            # data frame above and I want them in time sequence.
                xbouts$index <- c(seq(1:length(upper2)), seq(1:length(lower2)))
            
            # Determine the number of readings in each bout
                xbouts$count <- xbouts$end - xbouts$start
            
            # Reorder the bouts into temporal sequnce based on the start time.
                xbouts2 <- xbouts[order(xbouts$start),]
                xbouts2 <- subset(xbouts2, is.na(xbouts2$end) == FALSE & is.na(xbouts2$start) == FALSE)
                
            # Make a diagnostic plot of the identified on and off bout transition points. This can be used both to verify
            # that the code and settings used are working well and also as a way to identify sections of the records
            # that should be excluded and not scored (see filter_HOBO_code.R).
            
            # Set up to print to a pdf. This is VERY wide for large files but needs to be in order to be visible. Zoom in!
                pdf(here::here("2_plots/initial_diagnostic", paste0(mycsv[p], "trace", ".pdf")), width = nrow(data) / 2000 * 6, height = 5)
            # Set up the plot
                plot(data$time4, data$temp, type = "n", main = mycsv[p], xlab = "Time", ylab = "Temp (C)", xaxt = "n", ylim = c(14, 40))
            # Add an axis with reasonably spaced labels
                axis(1, seq(0, 1e9, 10000))
                axis(3, seq(0, 120 * 24 * 60 * 60 + 3600 * 2000, 3600), labels = rep("", 4881))  # adds tick marks for each hour
            # Add reference lines for temperature since you'll have to scroll way away from y axis.
                temp_lines <- seq(0, 50, 5)
                temp_cols <- c(rep("gray50", 6), "coral3", rep("gray50", 4))
                for(i in 1:length(temp_lines)){
                  abline(h = temp_lines[i], lty = 2, col = temp_cols[i])
                }
            # Add trace line of the raw tempearture data
                lines(data.s$time4, data.s$temp)
            # Add trace line of the loess smoothed line
                lines(smoothedM, x = data.s$time4[1:length(smoothedM)], col = "green")
            # Subset for plotting purposes to avoid duplicating points on start/end times
                xbouts3 <- subset(xbouts2, xbouts2$type == "on")
            # Add points for identified start and end times
                points(data$time4[xbouts3$start], 
                       data$temp[xbouts3$start], col = "blue", pch = 16)
                points(data$time4[xbouts3$end], 
                       data$temp[xbouts3$end], col = "red", pch = 16)
            ## Add all poitns of on and off
                points(data$time4, rep(15, nrow(data)), col = as.factor(data$on_off), pch = 16, cex = 0.8)
            # Add vertical lines for sunrise and sunset
                for(i in 1:nrow(sun)){
                  abline(v = sun$rise2[i], col = "goldenrod", lwd = 2)
                  abline(v = sun$set2[i], col = "gray40", lwd = 2)
                }
            # Turn off the pdf
                dev.off()
            
            ## Create a list that separates out each identified bout into one data frame and add a type
            # column to identify on vs. off bouts.
                bout_details <- list(NA)
                for(k in 1:nrow(xbouts2)){
                  bout_details[[k]] <- data[xbouts2[k, "start"] : xbouts2[k, "end"], ]
                  bout_details[[k]]$type <- rep(xbouts2$type[k], nrow(bout_details[[k]]))
                }  
            
            # This is identifying which elements are on or off bouts to make a plot
                store <- rep(NA, length(bout_details))
                for(i in 1:length(bout_details)){
                  store[i] <- bout_details[[i]][1, "type"]
                }
            
            # Now make subsets that include only on and off bouts, respectively. Only used for plotting.
                off_sub <- bout_details[c(which(store == "off"))]
                on_sub <- bout_details[c(which(store == "on"))]
                
            # Make a plot that shows lines for all on and off bouts
            # Set up the plot
                pdf(here::here("2_plots/initial_diagnostic", paste0(mycsv[p], "bouts", ".pdf")), width = 10, height = 7)
                par(mfrow = c(1, 2))
            # first on bouts
                plot(1, 1, type = "n", ylim = c(15, 40), xlim = c(0, 2100), bty = "n", xaxt = "n", yaxt = "n",
                     ylab = "Temperature (C)", xlab = "Minutes", main = "On Bouts", xaxs = "i")
                axis(1, seq(-600, 7200*10, 600), labels = seq(-10, 120*10, 10))
                axis(2, seq(0, 40, 5), las = 2)
                for(i in 1:length(on_sub)){
                  if(nrow(on_sub[[i]]) < 220){
                    lines(seq(1, nrow(on_sub[[i]]), 1) * 10, on_sub[[i]]$temp, col = col.alpha("black", 0.1))
                  }
                }
            # now off bouts
                plot(1, 1, type = "n", ylim = c(15, 40), xlim = c(0, 2100), bty = "n", xaxt = "n", yaxt = "n",
                     ylab = "Temperature (C)", xlab = "Minutes", main = "Off bouts", xaxs = "i")
                axis(1, seq(-600, 7200*10, 600), labels = seq(-10, 120*10, 10))
                axis(2, seq(0, 40, 5), las = 2)
                for(q in 1:length(off_sub)){
                  if(nrow(off_sub[[q]]) < 220){
                    lines(seq(1, nrow(off_sub[[q]]), 1) * 10, off_sub[[q]]$temp, col = col.alpha("black", 0.1))
                  }
                }
        # Turn off the pdf
            dev.off()
  
# Summarize ----  
        # Summarize the bout details object with one row per bout
              bout_out <- xbouts2[, c("type", "delta", "start", "end")]
              for(i in 1:nrow(bout_out)){
                sub <- bout_details[[i]]
                bout_out$mu_temp[i] <- mean(sub$temp)
                bout_out$sd_temp[i] <- sd(sub$temp)
                bout_out$JDate[i] <- sub$JDate[1]
                bout_out$min_temp[i] <- min(sub$temp)
                bout_out$max_temp[i] <- max(sub$temp)
                bout_out$time_start[i] <- sub$time2[1]
                bout_out$time_end[i] <- sub$time2[nrow(sub)]
                bout_out$duration[i] <- sub$time4[nrow(sub)] - sub$time4[1]
              } 
              
              unit <- strsplit(mycsv[p], "[.]")[[1]][1]
              nest <- strsplit(mycsv[p], "[.]")[[1]][2]
              
              bout_out$unit <- unit
              bout_out$nest <- nest
              bout_out$file <- mycsv[p]
        
        ## Join the summary data to sunrise and set to determine if it is at night
              bout_out <- join(bout_out, sun, "yr_doy")
              bout_out$is_night <- NA
              for(i in 1:nrow(bout_out)){
                if(bout_out$time_start[i] < bout_out$rise[i]){bout_out$is_night[i] <- 1}
                if(bout_out$time_start[i] > bout_out$set[i]){bout_out$is_night[i] <- 1}
                if(bout_out$time_start[i] > bout_out$rise[i] &
                   bout_out$time_start[i] < bout_out$set[i]){bout_out$is_night[i] <- 0}
              }
              
              all_bouts <- rbind(all_bouts, bout_out)
              
              
              print(paste(p, "out of", length(mycsv), sep = " "))
              
            }
# Save ----
    ## Save the final all bouts object. If no deletions are needed then the file can be used as is
    ## with the summarize_HOBO_bouts.R code. If there are sections of the code that you want to 
    ## exclude because the record is not good then you need to first construct an exclusion
    ## file and run the exclud_HOBO_sections.R code before summarizing.
        all_bouts$exclude_all <- 0
        all_bouts$exclude_temp <- 0
        saveRDS(all_bouts, file = here::here("processed_data_output/all_bouts.rds"))