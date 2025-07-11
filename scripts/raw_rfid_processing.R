# Script information ----
## This code is written to process RFID files used on tree swallow nest boxes.
## Last updated 9/16/2024
## Written by Conor Taff: cct663@gmail.com
##

# Load packages & settings ----
  # Load packages used in code. These may need to be installed the first time
    pacman::p_load(plyr, lubridate, reshape2, here)

# Read reference files ----

  ## Dates and info file. This has RFID info paired with nest dates (hatching,
  ## initiation, fledge/failure). The unit-box column must match the unitbox
  ## name used in the saved file name.
      RFIDRef <- read.delim(here("reference_files", "RFIDRef.txt"))
      
    # Do some wrangling to clean this up for nests with missing dates
      RFIDRef$j_hatch <- gsub("999", "NA", RFIDRef$j_hatch)
      RFIDRef$j_hatch <- as.numeric(RFIDRef$j_hatch)
      for(i in 1:nrow(RFIDRef)){
        if(is.na(RFIDRef$j_fate[i]) == TRUE){RFIDRef$j_fate[i] <- 220}
        if(RFIDRef$j_fate[i] == 999){RFIDRef$j_fate[i] <- 220}
        if(RFIDRef$j_first_egg[i] == 999){RFIDRef$j_first_egg[i] <- RFIDRef$j_incub[i] - 5}
      }
      RFIDRef$uby <- paste(RFIDRef$unitbox, RFIDRef$year, sep = "_")

  ## Full RFID reference. This is the full RFID reference library. It should list
  ## every possible RFID code that has been deployed to look for matches. 
      RFIDFULL <- read.delim(here("reference_files", "FULLRFIDReference.txt"))
      
  ## A list of all PIT Tags that were swapped on individuals. Used to switch all to original
      swap_PIT <- read.delim(here("reference_files", "swapped_pits.txt"))
      
  ## Swap rfid reference PIT tags to the original PIT for each bird
      for(i in 1:nrow(swap_PIT)){
        RFIDRef$f_rfid <- gsub(swap_PIT$new_rfid[i], swap_PIT$original_rfid[i], RFIDRef$f_rfid)
        RFIDRef$m_rfid <- gsub(swap_PIT$new_rfid[i], swap_PIT$original_rfid[i], RFIDRef$m_rfid)
      }

# Build breeding stage ----
  ## Construct a breeding stage reference file from the RFIDRef file. This just
  ## uses the dates already in the RFIDRef file to make the large breeding stage
  ## matrix that says what stage each nest was at on each day of the season.

      br_stg <- as.data.frame(matrix(nrow = 250, ncol = nrow(RFIDRef)))

      ## Now filling in full matrix of breeding stages
      for(i in 1:nrow(RFIDRef)){  
        
          vec <- rep(NA, 250) ## filling for doy 0 to 250
        ## For all nests
        # Fill in post fate known date
          vec[RFIDRef$j_fate[i]:250] <- seq(30, (30 + 250 - RFIDRef$j_fate[i]), 1)
        # Fill in pre-laying
          vec[1:(RFIDRef$j_first_egg[i] - 1)] <- 
            seq(-RFIDRef$j_first_egg[i] - 28, -30, 1)
        # Fill in for laying
          vec[RFIDRef$j_first_egg[i]:(RFIDRef$j_incub[i] - 1)] <-
            seq(RFIDRef$j_first_egg[i] - RFIDRef$j_incub[i] - 19, -20, 1)
        
        ## For nests that successfully hatched
        if(is.na(RFIDRef$j_hatch[i]) == FALSE){
          # Fill in from hatch to fate known
            vec[RFIDRef$j_hatch[i]:(RFIDRef$j_fate[i] - 1)] <- 
              seq(0, RFIDRef$j_fate[i] - RFIDRef$j_hatch[i] - 1, 1)
          # Fill in incubation
            vec[RFIDRef$j_incub[i]:(RFIDRef$j_hatch[i] - 1)] <-
              seq(RFIDRef$j_incub[i] - RFIDRef$j_hatch[i], -1, 1)
          
        }
        
        ## For nests that failed before hatching
          if(is.na(RFIDRef$j_hatch[i]) == TRUE){
          # Fill in from incubation to failure date
            vec[RFIDRef$j_incub[i]:(RFIDRef$j_fate[i] - 1)] <-
              seq(RFIDRef$j_incub[i] - RFIDRef$j_fate[i], -1, 1)
        }
        
        
        br_stg[, i] <- vec
      }

    br_stg2 <- as.data.frame(t(br_stg))
    colnames(br_stg2) <- seq(1, 250)
    br_stg2$uby <- RFIDRef$uby
    RFIDRef <- join(RFIDRef, br_stg2, "uby")

# Build other objects ---- 
    ## This section sets up a few objects that will be used later on.
    
    ## This reads in the list of file names (should be "unit.box.mm.dd.yyyy.txt" files).
    ## Naming convention is important because unit and box info will be stripped from 
    ## the file name. The date part actually doesn't matter because each read has a timestamp.
        DoList <- list.files(here("rfid_csvs")) ## Folder with only RFID files to process

# Cleaning Loop ----
      
      # Below here is where the file processing starts in earnest. Everything above was just setting up
      # objects and reference files to be used later on.
      
    ## This is the first big for loop that starts processing each file.
      # This one does a few different things, including:
      # 1. Read raw RFID text and use a regular expression to strip all lines before the hits.
      # 2. Add columns to the reads by stripping unit and box from file name and joining to 
        # reference tables that have additional columns.
      # 3. For unknown RFID numbers, calculate how similar they are to the known male and female
        # at this box (# of digits different) and if they are 1 or 2 digits off swap the RFID
        # number out for the focal bird (this corrects for mis-reads that are common).
      
        # simple object to fix some years that are recorded with 2 vs. 4 digits
          year_fix <- data.frame(yearX = c(seq(13, 24, 1), seq(2013, 2024, 1)),
                               year = rep(seq(2013, 2024, 1), 2))
      
      ### Note that R sometimes gives some warnings at the end of this loop about invalid factor levels or 
      ### incomplete final lines, etc. For the most part, these are fine (they might result in one lost line
      ### of RFID output). You might be able to resolve them through formatting of the raw text files, but 
      ### if the output is being produced correctly (check the modified files) they are probably safely ignored.
      
          for(j in 1:length(DoList)){
            
              Sys.setlocale('LC_ALL', 'C')
            
            # Set the first file to process
            
              file <- DoList[j]
            
            # Read in all the lines from the first file and then use a regular expression to find
            # the start of the actual reads by matching to '????Data'. Store the start position
            # returned and then read in the file again as an R object with just the relevant lines.
            
              Lines <- readLines(here("rfid_csvs", file), warn = FALSE)
              markers <- grepl("????DATA", Lines)
              starter <- rle(markers)$lengths[1] + 2
              c <- read.delim(here("rfid_csvs", file), sep = " ", skip = starter, header = FALSE)
              colnames(c) <- c("rfid", "Date", "Time")
              
            # swap any rfids for birds that have had their tags replaced
              for(i in 1:nrow(swap_PIT)){
                c$rfid <- gsub(swap_PIT$new_rfid[i], swap_PIT$original_rfid[i], c$rfid)
              }
            
            # Rename columns and join to reference sheets. Use the file name to create
            # columns that can be used to join to other information.		
            
              c$rfid <- as.character(c$rfid)
              RFIDFULL$rfid <- as.character(RFIDFULL$rfid)
              c <- join(c, RFIDFULL[, c("band", "rfid", "sex")], "rfid")
              NestDetail <- strsplit(DoList[j], "[.]")
              Box <- as.factor(paste(unlist(NestDetail)[1], unlist(NestDetail)[2], sep = "_"))
              c$unitbox <- rep(Box, nrow(c))
              c$unitbox <- as.character(c$unitbox)
              c$yearX <- NestDetail[[1]][5]
              c <- plyr::join(c, year_fix, "yearX")
              c$uby <- paste(c$unitbox, c$year, sep = "_")
              
              RFIDRef$uby <- as.character(RFIDRef$uby)
              c <- join(c, RFIDRef, "uby")
            
            # This pairs the RFID read to first the male and then the female
            # that own the box and compares each character one-by-one assigning a 1 if they 
            # match and a 0 if they don't. The 10 comparisons are summed so a 10 is a perfect
            # match. This value is saved in a new column. This is necessary because the readers
            # sometimes 'slip' and read a single digit wrong. So with 
            # a lot of reads for a focal bird you often get smaller peaks that are just one or
            # two characters different from the focal birds.
            
                c$MatchFemale <- ifelse(substr(c$rfid, 1, 1) == substr(c$f_rfid, 1, 1), 1, 0) +
                  ifelse(substr(c$rfid, 2, 2) == substr(c$f_rfid, 2, 2), 1, 0) +
                  ifelse(substr(c$rfid, 3, 3) == substr(c$f_rfid, 3, 3), 1, 0) +
                  ifelse(substr(c$rfid, 4, 4) == substr(c$f_rfid, 4, 4), 1, 0) +
                  ifelse(substr(c$rfid, 5, 5) == substr(c$f_rfid, 5, 5), 1, 0) +
                  ifelse(substr(c$rfid, 6, 6) == substr(c$f_rfid, 6, 6), 1, 0) +
                  ifelse(substr(c$rfid, 7, 7) == substr(c$f_rfid, 7, 7), 1, 0) +
                  ifelse(substr(c$rfid, 8, 8) == substr(c$f_rfid, 8, 8), 1, 0) +
                  ifelse(substr(c$rfid, 9, 9) == substr(c$f_rfid, 9, 9), 1, 0) +
                  ifelse(substr(c$rfid, 10, 10) == substr(c$f_rfid, 10, 10), 1, 0)
                c$MatchMale<-ifelse(substr(c$rfid, 1, 1) == substr(c$m_rfid, 1, 1), 1, 0) +
                  ifelse(substr(c$rfid, 2, 2) == substr(c$m_rfid, 2, 2), 1, 0) +
                  ifelse(substr(c$rfid, 3, 3) == substr(c$m_rfid, 3, 3), 1, 0) +
                  ifelse(substr(c$rfid, 4, 4) == substr(c$m_rfid, 4, 4), 1, 0) +
                  ifelse(substr(c$rfid, 5, 5) == substr(c$m_rfid, 5, 5), 1, 0) +
                  ifelse(substr(c$rfid, 6, 6) == substr(c$m_rfid, 6, 6), 1, 0) +
                  ifelse(substr(c$rfid, 7, 7) == substr(c$m_rfid, 7, 7), 1, 0) +
                  ifelse(substr(c$rfid, 8, 8) == substr(c$m_rfid, 8, 8), 1, 0) +
                  ifelse(substr(c$rfid, 9, 9) == substr(c$m_rfid, 9, 9), 1, 0) +
                  ifelse(substr(c$rfid, 10, 10) == substr(c$m_rfid, 10, 10), 1, 0)
            
            # This for loop goes through and looks at each RFID hit that did not match any known
            # birds in the full reference database. Presumably, these are mostly mis-reads.
            # It then asks whether these numbers are only 1 or 2 characters different from
            # the focal birds at the nest. If the answer is yes then it replaces the read
            # with that of the focal bird that matches.	
            
                for(r in 1:nrow(c)){	
                  ifelse(is.na(c$Band[r]) == TRUE,
                         ifelse(c$MatchFemale[r] > 7,
                                c$rfid[r] <- c$f_rfid[r],
                                ifelse(c$MatchMale[r] > 7,
                                       c$rfid[r] <- c$m_rfid[r],
                                       c$rfid[r] <- c$rfid[r])),
                         c$rfid[r] <- c$rfid[r])	
                }	
            
            # This strips out most of the (now extraneous) columns that we have made and then
            # joins anew to the full RFID reference sheet (results will differ now because
            # mis-reads have been subbed out).
            
              c <- c[, c("rfid", "Date", "Time", "year", "unitbox", "uby")]
              c <- join(c, RFIDFULL[, c("band", "rfid", "sex")], "rfid")
            
            # This gets rid of any more RFID reads that are still not matching. There should not
            # be many of these, but some could be mis-reads of other bands (rare) or reads of
            # the tag that we use for testing when deploying (more likely) or reads of birds
            # that never made it into the full RFID reference spreadsheet somehow (possible).
            
              c <- subset(c, is.na(c$band) == FALSE)	
            
            # Now that is all done, strip down to the essential info and join to the JDates
            # object to convert dates to meaningful values.
            
              c <- c[,c("rfid", "Date", "Time", "unitbox", "uby", "year")]
              c$yday <- yday(mdy(c$Date))
              c$year_rd <- year(mdy(c$Date))
              c <- subset(c, is.na(c$yday) == FALSE)
              c <- subset(c, c$year == c$year_rd) # gets rid of any old reads from previous years
            
            # now wrangle to the right time and date info  
              if(nrow(c) > 0){
                
                ##  Convert the time string to continuous seconds
                  c$Seconds <- period_to_seconds(hms(c$Time))
                
                ## Join the data object to the RFID reference sheet that has by nest info based on "UnitBox"
                  c <- join(c, RFIDRef[, c("uby", "f_band", "f_rfid", "m_band", "m_rfid", "j_first_egg",
                                             "j_incub", "j_hatch", "j_fate", "j_install")], 
                             "uby", type = "left", match = "first")
                
                ## Make a new column that indicates whether a read is from the focal male/female or a different bird.
                  c$focal <- ifelse(as.character(c$f_rfid) == as.character(c$rfid), "Yes",
                                   ifelse(as.character(c$m_rfid) == as.character(c$rfid), "Yes", "No"))
                  
                ## Calculate the offset between observation and hatch date (e.g., 2 days before hatch
                ## would yield a result of -2)
                    c$Offset <- c$yday - c$j_hatch
                  
                ## Take a subset of the file that is from provisioning (offset > -1)
                    c <- subset(c, c$Offset > -1)
                    c <- subset(c, c$focal == "Yes") # subst to only focal birds
                    
                ## Merge the version of d2 created so far into one giant ojbect that includes 
                # all observations recorded from the entire batch of files.
                    if(j == 1){all_prov_merged <- c}				
                    if(j > 1){all_prov_merged <- rbind(all_prov_merged, c)}
                    
                    ## End of the if loops from above
                  
                }  
            
            ## Print this at the end of each pass through the loop just to be able to see that R
            # is still working...
            
              print(paste("Cycle", j, "out of", length(DoList), sep = " "))
            
          }
      
      ## Now that all the files have been processed, the object that has all the reads
      # Get rid of duplicate observations that are introduced by reading in overlapping RFID files. First makes a unique
      # identifier column and then discards duplicates.
          all_prov_merged$UniqueRecord <- paste(all_prov_merged$rfid, all_prov_merged$yday,
                                                all_prov_merged$Seconds, sep = "_")
          all_prov_merged <- subset(all_prov_merged, !duplicated(UniqueRecord))
      
      # Make a column that records time in continuous seconds starting with the first of the year. This is kind of silly, 
      # but is easier to work with for calculating intervals, etc.
          all_prov_merged$FullTime <- all_prov_merged$yday*24*60*60 + all_prov_merged$Seconds
          write.table(all_prov_merged, here("processed_data_output", "all_prov_merged.txt"), sep = "\t")

      
## Calculate feeding rate for each hour ----
        d2 <- all_prov_merged # if this isn't loaded, read in from saved object
        d2 <- subset(d2, d2$Offset > -1 & d2$Offset < 21)  
    
    ## This applies time thresholds to feeding visits based on Vitousek et al. 2018 PRSB.
    # The thresholds were determined from video files. See supplement of that paper.
        thresholds <- as.data.frame(matrix(nrow = 18, ncol = 2))
        colnames(thresholds) <- c("Offset", "Limit")
        thresholds$Offset <- seq(from = 1, to = 18, by = 1)
        thresholds$Limit <- c(rep(136.5, 3), rep(55.5, 3), rep(36.5, 3),
                              rep(20, 3), rep(11, 3), rep(25.5, 3))
    
    ## Add the thresholds to the main data frame
        d2<-join(d2, thresholds, "Offset")
    
    ## Add the time of next reading and calculate time gap between reads
        d2$Time2 <- as.character(d2$Time)
        d2$Hour <- floor(d2$Seconds / 60 / 60)

    ## make more columns
        d2$ubyd <- paste(d2$uby, d2$yday, sep = "_")
        d2$ubydh <- paste(d2$ubyd, d2$Hour, sep = "_")
    
    ## Make separate male and female objects
        d2f <- subset(d2, as.character(d2$f_rfid) == as.character(d2$rfid))
        d2m <- subset(d2, as.character(d2$m_rfid) == as.character(d2$rfid))
    
    ## Add in the time of the next read to each row 
        d2f$FullTimeTneg1 <- c(0, d2f$FullTime[1:(nrow(d2f) - 1)])
        d2m$FullTimeTneg1 <- c(0, d2m$FullTime[1:(nrow(d2m) - 1)])
    
    ## Calculate time gap between subsequent reads   
        d2f$Gap <- d2f$FullTime - d2f$FullTimeTneg1
        d2m$Gap <- d2m$FullTime - d2m$FullTimeTneg1
    
    ## Remove reads that are too close together based on thresholds  
        d3f <- subset(d2f, d2f$Gap > d2f$Limit)
        d3m <- subset(d2m, d2m$Gap > d2m$Limit)
    
    ## Calculate feeding rate for each hour    
        f_feed <- d3f %>%
          group_by(unitbox, uby, ubyd, ubydh, year, yday, f_band, f_rfid, j_first_egg, j_incub, j_hatch, j_fate, j_install, Offset, Hour) %>%
          summarize(f_feeds = n())
        
        m_feed <- d3m %>%
          group_by(ubydh, m_band, m_rfid) %>%
          summarize(m_feeds = n())
    
    # join male feeding rate to female feeding rate
        feeds <- plyr::join(f_feed, m_feed, "ubydh", "left", "first")
    
    # add in the total number of male and female reads without any filtering by threshold    
        d4f <- d2f %>%
          group_by(ubydh) %>%
          summarize(f_reads = n())
        
        d4m <- d2m %>% 
          group_by(ubydh) %>%
          summarize(m_reads = n())
        
        feeds <- plyr::join(feeds, d4f, "ubydh", "left", "first")
        feeds <- plyr::join(feeds, d4m, "ubydh", "left", "first")
        
    # save the feeds object which now has hourly feeding rate for every nest
        write.table(feeds, here("processed_data_output", "feed_per_hour.txt"), sep = "\t")
        feeds <- read.delim(here::here("processed_data_output", "feed_per_hour.txt"))
    
        
        
        
        
        
# some exploratory plots and data wrangling ----        
    ggplot(feeds[feeds$Offset > 0 & feeds$Offset < 17, ], mapping = aes(x = Hour, y = f_feeds)) + 
      geom_smooth(color = "coral3") + 
      facet_wrap(~as.factor(Offset)) + 
      xlim(c(6, 20)) + 
      ylab("Feeding trips per hour") + 
      geom_smooth(mapping = aes(x = Hour, y = m_feeds), color = "slateblue")
    
    ggplot(feeds[feeds$Hour > 5 & feeds$Hour < 20, ], mapping = aes(x = Offset, y = f_feeds)) +
      geom_smooth(color = "coral3") +
      #facet_wrap(~cut(Hour, breaks = 4)) +
      xlim(c(0, 20)) +
      ylab("Feeding trips per hour") +
      geom_smooth(mapping = aes(x = Offset, y = m_feeds), color = "slateblue") +
      xlab("Nestling age")
    
    feeds$yr_day_hr <- paste(feeds$year, feeds$yday, feeds$Hour, sep = "_")
    feeds <- plyr::join(feeds, dw2, "yr_day_hr", "left", "first")
    
    feeds$unit <- sub("_.*", "", feeds$unitbox)
    
    f_plot <- feeds[feeds$Hour > 9 & feeds$Hour < 18 & feeds$Offset > 5 & feeds$Offset < 13 & is.na(feeds$ambient_C) == FALSE, ]
    f_plot$amb_rd <- round(f_plot$ambient_C, 0)
    f_plot2 <- f_plot %>%
      group_by(amb_rd) %>%
      summarise(f_feed_av = mean(na.omit(f_feeds)), m_feed_av = mean(na.omit(m_feeds)), f_fd_se = sd(na.omit(f_feeds)) / sqrt(n()),
                m_fd_se = sd(na.omit(m_feeds)) / sqrt(n()))
    
    ggplot(f_plot, 
           mapping = aes(x = ambient_C, y = f_feeds)) +
      geom_smooth(color = "coral3", fill = "coral3", alpha = 0.3) +
      #geom_point(color = "coral3", alpha = 0.3) +
      #facet_wrap(~cut(ambient_C, breaks = 4)) +
      #xlim(c(6, 12)) +
      ylab("Feeding trips per hour") +
      geom_smooth(mapping = aes(x = ambient_C, y = m_feeds), color = "slateblue", fill = "slateblue", alpha = 0.3) +
      xlab("Temperature (C)") +
      theme_bw() +
      theme(panel.grid = element_blank(), axis.text = element_text(size = 12), axis.title = element_text(size = 14)) +
      geom_segment(data = f_plot2, mapping = aes(x = amb_rd-0.1, xend = amb_rd-0.1, y = f_feed_av + f_fd_se, yend = f_feed_av - f_fd_se), color = "coral3") +
      geom_point(data = f_plot2, mapping = aes(x = amb_rd-0.1, y = f_feed_av), shape = 21, fill = "coral3", size = 1.5) +
      geom_segment(data = f_plot2, mapping = aes(x = amb_rd+0.1, xend = amb_rd+0.1, y = m_feed_av + m_fd_se, yend = m_feed_av - m_fd_se), color = "slateblue") +
      geom_point(data = f_plot2, mapping = aes(x = amb_rd+0.1, y = m_feed_av), shape = 21, fill = "slateblue", size = 1.5) +
      scale_color_manual(name = "Sex", values = c("Females" = "coral3", "Males" = "slateblue")) +
      scale_fill_manual(name = "Sex", values = c("Females" = "coral3", "Males" = "slateblue")) +
      guides(color = guide_legend(override.aes = list(fill = c("coral3", "slateblue"), shape = 21)), 
             fill = guide_legend(override.aes = list(shape = 21))) +
      # Manually add the legend labels
      labs(color = "Sex", fill = "Sex")
    
    
    f_plot_l <- f_plot %>%
      pivot_longer(cols = c(f_feeds, m_feeds), names_to = "sex", values_to = "f_hour")
    f_plot_l$sex <- gsub("f_feeds", "Females", f_plot_l$sex)
    f_plot_l$sex <- gsub("m_feeds", "Males", f_plot_l$sex)
    ggplot(f_plot_l, 
           mapping = aes(x = ambient_C, y = f_hour, by = as.factor(year))) +
      geom_smooth(alpha = 0.3, mapping = aes(color = sex, fill = sex), se = FALSE) +
      ylab("Feeding trips per hour") +
      xlab("Temperature (C)") +
      facet_wrap(~sex) +
      theme_bw() +
      theme(panel.grid = element_blank(), axis.text = element_text(size = 12), axis.title = element_text(size = 14)) +
      #geom_segment(data = f_plot2, mapping = aes(x = amb_rd-0.1, xend = amb_rd-0.1, y = f_feed_av + f_fd_se, yend = f_feed_av - f_fd_se), color = "coral3") +
      #geom_point(data = f_plot2, mapping = aes(x = amb_rd-0.1, y = f_feed_av), shape = 21, fill = "coral3", size = 1.5) +
      #geom_segment(data = f_plot2, mapping = aes(x = amb_rd+0.1, xend = amb_rd+0.1, y = m_feed_av + m_fd_se, yend = m_feed_av - m_fd_se), color = "slateblue") +
      #geom_point(data = f_plot2, mapping = aes(x = amb_rd+0.1, y = m_feed_av), shape = 21, fill = "slateblue", size = 1.5) +
      scale_color_manual(values = c("coral3", "slateblue")) +
      scale_fill_manual(values = c("coral3", "slateblue")) +
      theme(legend.title = element_blank())
      
    
    
    
    
    ggplot(feeds[feeds$Hour >9 & feeds$Hour < 17 & is.na(feeds$ambient_C) == FALSE, ],
           mapping = aes(x = ambient_C, y = f_feeds)) +
      geom_smooth(color = "coral3", fill = "coral3", alpha = 0.3) +
      ylab("Feeding trips per hour") +
      facet_wrap(~cut(Offset, breaks = c(0, 3, 6, 9, 12, 15, 18)), scales = "free_y") +
      geom_smooth(mapping = aes(x = ambient_C, y = m_feeds), color = "slateblue", fill = "slateblue", alpha = 0.3) +
      xlab("Temperature (C)")
      
    