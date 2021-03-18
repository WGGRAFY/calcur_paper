

### script to subset the WareHouse.All.Ages.Env object


# create temporary object
WareHouse.All.Ages.Env_2 <- WareHouse.All.Ages.Env

# add in month
WareHouse.All.Ages.Env_2$mon <- as.numeric(substr(WareHouse.All.Ages.Env_2$datetime_utc_iso,6,7))


# temporary object filtered by survey (and month and years for AFSC slope survey)
tmp  <- WareHouse.All.Ages.Env_2[WareHouse.All.Ages.Env_2$project %in% 
                         c("Groundfish Slope Survey","Groundfish Triennial Shelf Survey","Groundfish Slope and Shelf Combination Survey"),]
AFSC_slope <-   WareHouse.All.Ages.Env_2[WareHouse.All.Ages.Env_2$project == "AFSC/RACE Slope Survey" & WareHouse.All.Ages.Env_2$mon<=10 & WareHouse.All.Ages.Env_2$year>=1997,]
tmp <- merge(tmp,AFSC_slope, all=T)


#  getting the darkblotched data -- years vary by survey
darkblotched <- tmp[tmp$common_name=="darkblotched rockfish" &
              ((tmp$project=="AFSC/RACE Slope Survey" & tmp$year==2001) |
              (tmp$project=="Groundfish Slope and Shelf Combination Survey") |
              (tmp$project=="Groundfish Slope Survey" & tmp$year==2001) |
              (tmp$project=="Groundfish Triennial Shelf Survey" & tmp$year %in% c(2001,2004))),]
              
#  getting the sablefish data -- years vary by survey
sablefish <- tmp[tmp$common_name=="sablefish" &
              ((tmp$project=="AFSC/RACE Slope Survey" & tmp$year==2001) |
              (tmp$project=="Groundfish Slope and Shelf Combination Survey") |
              (tmp$project=="Groundfish Slope Survey" & tmp$year %in% c(1998,2001)) |
              (tmp$project=="Groundfish Triennial Shelf Survey" & tmp$year %in% c(1998,2001,2004))),]

# getting trhe shallow species (+shortbelly) -- use all years from two surveys
othspp_names <- c("petrale sole","lingcod","Pacific hake","shortbelly rockfish","Pacific sanddab")
othspp_dat <-  tmp[tmp$common_name %in% othspp_names &
                        tmp$project %in% c("Groundfish Slope and Shelf Combination Survey", "Groundfish Triennial Shelf Survey"),] 

# merging data
WareHouse.All.Ages.Env_filtered_Feb2021 <- merge(darkblotched,sablefish, all=T)
WareHouse.All.Ages.Env_filtered_Feb2021 <- merge(WareHouse.All.Ages.Env_filtered_Feb2021,othspp_dat, all=T)

### get the depth strata

WareHouse.All.Ages.Env_filtered_Feb2021$depth_bin[WareHouse.All.Ages.Env_filtered_Feb2021$depth_ftm*1.8288 <= 184] <- 1 
WareHouse.All.Ages.Env_filtered_Feb2021$depth_bin[WareHouse.All.Ages.Env_filtered_Feb2021$depth_ftm*1.8288 > 184  & WareHouse.All.Ages.Env_filtered_Feb2021$depth_ftm*1.8288 <= 550] <- 2
WareHouse.All.Ages.Env_filtered_Feb2021$depth_bin[WareHouse.All.Ages.Env_filtered_Feb2021$depth_ftm*1.8288 > 550] <- 3
WareHouse.All.Ages.Env_filtered_Feb2021$depth_bin[is.na(WareHouse.All.Ages.Env_filtered_Feb2021$depth_ftm)] <- -9

# ### get the area strata (INPFC areas) 

WareHouse.All.Ages.Env_filtered_Feb2021$area[WareHouse.All.Ages.Env_filtered_Feb2021$latitude<=36] <- "Conception"
WareHouse.All.Ages.Env_filtered_Feb2021$area[WareHouse.All.Ages.Env_filtered_Feb2021$latitude>36 & WareHouse.All.Ages.Env_filtered_Feb2021$latitude<=40.5] <- "Monterey"
WareHouse.All.Ages.Env_filtered_Feb2021$area[WareHouse.All.Ages.Env_filtered_Feb2021$latitude>40.5 & WareHouse.All.Ages.Env_filtered_Feb2021$latitude<=43] <- "Eureka"
WareHouse.All.Ages.Env_filtered_Feb2021$area[WareHouse.All.Ages.Env_filtered_Feb2021$latitude>43 & WareHouse.All.Ages.Env_filtered_Feb2021$latitude<=47.5] <- "Columbia"
WareHouse.All.Ages.Env_filtered_Feb2021$area[WareHouse.All.Ages.Env_filtered_Feb2021$latitude>47.5] <- "Vancouver"

# ### get the area strata (INPFC areas), but combining the Eureka, columbia, and Vancouver areas into ECV 

WareHouse.All.Ages.Env_filtered_Feb2021$area2[WareHouse.All.Ages.Env_filtered_Feb2021$latitude<=36] <- "Conception"
WareHouse.All.Ages.Env_filtered_Feb2021$area2[WareHouse.All.Ages.Env_filtered_Feb2021$latitude>36 & WareHouse.All.Ages.Env_filtered_Feb2021$latitude<=40.5] <- "Monterey"
WareHouse.All.Ages.Env_filtered_Feb2021$area2[WareHouse.All.Ages.Env_filtered_Feb2021$latitude>40.5] <- "ECV"

# remove temporary objects
rm(WareHouse.All.Ages.Env_2, tmp, AFSC_slope, darkblotched, sablefish,
   othspp_names,othspp_dat)




