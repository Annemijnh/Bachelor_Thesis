# Required packages
library(tidyverse)
library(data.table)
library(readr)
library(readxl)
library(fs)

# Import inclusion --------------------------------------------------------

# First, set working directory to one that includes the folder "DIALECT data"
# The first part of the script will check whether step and glucose data is in the same folder, if it is, it will be included for import. 
# It's possible to add diet info too, but for now I removed it.

file.list <- list.files(path=".", recursive = TRUE, pattern="\\glucose.txt")
activity <- list.files(path = '.', recursive= TRUE, pattern = '\\-steps.xlsx')
activity.list <- activity[lengths(strsplit(activity, "/")) < 4 ] # This makes sure that only the main directories are used
inclusion_cgm <- list()
inclusion_step <- list()
# The for loops below check whether the CGM and activity files are in the subject directory, and if so, include the full directory + filename to list. 
for (i in file.list){
  if (path_dir(i) %in% path_dir(activity.list)){
    inclusion_cgm <- append(inclusion_cgm,i)
  } 
}
for (i in activity.list){
  if (path_dir(i) %in% path_dir(file.list)){
    inclusion_step <- append(inclusion_step,i)
  }
}

# CGM data ----------------------------------------------------------------
# Importing only the included CGM data, setting it to a form that's useful for the other script
cgm_data <- lapply(inclusion_cgm, function(x) {
  read.table(x,sep ="", dec=",", header=T, na.strings = 'NA',fill = TRUE, skip = 2)|>
    select(c(2,3,5))|>  # only id, datetime, glucose value
    setnames(new = c('d', 't', 'g'))|>
    mutate(id = as.factor(gsub('[DIALECT data/]','',path_dir(x))) , .before='d',
           dx = parse_date_time(d,orders = c('dmy','ymd','mdy')),
           time = strptime(paste(dx,t),format = "%Y-%m-%d %H:%M", tz ='UTC'),
           time = time + 1, # R seems to mess up midnight and removes it from the datetime, so to fix it I add 1 second. Not using it in analysis anyway, since I mostly look at days, hours and minutes. 
           gl = as.numeric(g))|>
    select('id','time','gl')|>
    arrange(id,time)|>
    na.omit()
})

# This patient had their data in different columns
cgm_dataspecial <- read.table('DIALECT data/745/745-glucose.txt',sep ="", dec=",", header=T, na.strings = 'NA',fill = TRUE, skip = 2)|>
  select(c(4,5,7))|>
  setnames(new = c('d', 't', 'g'))|>
  mutate(id = as.factor(gsub('[DIALECT data/]','',path_dir('DIALECT data/745/745-glucose.txt'))) , .before='d',
         dx = parse_date_time(d,orders = c('dmy','ymd','mdy')),
         time = strptime(paste(dx,t),format = "%Y-%m-%d %H:%M", tz ='UTC'),
         time = time + 1, # R seems to mess up midnight and removes it from the datetime, so to fix it I add 1 second. Not using it in analysis anyway, since I mostly look at days, hours and minutes. 
         gl = as.numeric(g))|>
  select('id','time','gl')|>
  arrange(id,time)|>
  na.omit()

# These patients have extra data at the bottom that messes up the conversion of gl
extra_data <- c(633, 639, 654, 661, 674, 676, 677,678, 679, 680, 681, 682, 684, 685, 690, 691, 692, 696, 704, 705, 710, 714,734)
extra_data_idx <- list()
for (i in 1:length(cgm_data)){
  if (any(extra_data %in% gsub('[DIALECT data/]','',path_dir(inclusion_cgm[[i]])))){
    extra_data_idx <- append(extra_data_idx,i)
  }
}
extra_data_list <- unlist(inclusion_cgm[unlist(extra_data_idx)])
cgm_datadifferent <- lapply(extra_data_list, function(x) {
  read.table(x,sep ="", dec=",", header=T, na.strings = 'NA',fill = TRUE, skip = 2)|>
  select(c(2,3,5))|>
  head(-2)|>
  setnames(new = c('d', 't', 'g'))|>
  mutate(id = as.factor(gsub('[DIALECT data/]','',path_dir(x))) , .before='d',
           dx = parse_date_time(d,orders = c('dmy','ymd','mdy')),
           time = strptime(paste(dx,t),format = "%Y-%m-%d %H:%M", tz ='UTC'),
           time = time + 1, # R seems to mess up midnight and removes it from the datetime, so to fix it I add 1 second. Not using it in analysis anyway, since I mostly look at days, hours and minutes. 
           gl = as.numeric(gsub(',','.',g)))|>
  select('id','time','gl')|>
    arrange(id,time)|>
    na.omit()
})

# Seperate cleanup, some intervals need to be selected since they are viable. At the bottom of this file, I added the subject numbers to evaluate the data. 
# this is very specific and not very nice, but I dont see a way to automate it for now. Remove this bit for reproducability
cgm_data[[9]] <- cgm_data[[9]]|>tail(-2637)
cgm_data[[11]] <- cgm_data[[11]]|>head(16381) 
cgm_data[[12]] <- cgm_data[[12]]|>head(10718)
cgm_data[[17]] <- cgm_data[[17]]|>head(16695)
cgm_data[[21]] <- cgm_data[[21]]|>head(1063)
cgm_data[[24]] <- cgm_data[[24]]|>tail(-7681)
cgm_data[[26]] <- cgm_data[[26]]|>tail(-450)|> head(800)
cgm_data[[29]] <- cgm_data[[29]]|>tail(- 455)
cgm_data[[37]] <- cgm_data[[37]]|>head(1174)
cgm_data[[42]] <- cgm_data[[42]]|>head(1010)
cgm_data[[43]] <- cgm_data[[43]]|>head(811)
cgm_data[[44]] <- cgm_data[[44]]|>head(699)
cgm_data[[49]] <- cgm_data[[49]]|>tail(-501)
cgm_data[[59]] <- cgm_data[[59]]|>tail(-111)|> head(943)
cgm_data[[69]] <- cgm_data[[69]]|>tail(-85)
cgm_data[[73]] <- cgm_data[[73]]|>head(900)
cgm_data[[75]] <- cgm_data[[75]]|>tail(-383)
cgm_datadifferent[[10]] <- cgm_datadifferent[[10]]|> head(1308)
cgm_datadifferent[[11]] <- cgm_datadifferent[[11]]|> head(983)
cgm_data[[89]] <- cgm_data[[89]]|>head(988)
cgm_data[[92]] <- cgm_data[[92]]|>head(1003)
cgm_data[[99]] <- cgm_data[[99]]|>head(754)
cgm_data[[108]] <- cgm_data[[108]]|>head(815)
cgm_datadifferent[[19]] <- cgm_datadifferent[[19]]|> tail(-160)
cgm_data[[114]] <- cgm_data[[114]]|>head(999)
cgm_data[[116]] <- cgm_data[[116]]|>head(964)
cgm_data[[117]] <- cgm_data[[117]]|>head(908)
cgm_data[[124]] <- cgm_data[[124]]|>tail(-6)
cgm_data[[126]] <- cgm_data[[126]]|>tail(-62)|> head(718)
cgm_data[[130]] <- cgm_data[[130]]|>tail(-471)
cgm_data[[133]] <- cgm_data[[133]]|>tail(-72)
cgm_data[[134]] <- cgm_data[[134]]|>tail(-9899)|> head(728)
cgm_data[[152]] <- cgm_data[[152]]|>head(887)
cgm_data[[154]] <- cgm_data[[154]]|>tail(-112)|> head(764)
cgm_data[[155]] <- cgm_data[[155]]|>tail(-  34158)
cgm_data[[156]] <- cgm_data[[156]]|>head(758)

cgm_data <- bind_rows(cgm_data)
cgm_datadifferent <- bind_rows(cgm_datadifferent)
cgm_data <- rbind(cgm_data,cgm_datadifferent)
cgm_data <- rbind(cgm_data,cgm_dataspecial)

# Activity data -----------------------------------------------------------
# Before you proceed: id 752, 754, 755 have sums at the bottom? I removed manually, otherwise import gives NA 
# Same import here, importing and wrangling into a usable dataframe
activity_data <- lapply(inclusion_step,function(x){
  read_excel(x)|> 
    rename_at(1,~'time')|>
    mutate(id = as.factor(gsub('[DIALECT data/]','',path_dir(x))) , .before='time')
})
weird <- list()
strange <- list()
strange_dates <- c(745,746,748,749,750,751,752,753,754,755,756,757,758,759,760,761,762,765,766) # dates are also weird with this set
strange_data <- c(758,760,762,765,766) # this data was imported in a weird way (for some reason, it adds 1899-12-31 to the time), so for reproducability I search the subject id that needs adjustment
for (i in 1:length(activity_data)){
  if (any(strange_data %in% gsub('[DIALECT data/]','',path_dir(inclusion_step[[i]])))){
    weird <- append(weird,i)
  }
}
for (i in 1:length(activity_data)){
  if (any(strange_dates %in% gsub('[DIALECT data/]','',path_dir(inclusion_step[[i]])))){
    strange <- append(strange,i)
  }
}
# To fix the strange 1899-12-31 issue, this was set to just time with the loop below
for (i in weird){ 
activity_data[[i]]$time <- as.character(strftime(activity_data[[i]]$time, format="%H:%M:%S"))
}
# Columns indicate date, so change that to rows
for (i in 1:length(activity_data)){
activity_data[[i]] <- activity_data[[i]]|>
  pivot_longer(
    cols = -c(time,id),
    names_to = 'date',
    values_to = 'steps')
}
# Now dates are in rows, change the strange numerical dates to actual dates
for (i in strange){
  activity_data[[i]]$date <- as.character(as.Date(as.numeric(activity_data[[i]]$date), origin = "1899-12-30"))
}
# Now, parse all date and time as datetime
for (i in 1:length(activity_data)){
activity_data[[i]] <- activity_data[[i]]|>  
  mutate(date = parse_date_time(date,orders = c('dmy','ymd','mdy')),
         time = ymd_hms(paste(date,time),truncated = 1),
         time = time + 1)|> # To fix midnight measurements
  arrange(id,time)|>
  select(id,time,steps)
}
activity_data <- bind_rows(activity_data)

# Selection ---------------------------------------------------------------
# Check intervals of each dataset
activity_check <- activity_data|>
  group_by(id)|>
  mutate(start_date = min(time),
         end_date = max(time),
         interval_act = interval(start_date,end_date))|>
  select(id,interval_act)
activity_check <- distinct(activity_check)

cgm_check <- cgm_data|>
  na.omit()|>
  group_by(id)|>
  mutate(start_date = min(time),
         end_date = max(time),
         interval_cgm = interval(start_date,end_date))|>
  select(id,interval_cgm)
cgm_check <- distinct(cgm_check)

# Determine interval, determine overlap
df <- merge(cgm_check,activity_check, by = 'id')|>
  mutate(overlap = int_overlaps(interval_cgm,interval_act),
         overlap_start_date = case_when(int_start(interval_cgm) > int_start(interval_act) ~ int_start(interval_cgm),
                                        int_start(interval_act) > int_start(interval_cgm) ~ int_start(interval_act),
                                        TRUE~ int_start(interval_cgm)),
         overlap_end_date = case_when(int_end(interval_cgm) < int_end(interval_act) ~ int_end(interval_cgm),
                                      int_end(interval_act) < int_end(interval_cgm) ~ int_end(interval_act),
                                      TRUE ~ int_end(interval_cgm)),
         overlap_interval = interval(overlap_start_date,overlap_end_date),
         length_interval = int_length(overlap_interval)/60/60/24)

# Select overlap
selection <- df|>
  select(id, overlap,overlap_interval,length_interval)

# Merge
cgm_data_selection <- merge(cgm_data,selection,by = 'id')
cgm_data_selection <- subset(cgm_data_selection,overlap == TRUE)|>
  na.omit()|>
  group_by(id)|>
  group_split()
# Select the cgm data that is in the overlap
for (i in 1:length(cgm_data_selection)){
  cgm_data_selection[[i]] <- cgm_data_selection[[i]]|>
    filter(time %within% overlap_interval)
}
cgm_data_selection <- bind_rows(cgm_data_selection)|>
  filter(length_interval > 7) # Only intervals longer than 7 days are accepted
length(unique(cgm_data_selection$id))
cgm_data_selection <- cgm_data_selection|>
  group_by(id)|>
  arrange(id,time)|>
  mutate(diff = as.double(time - lag(time),units = 'hours'),
         max_diff = max(diff,na.rm = TRUE)) #Determine maximum gaps in hours

# The bit below shows the inclusion dependent on the maximum gaps. 
cgm_fun <- list()
for (i in 1:24){
  cgm_fun[[i]] <- subset(cgm_data_selection, max_diff < i)
}
cgm_inc <- list()
cgm_sub_inc <- list()
for (i in 1:24){
  cgm_inc <- append(cgm_inc,length(unique(cgm_fun[[i]]$id)))
  cgm_sub_inc <- append(cgm_sub_inc,list(unique(cgm_fun[[i]]$id)))
}

# The choice was made to do a gap of 4 hours max
cgm_selection <- subset(cgm_data_selection, max_diff < 4)|>
  select(id, time, gl)|>
  mutate(time = as.POSIXct(paste(time,'+0100'),format="%F %T %z", tz="Europe/Amsterdam")) # Setting UTC to CET/CEST, since it has daylight savings
inclusion <- unlist(unique(cgm_selection$id)) # List of ids that are included, needed later for subject data import

# Select activity data of ids that were determined previously
act_data_selection <- merge(activity_data,selection,by = 'id')
act_data_selection <- subset(act_data_selection,id %in% unique(cgm_selection$id))|>
  na.omit()|>
  group_by(id)|>
  group_split()

df_included <- subset(df, id %in% inclusion)|>
  select(id, length_interval)
# Only select time within interval
for (i in 1:length(act_data_selection)){
  act_data_selection[[i]] <- act_data_selection[[i]]|>
    filter(time %within% overlap_interval)
}
act_data_selection <- bind_rows(act_data_selection)|>
  select(id,time,steps)|>
  mutate(time = as.POSIXct(paste(time,'+0100'),format="%F %T %z", tz="Europe/Amsterdam")) # Setting UTC to CET/CEST, since it has daylight savings

# Write to csv ------------------------------------------------------------
# Write to files, to be used in other script
# Or save the workspace, or just run the other scripts after this one.
# write.table(cgm_selection, file = 'cgm_data.csv',sep =',',row.names=FALSE)
# write.table(act_data_selection, file = 'activity_data.csv',sep = ',',row.names = FALSE)
# write.csv(inclusion,file='inclusion.csv',row.names = FALSE)


# Evaluate excluded cgm data ----------------------------------------------
# This makes a list of the ids that are not included 
# sub_check <- unique(cgm_data$id)[which(!(unique(cgm_data$id) %in% cgm_sub_inc[[4]]))]

# seperate_cgm <- subset(cgm_data_selection,id ==  745)|> # This one I manually adjusted w.r.t. id, and the tail() and head() to see how much data to remove from the beginning and end. 
#   na.omit()|>
#   arrange(time)|>
#   # tail(-  160)|>
#   # head(657)|>
#   mutate(diff = as.double(time - lag(time),units = 'hours'),
#          max_diff = max(diff,na.rm = TRUE),
#          length = as.double(max(time)-min(time),units = 'days'))
#   
# plot(seperate_cgm$time,seperate_cgm$gl) # Plot the cgm data to visually inspect gaps

## EXCLUDED
# pt 1002: multiple gaps over 4 hours
# pt 595: removing large gaps leaves little data
# pt 596: removing gaps leaves too little data
# pt 597: too short
# pt 600: many large gaps in middle of measurements
# pt 605: removing gaps leaves too little data
# pt 608: removing gaps leaves too little data
# pt 614: too many gaps in middle of data
# pt 617: many gaps throughout data
# pt 619: too little data left after removing gap
# pt 620: many gaps through data
# pt 622: removing the largest gap did not help much, many gaps exist still
# pt 623: too little data after removal
# pt 625: removing large gaps still shows some gaps over 8h in the middle of the set
# pt 626: many gaps throughout data
# pt 627: removal leaves too little data
# pt 630: removal leaves too little data
# pt 631: removal leaves too little data
# pt 632: removal leaves too little data
# pt 637: removal leaves too little data
# pt 642: removing large gap leaves a 6h gap in middle of data
# pt 645: many gaps throughout data
# pt 646: gap in middle of data
# pt 647: many gaps throughout data
# pt 652: removing large gap leaves 5h gaps in middle of data
# pt 653: removing large gap makes length interval too short
# pt 654:many gaps throughout data
# pt 656: removing large gap leaves little data
# pt 657: removing large gap leaves too short intervals
# pt 659: many gaps throughout data
# pt 661: many gaps throughout data
# pt 663: many gaps throughout data
# pt 667: many gaps throughout data
# pt 670: many gaps throughout data
# pt 671: too little data left after removals
# pt 672: many gaps throughout data
# pt 673: too short
# pt 674: many gaps throughout data
# pt 685: too little data left after removal
# pt 689: too little data left after removal
# pt 690: too little data left after removal
# pt 691: large gap
# pt 692: too little data left after removal
# pt 694: too little data left after removal
# pt 696: same
# pt 697 : too little data left after removal
# pt 698: same
# pt 699: interval too short
# pt 700 : many gaps throughout data
# pt 703: too little data left after removal  
# pt 706: many gaps throughout data
# pt 710: too many gaps
# pt 715: many gaps throughout data
# pt 728: many gaps throughout data
# pt 733:many gaps throughout data
# pt 734: too little data left after removal  
# pt 745: too many gaps throughout data
# pt 746:many gaps throughout data
# pt 751: no overlap
# pt 752: many gaps throughout data
# pt 754: removal leaves too little data
# pt 755: gaps throughout data
# pt 765: many gaps throughout data
# pt 766: many gaps throughout data

## INCLUDED: need these adjustments 
# pt 598: gap removal leaves enough data: PASS tail(-2637)|>
# pt 601: gap removal leaves enough data: PASS head(16381)|>
# pt 604: gap removal leaves enough data: PASS head(10718)|>
# pt 609: PASS head(16695)|>
# pt 613: PASS head(1063)|>
# pt 616: PASS tail(-7681)|>
# pt 618: PASS   tail(-450)|> head(800)|>
# pt 621: PASS   tail(- 455)|>
# pt 629: PASS  head(1174)|>
# pt 634: PASS head(1010)|>
# pt 635: PASS head(811)|>
# pt 636: PASS head(699)|>
# pt 641: PASS tail(-501)|>
# pt 651: PASS   tail(-111)|> head(943)|>
# pt 662: CAN PASS -- need to remove first 85 values (tail(-85))
# pt 666: PASS  head(900)|>
# pt 668: PASS   tail(-383)|>
# pt 680: PASS   head(1308)|> 
# pt 681: PASS   head(983)|>
# pt 683: PASS   head(988)|>
# pt 686: PASS   head(1003)|>
# pt 693: PASS   head(754)|>
# pt 702: PASS head(815)|>
# pt 704: PASS   tail(-  160)|>
# pt 709: PASS head(999)|>
# pt 711: PASS head(964)|>
# pt 712: PASS   head(908)|>
# pt 719: PASS tail(- 6)|>
# pt 722: PASS   tail(- 62)|> head(718)|>
# pt 726: PASS   tail(- 471)|>
# pt 729: PASS   tail(-  72)|>
# pt 730: PASS   tail(-  9899)|> head(728)|>
# pt 757: PASS   head(887)|>
# pt 759: PASS   tail(-  112)|> head(764)|>
# pt 760: PASS   tail(-  34158)|>
# pt 761: PASS   head(758)|>
