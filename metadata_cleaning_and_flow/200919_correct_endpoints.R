#function to determine the true number of days on a ventilator
#counts only days when intubated   

correct_vent_days = function(endpoints, bal_data)
{
  require(tidyverse)
  #get colnames into standard form
  colnames(bal_data) = gsub("\\.", "_", colnames(bal_data)) #I like underscores
  colnames(bal_data) = gsub("_+$", "", colnames(bal_data)) # remove trailing
  colnames(bal_data) = gsub("^_+", "", colnames(bal_data)) # remove leading
  colnames(bal_data) = gsub("_+", "_", colnames(bal_data)) #remove dup underscores
  colnames(endpoints) = gsub("\\.", "_", colnames(endpoints))
  colnames(endpoints) = gsub("_+$", "", colnames(endpoints))
  colnames(endpoints) = gsub("^_+", "", colnames(endpoints)) 
  colnames(endpoints) = gsub("_+", "_", colnames(endpoints))
  
  #get endpoints into long-form, if necessary
  is_long_form = !("First_intub_start" %in% colnames(endpoints))
  if(is_long_form == FALSE)
  {
    endpoints = endpoints %>% 
      pivot_longer(cols = c(First_intub_start, First_intub_stop,
                            Second_intub_start, Second_intub_stop,
                            Third_intub_start, Third_intub_stop),
                   names_to = c("intub_index", "intub_action"),
                   names_pattern = "(.+)_intub_(.+)",
                   values_to = "intub_action_date", 
                   values_drop_na = TRUE) %>% 
      dplyr::mutate(intub_index = as.numeric(case_when(intub_index == "First" ~ 1,
                                                 intub_index == "Second" ~ 2,
                                                 intub_index == "Third" ~ 3))) %>% 
      #now make one row for each intubation
      pivot_wider(names_from = "intub_action",
                  values_from = "intub_action_date",
                  names_prefix = "intub_")
  }
  
  #remove times because accuracy is not good enough
  endpoints = mutate(endpoints, 
                     intub_stop = as.Date(intub_stop),
                     intub_start = as.Date(intub_start),
                     intub_duration = difftime(intub_stop, intub_start, units = "days")) 
  bal_data = mutate(bal_data, BAL_order_date = as.Date(BAL_order_date))
  
  #remove any duplicate entries
  endpoints = unique(endpoints)
  bal_data = unique(bal_data)
  
  corrected_vent_days = vector(mode = "numeric", length = nrow(bal_data))
  intub_index = vector(mode = "numeric", length = nrow(bal_data))
  between_intubs = vector(mode = "logical", length = nrow(bal_data))
  for(i in 1:nrow(bal_data))
  {
    cur_pt = bal_data[i, "ir_id", drop = T]
    cur_bal_date = bal_data[i, "BAL_order_date", drop = T]
    endpoints_sub = dplyr::filter(endpoints, ir_id == cur_pt)
    #find vent period this BAL falls in, get day
    current_vent_period = dplyr::filter(endpoints_sub, intub_start <= cur_bal_date & 
                                          intub_stop >= cur_bal_date)
    if(nrow(current_vent_period) > 0)
    {
      current_vent_period_id = current_vent_period$intub_index
      days_in_vent_period = as.numeric(difftime(cur_bal_date, current_vent_period$intub_start, units = "days"))
      between_intubs[i] = FALSE
    } else #BAL between vent periods, find last completed
    {
      current_vent_period = dplyr::filter(endpoints_sub, intub_stop <= cur_bal_date)
      current_vent_period_id = ifelse(nrow(current_vent_period) > 0,
                                      yes = max(current_vent_period$intub_index) + 1, #so we look at previous periods
                                      no = 0)
      days_in_vent_period = 0 #since we haven't entered it yet
      between_intubs[i] = TRUE
    }
    intub_index[i] = current_vent_period_id
    if(length(current_vent_period_id) > 1)
    {
      stop(paste("Error: duplicated patient / BAL pair at row ", i))
    }
    
    #add sum of previous vent periods
    if(current_vent_period_id > 1)
    {
      prior_vent_periods = dplyr::filter(endpoints_sub, intub_index < current_vent_period_id)
      prior_vent_days = sum(prior_vent_periods$intub_duration)
    } else
    {
      prior_vent_days = 0
    }
    
    corrected_vent_days[i] = as.numeric(days_in_vent_period + prior_vent_days)
  }
  
  bal_data$corrected_vent_days = corrected_vent_days
  bal_data$intub_index = intub_index
  bal_data$between_intubs = between_intubs
  
  #add reliable flag for true 48hr, first intubation (and first of those meeting these criteria)
  first_bals = bal_data %>% 
    group_by(ir_id) %>% 
    slice(which.min(BAL_order_date)) %>% 
    .$BAL_order_accession_num
  
  bal_data = bal_data %>% 
    mutate(true_early = (corrected_vent_days >= 0 &
                           corrected_vent_days <= 2 &
                           intub_index == 1 &
                           between_intubs == FALSE &
                           BAL_order_accession_num %in% first_bals))
  
  out = list(bal = bal_data, 
             endpoints = endpoints)
  return(out)
}