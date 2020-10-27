#function to determine detected pathogens in BAL data and code infections accordingly

code_pathogens = function(bal_data)
{
  require(tidyverse)
  #fix test columns
  test_cols = colnames(bal_data)[c(which(colnames(bal_data) == "ASPERGILLUS_GALACTOMANNAN_ANTIGEN_NMH_LFH"):which(colnames(bal_data) == "RESPIRATORY_SYNCYTIAL_VIRUS_RESPAN22"),
                              which(colnames(bal_data) == "STREPTOCOCCUS_PNEUMONIAE_ANTIGEN_URINE_1"),
                              which(colnames(bal_data) == "LEGIONELLA_ANTIGEN__EIA__URINE_1"))]
  
  bal_data = bal_data %>% mutate_at(.vars = test_cols, 
                          .funs = function(x){
                            x = factor(ifelse(is.na(x),
                                              yes = NA,
                                              no = ifelse((grepl("Not Detected", x, ignore.case = T) |
                                                             grepl("Negative", x, ignore.case = T)),
                                                          yes = "Negative",
                                                          no = "Positive")))
                            return(x) })
  
  # fix culture data 
  bal_data = bal_data %>% 
    pivot_longer(cols = contains("organism_"),
                 names_to = "microbiology_parameter",
                 values_to = "microbiology_value") %>% 
    mutate(microbiology_value = trimws(microbiology_value))

  bal_data$main_microbiology_parameter = factor(gsub("_*\\d", "", bal_data$microbiology_parameter))
  #flatten these parameters into lists of values
  bal_data = bal_data %>%
    group_by(ir_id, BAL_order_timestamp, main_microbiology_parameter) %>%
    mutate(microbiology_value_list = list(microbiology_value)) %>%
    ungroup() %>% 
    rowwise() %>% 
    mutate_at(.vars = "microbiology_value_list", .funs = function(x){ #remove NA
      cur = na.omit(x)
      if(length(cur) == 0)
      {
        return(list(NULL))
      } else
      {
        return(list(cur))
      }
    }) %>%
    select(-c(microbiology_parameter, microbiology_value)) %>% 
    ungroup()
  
  #pivot back into wide form for merging (list values get duplicated, need to fix)
  listcols = as.character(unique(bal_data$main_microbiology_parameter))
  bal_data = bal_data %>%
    pivot_wider(names_from = main_microbiology_parameter,
                values_from = microbiology_value_list) %>% 
    rowwise() %>% 
    #have to remove duplicated list vals
    mutate_at(.vars = listcols, .funs = function(x){ 
      return(list(x[[1]])) }) %>% 
    ungroup()
  bal_data$organism_quantity = lapply(bal_data$organism_quantity,
                                 function(x){
                                   x = gsub(">", "", x)
                                   x = gsub("<.+", "0", x)
                                   x = gsub(",", "", x)
                                   x = as.numeric(x)
                                   if(length(x) == 0)
                                   {
                                     return(NULL)
                                   }
                                   return(x)})
  
  #remove commensal
  bal_data = dplyr::select(bal_data, - Organism_ID) #really not needed
  all_pathogens = unique(unlist(bal_data$organism_name))
  commensal = c("Yeast, Not Cryptococcus Species", "Corynebacterium species", 
                "Stomatococcus species", "Yeast, Not Cryptococcus Species",
                "Lactobacillus species", "Yeast, Not Cryptococcus Species #2",
                "Staphylococcus coagulase negative", "Lactobacillus species #2", 
                all_pathogens[grepl("Candida", all_pathogens)])
  for(i in 1:nrow(bal_data))
  {
    cur_name = bal_data$organism_name[i][[1]]
    cur_quantity = bal_data$organism_quantity[i][[1]]
    if(is.null(cur_name))
    {
      next
    }
    #include "Stomatococcus species" if it's over 1000
    good_indices = which(!(cur_name %in% commensal) | 
                           (cur_name == "Stomatococcus species" & cur_quantity >= 1000))
    if(is.null(cur_name[good_indices]))
    {
      bal_data$organism_name[i][[1]] = NULL
      bal_data$organism_quantity[i][[1]] = NULL
    } else
    {
      bal_data$organism_name[i][[1]] = cur_name[good_indices]
      bal_data$organism_quantity[i][[1]] = cur_quantity[good_indices]
    }
  }
  
  # bin into broad categories   
  bal_data$any_iav = (bal_data$INFLUENZA_A == "Positive" | 
                   bal_data$INFLUENZA_A_RESPAN23 == "Positive" |
                   bal_data$INFLUENZA_A_H1_2009 == "Positive")
  bal_data$any_iav[is.na(bal_data$any_iav)] = FALSE
  
  bal_data$any_ibv = (bal_data$INFLUENZA_B == "Positive" | 
                   bal_data$INFLUENZA_B_RESPAN21 == "Positive")
  bal_data$any_ibv[is.na(bal_data$any_ibv)] = FALSE
  
  bal_data$any_influenza = (bal_data$INFLUENZA_A == "Positive" | 
                         bal_data$INFLUENZA_A_RESPAN23 == "Positive" |
                         bal_data$INFLUENZA_A_H1_2009 == "Positive" | 
                         bal_data$INFLUENZA_B_RESPAN21 == "Positive" |
                         bal_data$INFLUENZA_B == "Positive" |
                         bal_data$INFLUENZA_B_RESPAN21 == "Positive")
  bal_data$any_influenza[is.na(bal_data$any_influenza)] = FALSE
  
  bal_data$any_bacterial = (bal_data$ENTEROBACTER_AEROGENES == "Positive" |
                         bal_data$ENTEROBACTER_CLOACAE_COMPLEX == "Positive" |
                         bal_data$ESCHERICHIA_COLI == "Positive" |
                         bal_data$HAEMOPHILUS_INFLUENZAE == "Positive" |
                         bal_data$KLEBSIELLA_AEROGENES == "Positive" |
                         bal_data$KLEBSIELLA_PNEUMONIAE_GROUP == "Positive" |
                         bal_data$PSEUDOMONAS_AERUGINOSA == "Positive" |
                         bal_data$SERRATIA_MARCESCENS == "Positive" |
                         bal_data$STAPHYLOCOCCUS_AUREUS == "Positive" |
                         bal_data$STREPTOCOCCUS_AGALACTIAE == "Positive" |
                         bal_data$STREPTOCOCCUS_PNEUMONIAE_ANTIGEN_URINE_1 == "Positive" |
                         bal_data$STREPTOCOCCUS_PNEUMONIAE == "Positive" |
                         bal_data$STREPTOCOCCUS_PYOGENES == "Positive" |
                         bal_data$MECA_C_AND_MREJ == "Positive" |
                         bal_data$ACINETOBACTER_CALCOACETICUS_BAUMANNII_COMPLEX == "Positive" |
                         bal_data$KLEBSIELLA_OXYTOCA == "Positive" |
                         bal_data$PROTEUS_SPP == "Positive" |
                         bal_data$SERRATIA_MARCESCENS == "Positive" |
                         bal_data$PSEUDOMONAS_AERUGINOSA == "Positive" |
                         bal_data$HAEMOPHILUS_INFLUENZAE == "Positive" |
                         bal_data$CHLAMYDIA_PNEUMONIAE_lower_resp == "Positive" |
                         bal_data$MYCOPLASMA_PNEUMONAIE == "Positive" |
                         bal_data$MORAXELLA_CATARRHALIS == "Positive" |
                         bal_data$LEGIONALLA_PNEUMOPHILIA == "Positive" |
                         bal_data$BORDETELLA_PARAPERTUSSIS == "Positive" |
                         bal_data$BORDETELLA_PERTUSSIS == "Positive" |
                         bal_data$CHLAMYDIA_PNEUMONIAE == "Positive" |
                         bal_data$MYCOPLASMA_PNEUMONIAE == "Positive" |
                         bal_data$STREPTOCOCCUS_PNEUMONIAE_ANTIGEN_URINE_1 == "Positive")
  bal_data$any_bacterial[is.na(bal_data$any_bacterial)] = FALSE
  
  bal_data$any_parainfluenza = (bal_data$PARAINFLUENZAE_VIRUS == "Positive" |
                             bal_data$PARA_INFLU_1 == "Positive" |
                             bal_data$PARA_INFLU_2 == "Positive" |
                             bal_data$PARA_INFLU_3 == "Positive" |
                             bal_data$PARAINFLUENZA_VIRUS_4 == "Positive")
  bal_data$any_parainfluenza[is.na(bal_data$any_parainfluenza)] = FALSE
  
  bal_data$any_coronavirus_non = (bal_data$CORONAVIRUS_HKU1 == "Positive" |
                               bal_data$CORONAVIRUS_229E == "Positive" |
                               bal_data$CORONAVIRUS_NL63 == "Positive" |
                               bal_data$CORONAVIRUS_OC43 == "Positive")
  bal_data$any_coronavirus_non[is.na(bal_data$any_coronavirus_non)] = FALSE
  
  bal_data$any_adenovirus = (bal_data$ADENOVIRUS == "Positive")
  bal_data$any_adenovirus[is.na(bal_data$any_adenovirus)] = FALSE
  
  bal_data$any_rsv = (bal_data$RESPIRATORY_SYNCYTIAL_VIRUS == "Positive" |
                   bal_data$RESPIRATORY_SYNCYTIAL_VIRUS_RESPAN22 == "Positive")
  bal_data$any_rsv[is.na(bal_data$any_rsv)] = FALSE
  
  bal_data$any_metapneumonia = (bal_data$HUMAN_METAPNEUMOVIRUS_lower_resp == "Positive" |
                             bal_data$HUMAN_METAPNEUMOVIRUS == "Positive")
  bal_data$any_metapneumonia[is.na(bal_data$any_metapneumonia)] = FALSE
  
  bal_data$any_rhinovirus_enterovirus = (bal_data$HUMAN_RHINOVIRUS_ENTEROVIRUS_lower_resp == "Positive" |
                                      bal_data$HUMAN_RHINOVIRUS_ENTEROVIRUS == "Positive")
  bal_data$any_rhinovirus_enterovirus[is.na(bal_data$any_rhinovirus_enterovirus)] = FALSE
  
  bal_data$any_viral = (bal_data$any_influenza |
                     bal_data$any_adenovirus |
                     bal_data$any_parainfluenza |
                     bal_data$any_coronavirus_non |
                     bal_data$any_rsv |
                     bal_data$any_metapneumonia |
                     bal_data$any_rhinovirus_enterovirus)
  bal_data$any_viral[is.na(bal_data$any_viral)] = FALSE
  
  bal_data$any_viral_non_covid = (bal_data$any_influenza |
                         bal_data$any_adenovirus |
                         bal_data$any_parainfluenza |
                         bal_data$any_coronavirus_non |
                         bal_data$any_rsv |
                         bal_data$any_metapneumonia |
                         bal_data$any_rhinovirus_enterovirus)
  bal_data$any_viral_non_covid[is.na(bal_data$any_viral_non_covid)] = FALSE
  
  bal_data = bal_data %>%
    mutate(any_culture = lengths(organism_name) > 0) #length(NULL) = 0
  
  #for general "Stomatococcus species", keep only if N >= 1000
  for(i in 1:nrow(bal_data))
  {
    organisms = unlist(bal_data$organism_name[i])
    quantities = unlist(bal_data$organism_quantity[i])
    if(is.null(organisms))
    {
      next
    }
    if(length(organisms) == 1 && organisms[1] == "Stomatococcus species")
    {
      if(quantities[1] < 1000)
      {
        bal_data$organism_name[i] = list(NULL)
        bal_data$any_culture[i] = FALSE
      }
    }
  }
  
  bal_data$any_nonviral = (bal_data$any_bacterial |
                        bal_data$any_culture)
  bal_data$any_nonviral[is.na(bal_data$any_nonviral)] = FALSE
  
  return(bal_data)
}