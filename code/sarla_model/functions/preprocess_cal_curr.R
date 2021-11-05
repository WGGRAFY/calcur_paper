#' preprocess_cal_curr
#'A function that takes the input raw data, selects a species' data from it, and filters it.
#'
#' @param data__ A Warehouse data object containing survey data
#' @param common_ Common name of the species you are interested in
#' @param sex_ sex to look at, must be "M" or "F"
#' @param years_ The year at which the survey type changes, this can be NULL if you want the data from all or one survey across all years
#' @param survey_string A string in the `project` field that specifies which survey
preprocess_cal_curr <- function(data__, common_, sex_, years_, survey_string){
  #Select species data, add year column
  full_data <- data__ %>%
    filter(common_name==common_) %>%
    filter(sex==sex_) %>%
    mutate(year = substr(datetime_utc_iso, 1, 4))

  #Add project/date filter
  if(!is.na(years_)){
    spp_data <- full_data %>%
      filter((str_detect(project, survey_string[1]) & (year < years_))|(str_detect(project, survey_string[2]) & (year >= years_)))
  } else{
    spp_data <- full_data %>% filter(str_detect(project, survey_string))
  }

  spp_data <- spp_data %>%
    select(age_years, year, length_cm, project) %>%
    mutate(year = as.numeric(year))
  return(spp_data)
}
