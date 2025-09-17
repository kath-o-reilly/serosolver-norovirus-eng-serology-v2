
resample_to_target_age_distribution <- function(inf_chain, age_distribution) {
  # This function takes an infection chain dataframe and an age distribution dataframe,
  # and resamples individuals within each (samp_no, chain_no, time) group to match
  # the target age distribution.
  #
  # Args:
  #   inf_chain: DataFrame with columns samp_no, chain_no, time, individual, age
  #   age_distribution: DataFrame with columns time, Age, weighted (target proportions)
  #
  # Returns:
  #   DataFrame with resampled individuals to match target age distribution.
  
  # 0) Target proportions by time & age
  target_prop <- age_distribution %>%
    transmute(time, age = as.integer(Age), target_prop = weighted)
  
  # 1) Start from the set of individuals alive at each (samp_no, chain_no, time)
  #    If inf_chain can have multiple rows per individual-time (e.g., events),
  #    keep one row per individual so we sample persons, not occurrences.
  pop <- inf_chain %>%
    distinct(samp_no, chain_no, time, individual, age, .keep_all = TRUE)
  
  # 2) Compute the current age distribution within each (samp_no, chain_no, time)
  curr_age <- pop %>%
    count(samp_no, chain_no, time, age, name = "n_age") %>%
    group_by(samp_no, chain_no, time) %>%
    mutate(N_group = sum(n_age),
           prop_curr = n_age / N_group) %>%
    ungroup()
  
  # 3) Join target proportions and compute re-weighting factors per age
  weights_by_age <- curr_age %>%
    left_join(target_prop, by = c("time", "age")) %>%
    mutate(
      target_prop = coalesce(target_prop, 0),
      # If current prop is 0, set weight to 0 (no individuals of that age to sample)
      w_age = if_else(prop_curr > 0, target_prop / prop_curr, 0)
    ) %>%
    select(samp_no, chain_no, time, age, N_group, w_age)
  
  # 4) Attach the age-level weight to each individual, then normalize within group
  pop_weighted <- pop %>%
    left_join(weights_by_age, by = c("samp_no", "chain_no", "time", "age")) %>%
    group_by(samp_no, chain_no, time) %>%
    # If all w_age are 0 (e.g., target has zero mass on all available ages or missing),
    # fall back to equal weights; otherwise normalize w_age.
    mutate(
      w_age = replace_na(w_age, 0),
      sum_w = sum(w_age),
      prob = if_else(sum_w > 0, w_age / sum_w, 1 / n())
    ) %>%
    select(-sum_w) %>%
    ungroup()
  
  # 5) Sample WITH replacement to keep the same total size per group
  #    (rows are repeated according to how many times that individual is drawn)
  resampled <- pop_weighted %>%
    group_by(samp_no, chain_no, time) %>%
    group_modify(function(df, keys) {
      size_out <- df$N_group[match(TRUE, !is.na(df$N_group))]  # same N as original group
      idx <- sample.int(n = nrow(df), size = size_out, replace = TRUE, prob = df$prob)
      # Return the selected rows; add draw index if you want
      df[idx, , drop = FALSE] %>%
        mutate(draw_id = row_number())
    }) %>%
    ungroup() %>%
    arrange(samp_no, chain_no, time, draw_id)
  
  resampled_check <- resampled %>% group_by(samp_no,chain_no,time,age) %>% tally() %>% group_by(samp_no,chain_no,time) %>% mutate(total_n=sum(n)) %>% group_by(time,age) %>% summarize(med_n = median(n)) %>% group_by(time) %>% mutate(prop=med_n/sum(med_n)) %>%
    mutate(prop = prop/max(prop))
  
  
  resampled <- resampled %>% left_join(resampled_check)
  
  ## Setup alternative age groups, 0-4, 5-10, and 5+
  resampled <- resampled %>% mutate(age_group_wide = if_else(age <= 4, "0-4","5-10"),
                                    age_group_simple = if_else(age < 5, as.character(age), "5+"))
  return(resampled)
}