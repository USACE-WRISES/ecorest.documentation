# The purpose of this script is to test the functionality of ecorest functions to ensure they are working as expected
# Authors: Kiara Cushway
# Date: January 14th, 2026

## Clear workspace
rm(list = ls()) ## remove stored files and objects
gc(T) ## garbage collection
graphics.off() ## Turn off graphics

## Load library
library(ecorest)
library(backports)
library(stringr)
library(viridis)

################################################################################

# HSI functions

## Set seed for repeatability
set.seed(342)

## Create a function to test HSI-based functions in ecorest
HSI_tester = function(function_name, test_type = 'expected', iterations = 9999) {
  passes = 0 # set number of passes to zero initially
  ## Make sure the function name is correct
  if (!function_name %in% c('HSImin', 'HSIarimean', 'HSIgeomean', 'HSIwarimean')) {
    stop("Error: invalid function name. Please enter one of the following: 'HSImin', 'HSIarimean', 'HSIgeomean', 'HSIwarimean'.", call. = F)
  }
  ## Make sure the test_type is correct
  if (!test_type %in% c('NA', 'invalid', 'expected', 'invalid_weights', 'invalid_num_weights', 'dataframe')) {
    stop("Error: invalid test type. Please enter one of the following: 'NA', 'invalid', 'expected', 'invalid_weights', 'invalid_num_weights', 'dataframe'.", call. = F)
  }
  ## Make sure weighting is only compatible with HSIwarimean
  if (!function_name %in% 'HSIwarimean' & test_type %in% c('invalid_weights', 'invalid_num_weights')) {
    stop("Error: Invalid weights and invalid number of weights can only be tested with the HSIwarimean function.", call. = F)
  }
  ## Run specified number of iterations to test model
  for (i in 1:iterations) {
    num_inputs = sample(3:20, 1) # randomly generate number of inputs between 2 and 20
    input_vals = runif(n = num_inputs, min = 0, max = 1) # randomly generate specified number of input values
    num_NA = 0 # specify no NA values
    invalid = 0 # specify no invalid values
    ## If the test is testing inputs including NA values, randomly replace some values with NA
    if (test_type == 'NA') {
      num_NA = sample(1:(num_inputs - 2), 1) # randomly select number of NA values to be included
      input_vals[sample(length(input_vals), num_NA)] = NA # randomly replace some values with NA
    } else if (test_type == 'invalid') {
      # Inputs < 0 or > 1
        num_invalid = sample(1:(num_inputs - 2), 1) # randomly select number of invalid inputs to be included  
        invalid_type = sample(1:3, 1)
        if (invalid_type == 1) {
        invalid = sample(
          runif(1000, -10.4, 10.4) |> (\(x) x[x < 0 | x > 1])(),
          size = num_invalid,
          replace = TRUE
        ) # randomly generate invalid values
        } else if (invalid_type == 2) {
          # Infinite inputs
          invalid = rep(Inf, num_invalid)
        } else if (invalid_type == 3) {
          # Character inputs
          invalid = sample(letters, num_invalid)
        }
        
        input_vals[sample(length(input_vals), num_invalid)] = invalid 
    } else if (test_type == 'dataframe') {
      input_vals = data.frame(input_vals)
    }
    ## If the function being tested is HSImin, test limiting factor equation
    if (function_name == 'HSImin') {
      HSI = suppressWarnings(try(HSImin(input_vals), silent = T)) # calculate HSI
      ## Compare HSImin result with minimum value equation
      if (isTRUE(all.equal(HSI, min(input_vals, na.rm = T)))) {
        passes = passes + 1 # if values match, test passes
      } 
    } 
    ## If the function being tested is HSIarimean, compare arithmetic mean equation
    else if (function_name == 'HSIarimean') {
      HSI = suppressWarnings(try(HSIarimean(input_vals), silent = T)) # calculate HSI
      ## Compare HSIarimean to arithmetic mean output
      if (all(is.numeric(as.matrix(input_vals)))) {
        if (isTRUE(all.equal(HSI, suppressWarnings(sum(input_vals, na.rm = T) / (num_inputs - num_NA))))) {
          passes = passes + 1 # if values match, test passes
        }
      }
    } 
    ## If the function being tested is HSIgeomean, compare to geometric mean equation
    else if (function_name == 'HSIgeomean') {
      HSI = suppressWarnings(try(HSIgeomean(input_vals), silent = T)) # calculate HSI
      ## Compare HSIgeomean to geometric mean equation
      if (all(is.numeric(as.matrix(input_vals)))) {
        if (isTRUE(all.equal(HSI, (prod(input_vals, na.rm = TRUE) ^ (1 / (num_inputs - num_NA)))))) {
          passes = passes + 1 # if values match, test passes
        } 
      }
    
    }
    ## If the function being tested is HSIwarimean, compare to weighted arithmetic mean
    else if (function_name == 'HSIwarimean') {
      weights_raw = runif(num_inputs) # randomly generate weights
      ## If testing invalid weights, leave weights as raw so they don't add to one
      if (test_type == 'invalid_weights') {
        weights = weights_raw
      } 
      ## If testing incorrect number of weights, randomly add between 1 and 5 extra weights
      else if (test_type == 'invalid_num_weights') {
        weights = c(weights_raw, runif(sample(1:5, 1), min = 0, max = 1)) # add extra weights
        weights = weights / sum(weights) # normalize to sum to one
        weights[length(weights)] <- weights[length(weights)] + (1 - sum(weights)) # ensure that values sum to one
      } 
      ## If testing expected or NA values, calculate weights normally
      else {
        weights <- weights_raw
        
        weights[is.na(input_vals)] <- NA
        
        weights <- weights / sum(weights, na.rm = TRUE)
        
        weights[max(which(!is.na(weights)))] <-
          weights[max(which(!is.na(weights)))] + (1 - sum(weights, na.rm = TRUE))

      }
      if (test_type == 'dataframe') {
        df = data.frame(input_vals = input_vals, weights = weights)
        HSI = suppressWarnings(try(HSIwarimean(df$input_vals, df$weights), silent = T)) # calculate HSI
        ## Compare HSIwarimean to weighted arithmetic mean
        if (all(is.numeric(as.matrix(df$input_vals)))) {
          if (isTRUE(all.equal(HSI, suppressWarnings(sum(df$input_vals * df$weights, na.rm = TRUE)), tolerance = 1e-12))) {
            passes = passes + 1 # if values match, test passes
          }
        }
      } else {
        HSI = suppressWarnings(try(HSIwarimean(input_vals, weights), silent = T)) # calculate HSI
        ## Compare HSIwarimean to weighted arithmetic mean
        if (all(is.numeric(as.matrix(input_vals)))) {
          if (isTRUE(all.equal(HSI, suppressWarnings(sum(input_vals * weights, na.rm = TRUE)), tolerance = 1e-12))) {
            passes = passes + 1 # if values match, test passes
          }
        }
      }
    } 
  }
  # Pass rate
  if (test_type %in% c('expected', 'NA', 'dataframe')) {
  return(passes / iterations * 100) # if results should be valid, return the pass rate
  } else {
    return((iterations - passes) / iterations *100) # if results should be errors, return the fail rate
  }
}

## Apply function
HSI_tester('HSImin', 'expected', 9999)
HSI_tester('HSImin', 'NA', 9999)
HSI_tester('HSImin', 'invalid', 9999)
HSI_tester('HSImin', 'dataframe', 9999)

HSI_tester('HSIarimean', 'expected', 9999)
HSI_tester('HSIarimean', 'NA', 9999)
HSI_tester('HSIarimean', 'invalid', 9999)
HSI_tester('HSIarimean', 'dataframe', 9999)

HSI_tester('HSIgeomean', 'expected', 9999)
HSI_tester('HSIgeomean', 'NA', 9999)
HSI_tester('HSIgeomean', 'invalid', 9999)
HSI_tester('HSIgeomean', 'dataframe', 9999)

HSI_tester('HSIwarimean', 'expected', 9999)
HSI_tester('HSIwarimean', 'NA', 9999)
HSI_tester('HSIwarimean', 'invalid', 9999)
HSI_tester('HSIwarimean', 'invalid_weights', 9999)
HSI_tester('HSIwarimean', 'invalid_num_weights', 9999)
HSI_tester('HSIwarimean', 'dataframe', 9999)

################################################################################

# HSIeqtn
HSIeqtn_tester = function(test_type = 'expected', iterations = 100){
  passes = 0 # set initial pass rate to zero
  count_NA = 0 # set number of NA columns to zero
  num_invalid = 0 # set number of invalid inputs to zero
  num_extra = 0 # set number of extra values to zero
  exclude = NULL
  ## Make sure the test_type is correct
  if (!test_type %in% c('expected', 'exclude', 'invalid_name', 'invalid_SIV', 'wrong_num_SIV', 'dataframe')) {
    stop("Error: invalid test type. Please enter one of the following: 'expected', 'exclude', 'invalid_name', 'invalid_SIV', 'wrong_num_SIV', 'dataframe'.", call. = F)
  }

  ## For each model in HSImetadata, run the test 'iterations' number of times
   for (model in HSImetadata$model) {
     for (i in 1:iterations) {
       curves = data.frame(HSImodels[[model]]) # store curve breakpoints in data frame
       for (col in 1:length(curves)) {
         if (all(is.na(curves[,col])) == TRUE){
           count_NA = count_NA + 1 # identify and count columns with only NA values (i.e., not used in that model)
         } 
       }
       num_SIV = (ncol(curves) - count_NA)/2 # calculate the number of input variables
       input_vals = runif(num_SIV, 0, 1) # generate random test input from a uniform distribution between 0 and 1
       model.loc <- which(HSImetadata$model == model) # Find the location of the model in HSImetadata
       names(input_vals) = names(which(colSums(!is.na(HSImetadata[model.loc, 9:40])) > 0))
       ## If the test is for an invalid HSImodelname:
       if (test_type == 'invalid_name') {
         new_letters = str_flatten(sample(letters, 3)) # randomly choose three letters to tack on to the model name
         model = paste0(model, new_letters) # paste the model and letter together
       }
       ## If the test is for an invalid SIV value:
       else if (test_type == 'invalid_SIV') {
         # Inputs < 0 or > 1
         num_invalid = ifelse(num_SIV > 2, sample(1:(num_SIV - 2), 1), 1) # randomly select number of invalid inputs to be included  
         invalid_type = sample(1:3, 1)
         if (invalid_type == 1) {
           invalid = sample(
             runif(1000, -10.4, 10.4) |> (\(x) x[x < 0 | x > 1])(),
             size = num_invalid,
             replace = TRUE
           ) # randomly generate invalid values
         } else if (invalid_type == 2) {
           # Infinite inputs
           invalid = rep(Inf, num_invalid)
         } else if (invalid_type == 3) {
           # Character inputs
           invalid = sample(letters, num_invalid)
         }
         
         input_vals[sample(length(input_vals), num_invalid)] = invalid 
       } 
       ## If the test is for the wrong number of SIV values:
       else if (test_type == 'wrong_num_SIV') {
         num_extra = sample(1:20, 1) # randomly decide how many extra SIV values to include (between 1 and 20)
         extra_labels = paste0("SIV", 24:(24 + num_extra - 1))
         extra_SIV = runif(num_extra, min = 0, max = 1) # sample num_extra random values between zero and one
         input_vals = append(input_vals, extra_SIV)
         names(input_vals)[(num_SIV + 1):(num_SIV + num_extra)] = extra_labels
       } else if (test_type == 'dataframe') {
         input_vals = data.frame(t(input_vals))
       }
       HSI = suppressWarnings(try(HSIeqtn(model, input_vals, HSImetadata, exclude), silent = TRUE)) # calculate overall HSI if possible
       if (HSI >= 0 & HSI <= 1) {
         passes = passes + 1
         #print(model)
         #print(input_vals)
         #print(HSImetadata[HSImetadata$model == model, 'Eqtn'])
       }
       count_NA = 0 # Reset count_NA to zero
       num_SIV = 0 # Reset num_SIV to zero
       if (test_type == 'invalid_name') {
       model = substr(model, 1, nchar(model) - 3)
       }
     } 
   }
  # Pass rate
  if (test_type %in% c('expected', 'NA', 'dataframe', 'wrong_num_SIV')) {
    return(passes / (iterations * 349) * 100) # if results should be valid, return the pass rate
  } else {
    return(((iterations * 349) - passes) / (iterations * 349) * 100) # if results should be errors, return the fail rate
  }
}

HSIeqtn_tester('expected', 100)
HSIeqtn_tester('invalid_name', 100)
HSIeqtn_tester('invalid_SIV', 100)
HSIeqtn_tester('wrong_num_SIV', 100)
HSIeqtn_tester('dataframe', 100)

################################################################################

# SIcalc

SIcalc_tester = function(test_type, iterations) {
  passes = 0
  tests_run = 0
  count_NA = 0
  num_SIV = 0
  input = c()
  
  if (!test_type %in% c('manual', 'invalid', 'parameter_NA', 'HSImodels',
                        'mismatched_inputs', 'one_row_df', 'minmax_caps')) {
    stop(
      "Error: invalid test type. Please enter one of the following: 'manual', 'invalid', 'parameter_NA', 'HSImodels', 'mismatched_inputs', 'one_row_df', 'minmax_caps'.",
      call. = FALSE
    )
  }
  
  if (test_type == 'HSImodels') {
    for (model in HSImetadata$model) {
      for (i in 1:iterations) {
        input = c()
        count_NA = 0
        curves = HSImodels[[model]]
        
        for (col in 1:length(curves)) {
          if (all(is.na(curves[, col])) == TRUE) {
            count_NA = count_NA + 1
          }
          if (col %% 2 != 0) {
            n_rows = sum(!is.na(curves[, col]))
            input = c(input, as.character(sample(curves[1:n_rows, col], 1, FALSE)))
          }
        }
        
        SIs = try(SIcalc(curves, input), silent = TRUE)
        
        if (!inherits(SIs, "try-error") && all(SIs >= 0 & SIs <= 1 | is.na(SIs))) {
          passes = passes + 1
        }
      }
    }
    
  } else {
    for (test in 1:iterations) {
      input = c()
      model = c()
      num_vars = sample(1:10, 1)
      var_names = paste0('var_', 1:num_vars)
      siv_names = paste0('var_', 1:num_vars, '.SIV')
      var_is_cat = c()
      
      for (var in 1:num_vars) {
        num_SIV = sample(2:15, 1)
        cat = sample(0:1, 1)
        var_is_cat[var] = cat
        
        if (cat == 1) {
          categories = letters[1:num_SIV]
        } else {
          categories = as.numeric(sort(sample(0:1000, num_SIV)))
        }
        
        SIVs = runif(num_SIV, min = 0, max = 1)
        SIVs = c(SIVs, rep(NA, 15 - num_SIV))
        categories = c(categories, rep(NA, 15 - num_SIV))
        model = cbind(model, categories, SIVs)
      }
      
      colnames(model) = c(rbind(var_names, siv_names))
      model = data.frame(model, stringsAsFactors = FALSE)
      model = model[rowSums(is.na(model)) < ncol(model), ]
      
      model[] <- lapply(model, function(x) {
        num <- suppressWarnings(as.numeric(x))
        if (all(is.na(x) | !is.na(num))) {
          num
        } else {
          x
        }
      })
      
      if (test_type == 'parameter_NA') {
        odd_indices = seq(1, length(model), by = 2)
        na_col = sample(odd_indices, 1)
        model[, na_col:(na_col + 1)] = NA
      }
      
      for (col in 1:length(model)) {
        if (col %% 2 != 0) {
          if (all(is.na(model[, col]))) {
            input = append(input, NA)
          } else if (is.numeric(model[, col])) {
            min_val = min(model[, col], na.rm = TRUE)
            max_val = max(model[, col], na.rm = TRUE)
            input = append(input, runif(1, min = min_val, max = max_val))
          } else {
            min_val = min(match(model[, col], letters), na.rm = TRUE)
            max_val = max(match(model[, col], letters), na.rm = TRUE)
            input = append(input, letters[sample(min_val:max_val, 1)])
          }
        }
      }
      
      if (test_type == 'invalid') {
        inv_input = sample(1:length(input), 1)
        invalid_type = sample(1:3, 1)
        
        if (invalid_type == 1) {
          if (input[inv_input] %in% letters) {
            input[inv_input] = letters[sample(16:26, 1)]
          } else {
            input[inv_input] = letters[sample(1:26, 1)]
          }
        } else if (invalid_type == 2) {
          if (input[inv_input] %in% letters) {
            input[inv_input] = sample(c(1500, -1500), 1)
          } else {
            input[inv_input] = letters[sample(1:26, 1)]
          }
        } else {
          input[inv_input] = Inf
        }
      }
      
      if (test_type == 'minmax_caps') {
        continuous_vars <- which(var_is_cat == 0)
        
        if (length(continuous_vars) == 0) {
          next
        }
        
        chosen_var <- sample(continuous_vars, 1)
        cap_col <- 2 * chosen_var - 1
        input_index <- chosen_var
        
        x_num <- suppressWarnings(as.numeric(as.character(model[, cap_col])))
        x_num <- x_num[!is.na(x_num)]
        
        if (length(x_num) == 0) {
          next
        }
        
        tests_run = tests_run + 1
        
        min_val <- min(x_num)
        max_val <- max(x_num)
        
        if (sample(c(TRUE, FALSE), 1)) {
          input[input_index] <- min_val - runif(1, 1, 100)
        } else {
          input[input_index] <- max_val + runif(1, 1, 100)
        }
      }
      
      if (test_type == 'mismatched_inputs') {
        num_extra = sample(1:20, 1)
        input = append(input, sample(0:1000, num_extra))
      }
      
      if (test_type == 'one_row_df') {
        input <- as.data.frame(as.list(input), stringsAsFactors = FALSE)
      }
      
      SIs = suppressWarnings(try(SIcalc(model, input), silent = TRUE))
      
      if (!inherits(SIs, "try-error") && all(SIs >= 0 & SIs <= 1 | is.na(SIs))) {
        passes = passes + 1
      }
    }
  }
  
  if (test_type == 'HSImodels') {
    return(passes / (iterations * 349) * 100)
  } else if (test_type %in% c('manual', 'parameter_NA', 'one_row_df')) {
    return(passes / iterations * 100)
  } else if (test_type == 'minmax_caps') {
    return(passes / tests_run * 100)
  } else {
    return((iterations - passes) / iterations * 100)
  }
}

SIcalc_tester('HSImodels', 100)
SIcalc_tester('manual', 9999)
SIcalc_tester('parameter_NA', 9999)
SIcalc_tester('invalid', 9999) #now tests only wrong categorical inputs
SIcalc_tester('mismatched_inputs', 9999) 
SIcalc_tester('one_row_df', 9999)
SIcalc_tester('minmax_caps', 9999)
################################################################################

# HUcalc

HUcalc_tester = function(test_type, iterations) {
  passes = 0 # set initial passes to zero
  
  # Make sure test_types are correct
  if (!test_type %in% c('expected', 'invalid', 'NA', 'invalid_area', 'dataframe')) {
    stop("Error: invalid test type. Please enter one of the following: 'NA', 'expected', 'invalid', 'invalid_area', 'dataframe'.", call. = F)
  }
  
  functions = c('HSIarimean', 'HSIgeomean', 'HSImin', 'HSIwarimean', 'mean', 'max', 'min') # functions to use
  
  for (i in 1:iterations) {
    inputs = c() # reset inputs
    num_SIV = sample(2:20, 1) # randomly select a number of inputs
    inputs = runif(num_SIV, min = 0, max = 1) # randomly generate inputs
    area = sample(1:500, 1) # randomly generate area
    use = sample(1:7, 1)
    if (test_type == 'NA') {
      num_NA = sample(1:(num_SIV - 1), 1)
      inputs[sample(length(inputs), num_NA)] = NA # randomly replace some values with NA
    } else if (test_type == 'invalid') {
      num_invalid = sample(1:num_SIV, 1) # randomly decide number of invalid
      type_invalid = sample(1:3, 1) # randomly determine the type of invalid input
      if (type_invalid == 1) {
        invalid = sample(
          runif(1000, -10.4, 10.4) |> (\(x) x[x < 0 | x > 1])(),
          size = num_invalid,
          replace = TRUE
        ) # randomly generate invalid values
      } else if (type_invalid == 2) {
        invalid = sample(letters, num_invalid)
      } else {
        invalid = rep(Inf, num_invalid)
      }
    
      inputs[sample(length(inputs), num_invalid)] = invalid 
    } else if (test_type == 'invalid_area') {
        invalid_type = sample(1:3, 1)
        if (invalid_type == 1) {
          area = sample(-1000:-0.01, 1)
        } else if (invalid_type == 2) {
          area = sample(letters, 1)
        } else {
          area = Inf
        }
    } else if (test_type == 'dataframe') {
      inputs = data.frame(inputs)
    }
    if (functions[use] == 'HSIwarimean') {
      weights_raw = runif(num_SIV) # randomly generate weights
      
      weights <- weights_raw
      
      weights[is.na(inputs)] <- NA
      
      weights <- weights / sum(weights, na.rm = TRUE)
      
      weights[max(which(!is.na(weights)))] <-
        weights[max(which(!is.na(weights)))] + (1 - sum(weights, na.rm = TRUE))
      
      hu = suppressWarnings(try(HUcalc(inputs, area, eval(parse(text = functions[use])), weights), silent = T))
    } else {
      hu = suppressWarnings(try(HUcalc(inputs, area, eval(parse(text = functions[use]))), silent = T))
    }
    if (!inherits(hu, 'try-error') & is.data.frame(hu)) {
      passes = passes + 1
    } 
  }
  # Pass rate
  if (test_type %in% c('expected', 'NA', 'dataframe')) {
    return(passes / iterations * 100) # if results should be valid, return the pass rate
  } else {
    return((iterations - passes) / iterations * 100) # if results should be errors, return the fail rate
  }
}

HUcalc_tester('expected', 9999)
HUcalc_tester('NA', 9999)
HUcalc_tester('invalid', 9999)
HUcalc_tester('invalid_area', 9999)
HUcalc_tester('dataframe', 9999)

################################################################################

# HSIplotter

HSIplotter_tester = function(test_type, iterations = 349) {
  # Make sure test_type is correct
  if (!test_type %in% c('blue-book', 'manual', 'invalid')) {
    stop("Invalid test type. Please enter one of the following: 'blue-book', 'manual', 'invalid'.", call. = F)
  }
  
  # If testing blue book models:
  if (test_type == 'blue-book') {
    if(!dir.exists('Bluebook_tests')) {
      dir.create("Bluebook_tests")
    }
    for (i in 1:nrow(HSImetadata)) {
      mod=HSImodels[i]
      name=names(mod[1])
      model <- get(name,mod)
      try(HSIplotter(model, paste(getwd(), '/Bluebook_tests', '/', name, sep = "",".png")))
    }
  }
  # If testing manual inputs
  else if (test_type == 'manual') {
    if(!dir.exists('Manual_tests')) {
      dir.create("Manual_tests")
    }
    for (test in 1:iterations) {
      model = c() # reset model
      num_vars = sample(1:10, 1) # randomly set the number of variables to a value between one and ten
      var_names = paste0('var_', 1:num_vars)
      siv_names = paste0('var_', 1:num_vars, '.SIV')
      for (var in 1:num_vars) {
        num_SIV = sample(2:15, 1) # randomly set the number of SIVs between 1 and 15
        cat = sample(0:1, 1) # randomly choose if variable is categorical (1) or not (0)
      
        if (cat == 1) {
        categories = letters[1:num_SIV] # set categories to equal letters
        } else {
          categories = as.numeric(sort(sample(0:1000, num_SIV))) # if continuous, randomly select values between 0 and 1000
        }
        ## Randomly generate SIVs
        SIVs = runif(num_SIV, min = 0, max = 1)
        ## Add NAs to equal 15
        SIVs = c(SIVs, rep(NA, 15 - num_SIV)) # add NAs if needed
        categories = c(categories, rep(NA, 15 - num_SIV)) # add NAs if needed
        model = cbind(model, categories, SIVs) # bind data
      }
        colnames(model) = c(rbind(var_names, siv_names)) # name columns
        model = data.frame(model) # convert to dataframe
        model = model[rowSums(is.na(model)) < ncol(model), ] # Remove rows that are all NA
        model[] <- lapply(model, function(x) {
          num <- suppressWarnings(as.numeric(x))
          
          # convert only if every non-NA value successfully became numeric
          if (all(is.na(x) | !is.na(num))) {
            num
          } else {
            x
          }
        })
        try(HSIplotter(model, paste(getwd(), '/Manual_tests/model', test, sep = "",".png"))) # try graphing
        model = c() # reset model
    }
    dev.off()
  } else if (test_type == 'invalid') {
      if(!dir.exists('Invalid_tests')) {
        dir.create("Invalid_tests")
      }
      for (test in 1:iterations) {
        model = c() # reset model
        num_vars = sample(1:10, 1) # randomly set the number of variables to a value between one and ten
        var_names = paste0('var_', 1:num_vars)
        siv_names = paste0('var_', 1:num_vars, '.SIV')
        for (var in 1:num_vars) {
          num_SIV = sample(2:15, 1) # randomly set the number of SIVs between 1 and 15
          cat = sample(0:1, 1) # randomly choose if variable is categorical (1) or not (0)
        
          if (cat == 1) {
            categories = letters[1:num_SIV] # set categories to equal letters
          } else {
            categories = as.numeric(sort(sample(0:1000, num_SIV))) # if continuous, randomly select values between 0 and 1000
          }
          ## Randomly generate SIVs
          SIVs = runif(num_SIV, min = 0, max = 1)
          ## Add NAs to equal 15
          SIVs = c(SIVs, rep(NA, 15 - num_SIV)) # add NAs if needed
          categories = c(categories, rep(NA, 15 - num_SIV)) # add NAs if needed
          model = cbind(model, categories, SIVs) # bind data
        }
        colnames(model) = c(rbind(var_names, siv_names)) # name columns
        model = data.frame(model) # convert to dataframe
        model = model[rowSums(is.na(model)) < ncol(model), ] # Remove rows that are all NA
        model[] <- lapply(model, function(x) {
          num <- suppressWarnings(as.numeric(x))
        
          # convert only if every non-NA value successfully became numeric
          if (all(is.na(x) | !is.na(num))) {
            num
          } else {
            x
          }
        })
        invalid_loc = sample(1:length(model), 1) # randomly select a column to mess up
        invalid_type = sample(1:2, 1) # randomly select the type of invalid input
        
        if (invalid_type == 1) {
          if (invalid_loc %% 2 == 0) {
            model[1, invalid_loc] = sample(
              runif(1000, -10.4, 10.4) |> (\(x) x[x < 0 | x > 1])(),
              size = 1,
              replace = TRUE
            ) # randomly generate invalid values
          } else {
            invalid_loc = invalid_loc + 1
            model[1, invalid_loc] = sample(
              runif(1000, -10.4, 10.4) |> (\(x) x[x < 0 | x > 1])(),
              size = 1,
              replace = TRUE
            ) # randomly generate invalid values
          }
        } else {
          model[1, invalid_loc] = Inf
        }
        
        try(HSIplotter(model, paste(getwd(), '/Invalid_tests/model', test, sep = "",".png")), silent = T) # try graphing
        model = c() # reset model
      }
      dev.off()
  }
}

HSIplotter_tester('blue-book')
HSIplotter_tester('manual', 349)
HSIplotter_tester('invalid', 349)

################################################################################

# annualizer

annualizer_test = function(test_type = 'expected', iterations) {
  passes = 0 # Set number of passes to zero
  if (!test_type %in% c('NA', 'expected', 'wrong_time', 'wrong_benefit', 'wrong_order', 'duplicate_time', 'invalid_time', 'invalid_benefit', 'dataframe')) {
    stop("Error: invalid test type. Please enter one of the following: 'NA', 'expected', 'wrong_time', 'wrong_benefit', 'wrong_order, 'duplicate_time', 'invalid_time', 'invalid_benefit', 'dataframe'.", call. = F)
  }
  
  for (i in 1:iterations) {
    num_inputs = sample(2:100, 1) # randomly generate number of inputs between 2 and 100
    timevec = sort(sample(0:100, size = num_inputs)) # Randomly select time values with size num_inputs 
    benefits = sample(0:100, size = num_inputs) # Randomly generate benefit values with size num_inputs
    num_NA = 0 # specify no NA values
    
  if (test_type == 'NA') {
    num_NA = sample(1:(num_inputs - 1), 1) # randomly select number of NA values to be included
    timevec[sample(length(timevec), num_NA)] = NA # randomly replace some values with NA
    na_locations <- is.na(timevec) # store where NA values occur
    benefits[na_locations] = NA # replace same locations in benefits vector with NA
  }
  
  else if (test_type == 'wrong_time') {
    num_extra = sample(1:20, 1) # decide how many extra values to provide as inputs
    timevec = append(timevec, sample(1:500, num_extra)) # append to timevec vector
  }
  
  else if (test_type == 'wrong_benefit') {
    num_extra = sample(1:20, 1) # decide how many extra values to provide as inputs
    benefits = append(benefits, sample(1:500, num_extra)) # append to benefits vector
  } 
  
  else if (test_type == 'wrong_order') {
    # Shuffle timevec until it is out of order
    while (identical(timevec, sort(timevec))) {
      timevec = sample(timevec)
    }
  }
    
  else if (test_type == 'duplicate_time') {
    dup_position = sample(1:length(timevec), 1) # select time value to duplicate
    rep_position = sample(1:length(timevec), 1) # select value to replace with duplicate
    # Make sure the two positions are not the same
    while (dup_position == rep_position) {
      rep_position = sample(1:length(timevec), 1) # select value to replace with duplicate
    }
    timevec[rep_position] = timevec[dup_position]
  }
    
  else if (test_type == 'invalid_time') {
    invalid_type = sample(1:3, 1) # randomly select what type of invalid input
    if (invalid_type == 1) {
      timevec[1] = sample(-1000:-0.01, 1)
    } else if (invalid_type == 2) {
      timevec[1] = sample(letters, 1)
    } else {
      timevec[1] = Inf
    }
  }
    
  else if (test_type == 'invalid_benefit') {
    invalid_type = sample(1:3, 1) # randomly select what type of invalid input
    if (invalid_type == 1) {
      benefits[1] = sample(-1000:-0.01, 1)
    } else if (invalid_type == 2) {
      benefits[1] = sample(letters, 1)
    } else {
      benefits[1] = Inf
    }
  }
    
  else if (test_type == 'dataframe') {
    timevec = data.frame(timevec)
    benefits = data.frame(benefits)
  }
    
  
  result = suppressWarnings(try(annualizer(timevec, benefits), silent = T)) # Calculate time-averaged value
  
  ## If the result is numeric and is not NA or NA_real, pass
  if (is.numeric(result) & !is.na(result))
    passes = passes + 1
  }
  # Pass rate
  if (test_type %in% c('expected', 'wrong_order', 'dataframe')) {
    return(passes / iterations * 100) # if results should be valid, return the pass rate
  } else {
    return((iterations - passes) / iterations * 100) # if results should be errors, return the fail rate
  }
}

annualizer_test(test_type = 'expected', iterations = 9999)
annualizer_test(test_type = 'NA', iterations = 9999)
annualizer_test(test_type = 'wrong_time', iterations = 9999)
annualizer_test(test_type = 'wrong_benefit', iterations = 9999)
annualizer_test(test_type = 'wrong_order', iterations = 9999)
annualizer_test(test_type = 'duplicate_time', iterations = 9999)
annualizer_test(test_type = 'invalid_time', iterations = 9999)
annualizer_test(test_type = 'invalid_benefit', iterations = 9999)
annualizer_test(test_type = 'dataframe', iterations = 9999)

# Test specific inputs to ensure math calculations are correct
timevec=c(0, 10) 
benefits = c(2, 10)

annualizer(timevec, benefits)

################################################################################

# CEfinder

CEfinder_tester = function(test_type = 'expected', iterations) {
  passes = 0 # Set number of passes to zero
  ## Make sure the correct tests are being run
  if (!test_type %in% c('NA', 'expected', 'wrong_cost', 'wrong_benefit', 'invalid_cost', 'invalid_benefit', 'dataframe')) {
    stop("Error: invalid test type. Please enter one of the following: 'NA', 'expected', 'wrong_cost', 'wrong_benefit', 'invalid_cost', 'invalid_benefit', 'dataframe'.", call. = F)
  }
  for (i in 1:iterations) {
    num_inputs = sample(2:20, 1) # randomly select a number of samples
    benefit = sample(0:100, num_inputs) # randomly sample num_inputs benefits
    cost = sample(0:500, num_inputs) # randomly sample num_inputs costs
    
    if (test_type == 'NA') {
      num_NA = sample(1:(num_inputs - 1), 1) # randomly select number of NA values to be included
      benefit[sample(length(benefit), num_NA)] = NA # randomly replace some values with NA
      na_locations <- is.na(benefit) # store where NA values occur
      cost[na_locations] = NA # replace same locations in benefits vector with NA
    } else if (test_type == 'wrong_benefit') {
      num_extra = sample(1:20, 1) # decide how many extra values to provide as inputs
      benefit = append(benefit, sample(1:500, num_extra)) # append to benefits vector
    } else if (test_type == 'wrong_cost') {
      num_extra = sample(1:20, 1) # decide how many extra values to provide as inputs
      cost = append(cost, sample(1:500, num_extra)) # append to benefits vector
    }
    
    else if (test_type == 'invalid_cost') {
      invalid_type = sample(1:3, 1) # randomly select what type of invalid input
      if (invalid_type == 1) {
        cost[1] = sample(-1000:-0.01, 1)
      } else if (invalid_type == 2) {
        cost[1] = sample(letters, 1)
      } else {
        cost[1] = Inf
      }
    }
    
    else if (test_type == 'invalid_benefit') {
      invalid_type = sample(1:3, 1) # randomly select what type of invalid input
      if (invalid_type == 1) {
        benefit[1] = sample(-1000:-0.01, 1)
      } else if (invalid_type == 2) {
        benefit[1] = sample(letters, 1)
      } else {
        benefit[1] = Inf
      }
    }
    
    else if (test_type == 'dataframe') {
      benefit = data.frame(benefit)
      cost = data.frame(cost)
    }
    
    CE = try(CEfinder(benefit, cost), silent = T) # calculate cost effectiveness of each plan
    
    ## If all outputs are binary, pass
    if (all(unique(CE) %in% c(0, 1))) {
      passes = passes + 1
    }
  }
  # Pass rate
  if (test_type %in% c('expected', 'dataframe')) {
    return(passes / iterations * 100) # if results should be valid, return the pass rate
  } else {
    return((iterations - passes) / iterations * 100) # if results should be errors, return the fail rate
  }
}


CEfinder_tester('expected', 9999)
CEfinder_tester('NA', 9999)
CEfinder_tester('wrong_benefit', 9999)
CEfinder_tester('wrong_cost', 9999)
CEfinder_tester('invalid_cost', 9999)
CEfinder_tester('invalid_benefit', 9999)
CEfinder_tester('dataframe', 9999)

################################################################################


# BBfinder

BBfinder_tester = function(test_type = 'expected', iterations) {
  passes = 0 # Set number of passes to zero
  ## Make sure the correct tests are being run
  if (!test_type %in% c('NA', 'expected', 'wrong_cost', 'wrong_benefit', 'wrong_plans', 'invalid_benefit', 'invalid_cost', 'invalid_CE', 'duplicate_benefit', 'dataframe')) {
    stop("Error: invalid test type. Please enter one of the following: 'NA', 'expected', 'wrong_cost', 'wrong_benefit', 'wrong_plans', 'invalid', 'invalid_benefit', 'invalid_cost', 'invalid_CE', 'duplicate_benefit, 'dataframe'.", call. = F)
  }
    for (i in 1:iterations) {
      num_inputs = sample(2:20, 1) # randomly select a number of samples
      benefit = sample(0:100, num_inputs) # randomly sample num_inputs benefits
      cost = sample(0:500, num_inputs) # randomly sample num_inputs costs
      CE = CEfinder(benefit, cost)
      
      while (sum(CE, na.rm = T) <= 1) {
        num_inputs = sample(2:20, 1) # randomly select a number of samples
        benefit = sample(0:100, num_inputs) # randomly sample num_inputs benefits
        cost = sample(0:500, num_inputs) # randomly sample num_inputs costs
        CE = CEfinder(benefit, cost)
      }
       ## If NAs are included, add NAS
      if (test_type == 'NA') {
        num_NA = sample(1:(num_inputs - 1), 1) # randomly select number of NA values to be included
        benefit[sample(length(benefit), num_NA)] = NA # randomly replace some values with NA
        na_locations <- is.na(benefit) # store where NA values occur
        cost[na_locations] = NA # replace same locations in cost vector with NA
        CE[na_locations] = NA # replace same locations in CE vector with NA
        
        while (sum(CE, na.rm = T) <= 1) {
          num_inputs = sample(2:20, 1) # randomly select a number of samples
          benefit = sample(0:100, num_inputs) # randomly sample num_inputs benefits
          cost = sample(0:500, num_inputs) # randomly sample num_inputs costs
          CE = CEfinder(benefit, cost)
          
          num_NA = sample(1:(num_inputs - 1), 1) # randomly select number of NA values to be included
          benefit[sample(length(benefit), num_NA)] = NA # randomly replace some values with NA
          na_locations <- is.na(benefit) # store where NA values occur
          cost[na_locations] = NA # replace same locations in cost vector with NA
          CE[na_locations] = NA # replace same locations in CE vector with NA
        }
      } else if (test_type == 'wrong_benefit') {
        num_extra = sample(1:20, 1) # decide how many extra values to provide as inputs
        benefit = append(benefit, sample(1:500, num_extra)) # append to benefit vector
      } else if (test_type == 'wrong_cost') {
        num_extra = sample(1:20, 1) # decide how many extra values to provide as inputs
        cost = append(cost, sample(1:500, num_extra)) # append to benefits vector
      } else if (test_type == 'wrong_plans') {
        num_extra = sample(1:20, 1) # decide how many extra values to provide as inputs
        CE = append(CE, sample(0:1, num_extra, replace = T)) # append to benefits vector
      } else if (test_type == 'invalid_cost') {
          invalid_type = sample(1:3, 1) # randomly select what type of invalid input
          if (invalid_type == 1) {
            cost[1] = sample(-1000:-0.01, 1)
          } else if (invalid_type == 2) {
            cost[1] = sample(letters, 1)
          } else {
            cost[1] = Inf
          }
      } else if (test_type == 'invalid_benefit') {
          invalid_type = sample(1:3, 1) # randomly select what type of invalid input
          if (invalid_type == 1) {
            benefit[1] = sample(-1000:-0.01, 1)
          } else if (invalid_type == 2) {
            benefit[1] = sample(letters, 1)
          } else {
            benefit[1] = Inf
          }
      } else if (test_type == 'invalid_CE') {
          invalid_type = sample(1:4, 1) # randomly select what type of invalid input
          if (invalid_type == 1) {
            CE[1] = sample(-1000:-0.01, 1)
          } else if (invalid_type == 2) {
            CE[1] = sample(letters, 1)
          } else if (invalid_type == 3) {
            CE[1] = Inf
          } else {
            CE[1] = sample(2:200, 1)
          } 
      } else if (test_type == 'duplicate_benefit') {
          dup_position = sample(1:length(benefit), 1) # select time value to duplicate
          rep_position = sample(1:length(benefit), 1) # select value to replace with duplicate
          # Make sure the two positions are not the same
          while (dup_position == rep_position) {
            rep_position = sample(1:length(benefit), 1) # select value to replace with duplicate
          }
          benefit[rep_position] = benefit[dup_position]
      } else if (test_type == 'dataframe') {
        benefit = data.frame(benefit)
        cost = data.frame(cost)
        CE = data.frame(CE)
      }
    
      
      
      result = suppressWarnings(try(BBfinder(benefit, cost, CE), silent = TRUE))
      
      if (!inherits(result, 'try-error') & is.list(result)) {
        passes = passes + 1
      }
    }
    # Pass rate
    if (test_type %in% c('expected', 'duplicate_benefit', 'dataframe')) {
      return(passes / iterations * 100) # if results should be valid, return the pass rate
    } else {
      return((iterations - passes) / iterations * 100) # if results should be errors, return the fail rate
    }
}

BBfinder_tester('expected', 9999)
BBfinder_tester('NA', 9999)
BBfinder_tester('wrong_benefit', 9999)
BBfinder_tester('wrong_cost', 9999)
BBfinder_tester('wrong_plans', 9999)
BBfinder_tester('invalid_benefit', 9999)
BBfinder_tester('invalid_cost', 9999)
BBfinder_tester('invalid_CE', 9999)
BBfinder_tester('duplicate_benefit', 9999)
BBfinder_tester('dataframe', 9999)

################################################################################

## CEICAplotter

CEICAplotter_tester = function(test_type = 'expected', iterations) {
  ## Make sure the correct tests are being run
  if (!test_type %in% c('expected', 'dataframe')) {
    stop("Error: invalid test type. Please enter one of the following: 'expected', 'dataframe.", call. = F)
  }
  
  if(!dir.exists('CEICA_tests')) {
    dir.create("CEICA_tests")
  }
  
  for (i in 1:iterations) {
    num_inputs = sample(1:20, 1) # randomly select a number of samples
    benefit = sample(0:100, num_inputs) # randomly sample num_inputs benefits
    cost = sample(0:500, num_inputs) # randomly sample num_inputs costs
    CE <- CEfinder(benefit, cost)
    
    while (sum(CE, na.rm = T) <= 1) {
      num_inputs = sample(2:20, 1) # randomly select a number of samples
      benefit = sample(0:100, num_inputs) # randomly sample num_inputs benefits
      cost = sample(0:500, num_inputs) # randomly sample num_inputs costs
      CE = CEfinder(benefit, cost)
    }
    
    BB <- BBfinder(benefit, cost, CE)[[1]][,4]
    altnames = paste("Alt", seq(1:num_inputs), sep = "")
    
    if (test_type == 'dataframe') {
      altnames = data.frame(altnames)
      benefit = data.frame(benefit)
      cost = data.frame(cost)
      CE = data.frame(CE)
      BB = data.frame(BB)
    }
    try(CEICAplotter(altnames, benefit, cost, CE, BB, paste(getwd(), '/CEICA_tests/test', i, sep = "",".jpeg")))
  }
}

CEICAplotter_tester('expected', 300)
CEICAplotter_tester('dataframe', 300)
