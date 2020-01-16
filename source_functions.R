## S2MET Functions
## 
## A script with useful functions used in the prediction analysis of the S2MET
## 



# ## A generic prediction function that takes training and test data and returns
# ## PGVs and accuracy
# gblup <- function(formula, random, K, train, test, fun = c("rrblup", "sommer"), fit.env = TRUE, bootreps = NULL, add.mu = FALSE) {
#   
#   if (missing(formula)) formula <- value ~ 1 + environment
#   if (missing(random)) random <- ~ line_name
#   
#   ## separate random and fixed terms
#   fixed_terms <- attr(terms(formula), "term.labels")
#   random_terms <- attr(terms(random), "term.labels")
# 
#   ## Combine fixed and random to one formula
#   formula1 <- as.formula(paste(paste(as.character(formula)[c(2,1,3)], collapse = " "), as.character(random)[-1], sep = " + "))
#   
#   # Create a model.frame
#   mf <- model.frame(formula1, weights = std_error, data = train)
#   
#   
#   ## If the number of random terms is > 1, K must be a list
#   if (length(random_terms) == 1) {
#     stopifnot(levels(mf[[random_terms]]) %in% colnames(K))
#     
#     K <- list(K)
#     names(K) <- random_terms
#     
#   } else {
#     if (!is.list(K)) stop("If the number of random terms is > 1, K must be a list of relationship matrices.")
#     
#     ## Test names
#     name_test <- map2_lgl(.x = random_terms, .y = names(K1), ~.x == .y)
#     if (!all(name_test)) stop("If K is a list, the names of the list must match the random terms")
#     
#   }
#    
#   
#   fun <- match.arg(fun)
#   
#   
#   
#   # Vectors and matrices
#   y <- model.response(mf)
#   
#   
#   if (nlevels(mf[[fixed_terms]]) <= 1 | !fit.env) {
#     X <- model.matrix(~ 1, droplevels(mf))
#   
#     } else {
#     X <- model.matrix(formula, droplevels(mf))
#     
#   }
#   
#   ## Random effects
#   Z_list <- list()
#   
#   for (term in random_terms) {
#     Zi <- model.matrix(as.formula(paste0("~ -1 + ", term)), mf)
#     colnames(Zi) <- colnames(K[[term]])
#     Z_list[[term]] <- Zi
#     
#   }
#   
#   
#   # Split on function
#   if (fun == "rrblup") {
#     # Can't run with > 1 random term
#     stopifnot(length(random_terms) == 1)
#     
#     fit <- mixed.solve(y = y, Z = Z_list[[1]], K = K[[1]], X = X)
#     
#     # Extract PGVs
#     pgv <- fit$u %>% 
#       data.frame(line_name = names(.), pred_value = ., row.names = NULL, stringsAsFactors = FALSE)
#     
#     beta <- fit$beta[1]
#     
#     
#   } else if (fun == "sommer") {
#     
#     # R <- solve(diag(mf$`(weights)`^2))
#     # fit <- sommer::mmer(Y = y, X = X, Z = list(g = list(Z = Z, K = K)), R = list(res = R), silent = TRUE)
#     
#     random_list <- map2(Z_list, K, ~list(Z = .x, K = .y))
#     fit <- sommer::mmer(Y = y, X = X, Z = random_list, silent = TRUE)
#   
#     # Extract PGVs
#     pgv <- fit$u.hat
#     pgv <- data.frame(line_name = row.names(pgv), pred_value = pgv[,1], row.names = NULL, stringsAsFactors = FALSE)
#     
#     beta <- fit$beta.hat[1]
#     
#   }
#   
#   if (add.mu) pgv$pred_value <- pgv$pred_value + beta
#   
#   # If test is missing, just return the predictions
#   if (missing(test)) {
#     
#     comb <- pgv
#     acc <- boot <- NA
#     
#   } else {
#   
#     # Combine the PGVs with the phenotypic observations and calculate accuracy
#     comb <- left_join(test, pgv, by = "line_name")
#     acc <- cor(comb$value, comb$pred_value)
#     
#     # Bootstrap if replicates are provided
#     if (!is.null(bootreps)) {
#       boot <- bootstrap(x = comb$value, y = comb$pred_value, fun = "cor", boot.reps = bootreps)
#     } else {
#       boot <- NA
#     }
#     
#   }
#   
#   # Return a list
#   list(accuracy = acc, pgv = comb, boot = as_data_frame(boot))
# }




## A model fitting function
## 
## All models are random
## 
## Model1: y = G
## Model2: y = G + E
## Model3: y = G + E + GE
## 
## Note: for single traits only
## 

predict_gv <- function(train, test, model = c("model1", "model2", "model3"), relMat.list, covariate.list = list(NULL),
                       object, add.fixed = FALSE, verbose = FALSE) {

  ## Depending on the model, make sure correct relmats are present
  model <- match.arg(model)

  ## Assign required relmats
  req_relmat <- switch(
    model,
    model1 = c("G"),
    model2 = c("G", "E"),
    model2a = c("G"),
    model2b = c("G"),
    model3 = c("G", "E", "GE"),
    model3a = c("G", "GE"),
    model3b = c("G", "GE")
  )

  ## Assign ranefs
  ranefs <- switch(
    model,
    model1 = c("line_name"),
    model2 = c("line_name", "environment"),
    model2a = c("line_name"),
    model2b = c("line_name"),
    model3 = c("line_name", "environment", "line_name:environment"),
    model3a = c("line_name", "line_name:environment"),
    model3b = c("line_name", "line_name:environment")
  )

  ## Assign fixefs
  fixef <- switch(
    model,
    model1 = NULL,
    model2 = NULL,
    model2a = c("environment"),
    model2b = NULL,
    model3 = NULL,
    model3a = c("environment"),
    model3b = NULL
  )

  ## Assign covariates
  covariates <- switch(
    model,
    model1 = NULL,
    model2 = NULL,
    model2a = NULL,
    model2b = "E",
    model3 = NULL,
    model3a = NULL,
    model3b = "E"
  )


  ## Error handling

  ## Check for presence of relmats
  stopifnot(is.list(relMat.list))
  stopifnot(all(req_relmat %in% names(relMat.list)))

  ## Check for presence of covariates
  stopifnot(is.list(covariate.list))
  stopifnot(all(covariates %in% names(covariate.list)))


  ## Extract covariates
  covariates <- unlist(covariate.list[covariates])


  ## Create fixed/random formula
  fixed_formula <- reformulate(termlabels = c("1", c(fixef, covariates)), response = "value")


  # If object is missing, do not include contraints on the random effects
  if (missing(object)) {
    random_formula <- reformulate(termlabels = mapply(req_relmat, ranefs, FUN = function(.x, .y) {
      paste0("vs(", .y, ", Gu = relMat.list$", .x, ")") }))


    # Residual formula
    rcov_form <- ~ units

    # If not missing, extract variance components and impose restrictions
  } else {
    # Error handling
    # 'object' must have element called sigma
    stopifnot("sigma" %in% names(object))

    # Extract varcomps
    object_varcomps <- object$sigma
    names(object_varcomps) <- gsub(pattern = "u:", replacement = "", x = names(object_varcomps))
    object_varcomps <- object_varcomps[c(ranefs, "units")]

    # Create formula
    random_formula <- reformulate(termlabels = mapply(req_relmat, ranefs, FUN = function(.x, .y) {
      paste0("vs(", .y, ", Gu = relMat.list$", .x, ", Gt = object_varcomps[['", .y, "']], Gtc = fixm(1))") }) )

    # Residual formula
    rcov_form <- ~ vs(units, Gt = object_varcomps$units, Gtc = fixm(1))
  }


  ## Create a model frame
  full_formula <- add_predictors(fixed_formula, reformulate(ranefs))
  mf <- model.frame(formula = full_formula, train, drop.unused.levels = FALSE)

  # X matrix
  X <- model.matrix(fixed_formula, train)

  # List of Z matrices
  Zlist <- lapply(ranefs, function(term) model.matrix(reformulate(c(-1, term)), mf))
  names(Zlist) <- ranefs

  ## Fit the model
  fit <- mmer(fixed = fixed_formula, random = random_formula, rcov = rcov_form, data = mf, date.warning = FALSE, verbose = verbose)

  # If interaction present, split and take those ranefs
  if (str_detect(ranefs, ":")) {
    ranefs_no_int <- str_split(str_subset(ranefs, ":"), pattern = ":")[[1]]
    
  } else {
    ranefs_no_int <- str_subset(ranefs, ":", negate = TRUE)
    
  }

  ## Figure out what to return
  # If test is missing, return the object
  # Create a list to return
  if (missing(test)) {
    to_return <- list(object = fit)

  } else {

    ## Get the random effects
    U <- lapply(X = randef(fit), FUN = "[[", "value")
    # Rename
    names(U) <- str_remove(string = names(U), pattern = "u:")

    ## Use incidence matrices to predict random effs
    Zu_list <- mapply(Zlist, U, FUN = function(Z, u) {
      Zu <- Z %*% u
      `row.names<-`(Zu, names(u)[apply(X = Z, MARGIN = 1, FUN = function(row) which(row == 1))])
    }, SIMPLIFY = FALSE)

    ## Find unique values and convert to df
    Zu_df <- lapply(Zu_list, function(x) as.data.frame(unique(x))) %>%
      map(~rownames_to_column(., "term")) %>%
      list(., names(.), seq_along(.)) %>%
      pmap(., ~`names<-`(..1, c(..2, paste0("pred", ..3)))) %>%
      modify_if(.x = ., .p = ~str_detect(names(.)[1], ":"), ~separate(data = .x, col = 1, into = c("line_name", "environment"), sep = ":"))
    
    # If any dfs in the list have > 2 columns, sort and merge. 
    # Otherwise combine and full join
    if (any(map_lgl(Zu_df, ~ncol(.) > 2))) {
      Zu_join <- reduce(.x = Zu_df[order(map_dbl(Zu_df, ncol), decreasing = TRUE)], .f = left_join)
      
    } else {
      combn <- Zu_df[ranefs_no_int] %>% 
        map(select, 1) %>% 
        reduce(crossing)
      Zu_join <- reduce(c(list(combn), Zu_df), left_join)
      
    }
    
    ## Sum for predictions
    rand_pred <- Zu_join %>%
      mutate(prediction = rowMeans(select(., contains("pred")))) %>%
      select(ranefs_no_int, prediction)
    
    
    
    #### Fixed effects ####
    beta <- as.matrix(column_to_rownames(fit$Beta[,-1], "Effect"))
    Xb <- X %*% beta
    
    
    
    

    ## Extract predictions
    fit_pred <- randef(fit) %>%
      map("value") %>%
      map(~tibble(term = names(.x), pred = .x), .vars = vars(term)) %>%
      imap(~`names<-`(.x, c(.y, "pred"))) %>%
      unname() %>%
      imap(~`names<-`(.x, c(names(.x)[1], paste0("pred", .y)))) %>%
      map(~rename_at(.x, vars(1), ~str_remove(., "u:"))) %>%
      modify_if(.x = ., .p = ~str_detect(names(.)[1], ":"), ~separate(data = .x, col = 1, into = c("line_name", "environment"), sep = ":")) %>%
      # Include the testing df, but only the columns (ranefs no int). This makes merging easier
      c(list( distinct(select(test, ranefs_no_int)) ), .) %>%
      reduce(.x = ., .f = left_join) %>%
      ## Sum predictions
      mutate(prediction = rowMeans(select(., contains("pred")))) %>%
      select(ranefs_no_int, prediction)


    ## Merge with test
    test1 <- left_join(x = test, y = fit_pred, by = ranefs_no_int)

    # Add intercept?
    if (add.intercept) {
      test1$prediction <- test1$prediction + intercept
    }

    to_return <- list(predictions = test1)

  }

  # Return
  return(to_return)

}
  





####
## Function to analyze GxE
## 
## Used variance components and comparison of lack of correlation versus variance heterogeneity
####


analyze_gxe <- function(data) {
  
  
  
}
  






## Functions
# A function to generate LOO resamples using a grouped data.frame
crossv_loo2 <- function(data, id = ".id") {
  
  # Get the group keys
  grp_keys <- group_keys(data)
  # Get group indices
  grp_ind <- group_indices(data)
  
  ## Ungroup the data frame
  df <- ungroup(data)
  # Generate integer vector of rows
  rows <- seq(nrow(df))
  
  # For each key, outer join for train and inner join for test
  train_list <- map(seq(unique(grp_ind)), ~resample(data = df, idx = rows[! grp_ind %in% .]))
  test_list <- map(seq(unique(grp_ind)), ~resample(data = df, idx = rows[grp_ind %in% .]))
  
  # Package into tibble
  grp_keys[["train"]] <- train_list
  grp_keys[["test"]] <- test_list
  
  return(grp_keys)
}

# Create a function to generate these relationship matrices
Env_mat <- function(x, method = c("Jarquin2014", "Malosetti2016", "Rincent2019")) {
  
  method <- match.arg(method)
  
  if (method == "Jarquin2014") {
    EMAT <- (tcrossprod(x) / ncol(x))
    
  } else if (method == "Malosetti2016") {
    # Calculate pairwise euclidian distance
    euc_dist <- apply(X = x, MARGIN = 2, FUN = function(z) list( as.matrix(dist(z)) / ( max(z) - min(z) ) ) )
    # Get the first element of each list
    euc_dist <- lapply(euc_dist, "[[", 1)
    
    ## Sum all matrices
    EMAT <- Reduce(f = `+`, x = euc_dist)
    
  } else if (method == "Rincent2019") {
    D <- as.matrix(dist(x))
    EMAT <- 1 - (D / max(D))
    
  }
  
  return(EMAT)
  
}




# ### A function to perform forward stepwise regression
# ### while consider VIF
# step_fwd <- function(object, scope, vif.cutoff = 2) {
#   
#   # New vif function
#   vif1 <- function(object) {
#     if (length(attr(terms(formula(object)), "term.labels")) < 2) {
#       0
#     } else {
#       vif(object)
#     }}
#     
#   
#   ## Copy items
#   objecti <- object
#   add <- TRUE
#   
#   while(add) {
#   
#     # Perform the initial add1
#     add1_i <- add1(object = objecti, scope = scope)
#     add1_i_df <- as.data.frame(add1_i)
#     
#     ## Calculate vif for each of those add1 models
#     add1_models <- str_subset(string = row.names(add1_i), pattern = "none", negate = TRUE) %>% 
#       map(~update(objecti, formula = add_predictors(f = formula(objecti), reformulate(.))))
#     add1_vif <- map(add1_models, vif1)
#     
#     # Determine if there are any VIF > cutoff
#     any_large_vif <- map_lgl(add1_vif, ~any(. > vif.cutoff))
#     # Add this to thhe add1
#     add1_i_df$vif <- c(FALSE, any_large_vif)
#     
#     # Determine the regressor that minimizes AIC without introducting inflated vif
#     to_add <- row.names(add1_i)[which.min(subset(add1_i_df, !vif, AIC, drop = TRUE))]
#   
#     if (length(to_add) == 0 | to_add == "<none>") {
#       add <- FALSE
#       
#     } else {
#       # Add that regressor to the model
#       objecti_temp <- update(objecti, formula = add_predictors(f = formula(objecti), reformulate(to_add)))
#       objecti <- objecti_temp
# 
#     }
#     
#   }
#   
#   
# }
  
    











