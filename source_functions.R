## S2MET Functions
## 
## A script with useful functions used in the prediction analysis of the S2MET
## 



## Other/utility functions
# A function to assign cores to a data.frame
assign_cores <- function(df, n_core) {
  df$core <- sort(rep(seq(n_core), length.out = nrow(df)))
  return(df)
}

# A function to replace a character vector and convert it to a factor
as_replaced_factor <- function(x, replacement) {
  x_repl <- str_replace_all(string = x, pattern = replacement)
  factor(x_repl, levels = replacement)
}




## A generic prediction function that takes training and test data and returns
## PGVs and accuracy
gblup <- function(formula, random, K, train, test, fun = c("rrblup", "sommer"), fit.env = TRUE, bootreps = NULL, add.mu = FALSE) {
  
  if (missing(formula)) formula <- value ~ 1 + environment
  if (missing(random)) random <- ~ line_name
  
  ## separate random and fixed terms
  fixed_terms <- attr(terms(formula), "term.labels")
  random_terms <- attr(terms(random), "term.labels")

  ## Combine fixed and random to one formula
  formula1 <- as.formula(paste(paste(as.character(formula)[c(2,1,3)], collapse = " "), as.character(random)[-1], sep = " + "))
  
  # Create a model.frame
  mf <- model.frame(formula1, weights = std_error, data = train)
  
  
  ## If the number of random terms is > 1, K must be a list
  if (length(random_terms) == 1) {
    stopifnot(levels(mf[[random_terms]]) %in% colnames(K))
    
    K <- list(K)
    names(K) <- random_terms
    
  } else {
    if (!is.list(K)) stop("If the number of random terms is > 1, K must be a list of relationship matrices.")
    
    ## Test names
    name_test <- map2_lgl(.x = random_terms, .y = names(K1), ~.x == .y)
    if (!all(name_test)) stop("If K is a list, the names of the list must match the random terms")
    
  }
   
  
  fun <- match.arg(fun)
  
  
  
  # Vectors and matrices
  y <- model.response(mf)
  
  
  if (nlevels(mf[[fixed_terms]]) <= 1 | !fit.env) {
    X <- model.matrix(~ 1, droplevels(mf))
  
    } else {
    X <- model.matrix(formula, droplevels(mf))
    
  }
  
  ## Random effects
  Z_list <- list()
  
  for (term in random_terms) {
    Zi <- model.matrix(as.formula(paste0("~ -1 + ", term)), mf)
    colnames(Zi) <- colnames(K[[term]])
    Z_list[[term]] <- Zi
    
  }
  
  
  # Split on function
  if (fun == "rrblup") {
    # Can't run with > 1 random term
    stopifnot(length(random_terms) == 1)
    
    fit <- mixed.solve(y = y, Z = Z_list[[1]], K = K[[1]], X = X)
    
    # Extract PGVs
    pgv <- fit$u %>% 
      data.frame(line_name = names(.), pred_value = ., row.names = NULL, stringsAsFactors = FALSE)
    
    beta <- fit$beta[1]
    
    
  } else if (fun == "sommer") {
    
    # R <- solve(diag(mf$`(weights)`^2))
    # fit <- sommer::mmer(Y = y, X = X, Z = list(g = list(Z = Z, K = K)), R = list(res = R), silent = TRUE)
    
    random_list <- map2(Z_list, K, ~list(Z = .x, K = .y))
    fit <- sommer::mmer(Y = y, X = X, Z = random_list, silent = TRUE)
  
    # Extract PGVs
    pgv <- fit$u.hat
    pgv <- data.frame(line_name = row.names(pgv), pred_value = pgv[,1], row.names = NULL, stringsAsFactors = FALSE)
    
    beta <- fit$beta.hat[1]
    
  }
  
  if (add.mu) pgv$pred_value <- pgv$pred_value + beta
  
  # If test is missing, just return the predictions
  if (missing(test)) {
    
    comb <- pgv
    acc <- boot <- NA
    
  } else {
  
    # Combine the PGVs with the phenotypic observations and calculate accuracy
    comb <- left_join(test, pgv, by = "line_name")
    acc <- cor(comb$value, comb$pred_value)
    
    # Bootstrap if replicates are provided
    if (!is.null(bootreps)) {
      boot <- bootstrap(x = comb$value, y = comb$pred_value, fun = "cor", boot.reps = bootreps)
    } else {
      boot <- NA
    }
    
  }
  
  # Return a list
  list(accuracy = acc, pgv = comb, boot = as_data_frame(boot))
}





## Prediction functions for cross-validation
## Functions for each model
## 
## Model 1 - fixed environment, random genotype
model1 <- function(train, test, Kg) {
  
  mf <- model.frame(value ~ line_name + environment, data = train)
  
  # Matrices
  y <- model.response(mf)
  X <- model.matrix(~ 1 + environment, data = mf)
  Zg <- model.matrix(~ -1 + line_name, mf)
  colnames(Zg) <- colnames(Kg)
  
  fit <- mixed.solve(y = y, Z = Zg, K = Kg, X = X)

  # Extract PGVs
  pgv <- fit$u %>% 
    data.frame(line_name = names(.), pred_value = ., row.names = NULL)
  
  ## Measure accuracy
  # Combine the PGVs with the phenotypic observations and calculate accuracy
  comb <- left_join(test, pgv, by = "line_name") %>%
    select(environment, line_name, value, pred_value)
  acc <- cor(comb$value, comb$pred_value)
  
  # Return a list
  list(accuracy = acc, pgv = comb)
  
}

# Model 2 - fixed, environment, random genotype, random GxE
model2 <- function(train, test, Kg) {
  
  mf <- train %>%
    mutate(environment = factor(environment, levels = unique(c(unique(train$environment), unique(test$environment))))) %>%
    model.frame(value ~ line_name + environment, data = .)
  
  # Matrices
  y <- model.response(mf)
  X <- model.matrix(~ 1 + environment, data = droplevels(mf))
  Zg <- model.matrix(~ -1 + line_name, mf)
  colnames(Zg) <- colnames(Kg)
  Ze <- model.matrix(~ -1 + environment, mf)
  Zge <- model.matrix(~ -1 + line_name:environment, mf)
  
  ## Covariance matrix for E
  E <- diag(nlevels(mf$environment)); dimnames(E) <- replicate(2, levels(mf$environment), simplify = FALSE)
  Kge <- kronecker(E, Kg, make.dimnames = TRUE)
  colnames(Zge) <- colnames(Kge)
  
  # Fit
  fit <- sommer::mmer(Y = y, X = X, Z = list(g = list(Z = Zg, K = Kg), ge = list(Z = Zge, K = Kge)), silent = TRUE)

  
  # Extract PGVs
  g <- fit$u.hat$g %>%
    as.data.frame() %>%
    rownames_to_column("line_name") %>%
    rename(g = T1)
  
  ge <- fit$u.hat$ge %>% 
    as.data.frame() %>%
    rownames_to_column("term") %>%
    separate(col = term, c("environment", "line_name"), sep = ":") %>%
    rename(ge = T1)
  
  pgv <- full_join(g, ge, by = "line_name") %>% 
    mutate(pred_value = g + ge)
  
  ## Measure accuracy
  # Combine the PGVs with the phenotypic observations and calculate accuracy
  comb <- left_join(test, pgv, by = c("environment", "line_name")) %>%
    select(environment, line_name, value, pred_value)
  acc <- cor(comb$value, comb$pred_value)
  
  # Return a list
  list(accuracy = acc, pgv = comb)
  
}





  
