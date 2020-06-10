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




# ## A model fitting function
# ## 
# ## All models are random
# ## 
# ## Model1: y = G
# ## Model2: y = G + E
# ## Model3: y = G + E + GE
# ## 
# ## Note: for single traits only
# ## 
# 
# predict_gv <- function(train, test, model = c("model1", "model2", "model3"), relMat.list, covariate.list = list(NULL),
#                        object, add.fixed = FALSE, verbose = FALSE) {
# 
#   ## Depending on the model, make sure correct relmats are present
#   model <- match.arg(model)
# 
#   ## Assign required relmats
#   req_relmat <- switch(
#     model,
#     model1 = c("G"),
#     model2 = c("G", "E"),
#     model2a = c("G"),
#     model2b = c("G"),
#     model3 = c("G", "E", "GE"),
#     model3a = c("G", "GE"),
#     model3b = c("G", "GE")
#   )
# 
#   ## Assign ranefs
#   ranefs <- switch(
#     model,
#     model1 = c("line_name"),
#     model2 = c("line_name", "environment"),
#     model2a = c("line_name"),
#     model2b = c("line_name"),
#     model3 = c("line_name", "environment", "line_name:environment"),
#     model3a = c("line_name", "line_name:environment"),
#     model3b = c("line_name", "line_name:environment")
#   )
# 
#   ## Assign fixefs
#   fixef <- switch(
#     model,
#     model1 = NULL,
#     model2 = NULL,
#     model2a = c("environment"),
#     model2b = NULL,
#     model3 = NULL,
#     model3a = c("environment"),
#     model3b = NULL
#   )
# 
#   ## Assign covariates
#   covariates <- switch(
#     model,
#     model1 = NULL,
#     model2 = NULL,
#     model2a = NULL,
#     model2b = "E",
#     model3 = NULL,
#     model3a = NULL,
#     model3b = "E"
#   )
# 
# 
#   ## Error handling
# 
#   ## Check for presence of relmats
#   stopifnot(is.list(relMat.list))
#   stopifnot(all(req_relmat %in% names(relMat.list)))
# 
#   ## Check for presence of covariates
#   stopifnot(is.list(covariate.list))
#   stopifnot(all(covariates %in% names(covariate.list)))
# 
# 
#   ## Extract covariates
#   covariates <- unlist(covariate.list[covariates])
# 
# 
#   ## Create fixed/random formula
#   fixed_formula <- reformulate(termlabels = c("1", c(fixef, covariates)), response = "value")
# 
# 
#   # If object is missing, do not include contraints on the random effects
#   if (missing(object)) {
#     random_formula <- reformulate(termlabels = mapply(req_relmat, ranefs, FUN = function(.x, .y) {
#       paste0("vs(", .y, ", Gu = relMat.list$", .x, ")") }))
# 
# 
#     # Residual formula
#     rcov_form <- ~ units
# 
#     # If not missing, extract variance components and impose restrictions
#   } else {
#     # Error handling
#     # 'object' must have element called sigma
#     stopifnot("sigma" %in% names(object))
# 
#     # Extract varcomps
#     object_varcomps <- object$sigma
#     names(object_varcomps) <- gsub(pattern = "u:", replacement = "", x = names(object_varcomps))
#     object_varcomps <- object_varcomps[c(ranefs, "units")]
# 
#     # Create formula
#     random_formula <- reformulate(termlabels = mapply(req_relmat, ranefs, FUN = function(.x, .y) {
#       paste0("vs(", .y, ", Gu = relMat.list$", .x, ", Gt = object_varcomps[['", .y, "']], Gtc = fixm(1))") }) )
# 
#     # Residual formula
#     rcov_form <- ~ vs(units, Gt = object_varcomps$units, Gtc = fixm(1))
#   }
# 
# 
#   ## Create a model frame
#   full_formula <- add_predictors(fixed_formula, reformulate(ranefs))
#   mf <- model.frame(formula = full_formula, train, drop.unused.levels = FALSE)
# 
#   # X matrix
#   X <- model.matrix(fixed_formula, train)
# 
#   # List of Z matrices
#   Zlist <- lapply(ranefs, function(term) model.matrix(reformulate(c(-1, term)), mf))
#   names(Zlist) <- ranefs
# 
#   ## Fit the model
#   fit <- mmer(fixed = fixed_formula, random = random_formula, rcov = rcov_form, data = mf, date.warning = FALSE, verbose = verbose)
# 
#   # If interaction present, split and take those ranefs
#   if (str_detect(ranefs, ":")) {
#     ranefs_no_int <- str_split(str_subset(ranefs, ":"), pattern = ":")[[1]]
#     
#   } else {
#     ranefs_no_int <- str_subset(ranefs, ":", negate = TRUE)
#     
#   }
# 
#   ## Figure out what to return
#   # If test is missing, return the object
#   # Create a list to return
#   if (missing(test)) {
#     to_return <- list(object = fit)
# 
#   } else {
# 
#     ## Get the random effects
#     U <- lapply(X = randef(fit), FUN = "[[", "value")
#     # Rename
#     names(U) <- str_remove(string = names(U), pattern = "u:")
# 
#     ## Use incidence matrices to predict random effs
#     Zu_list <- mapply(Zlist, U, FUN = function(Z, u) {
#       Zu <- Z %*% u
#       `row.names<-`(Zu, names(u)[apply(X = Z, MARGIN = 1, FUN = function(row) which(row == 1))])
#     }, SIMPLIFY = FALSE)
# 
#     ## Find unique values and convert to df
#     Zu_df <- lapply(Zu_list, function(x) as.data.frame(unique(x))) %>%
#       map(~rownames_to_column(., "term")) %>%
#       list(., names(.), seq_along(.)) %>%
#       pmap(., ~`names<-`(..1, c(..2, paste0("pred", ..3)))) %>%
#       modify_if(.x = ., .p = ~str_detect(names(.)[1], ":"), ~separate(data = .x, col = 1, into = c("line_name", "environment"), sep = ":"))
#     
#     # If any dfs in the list have > 2 columns, sort and merge. 
#     # Otherwise combine and full join
#     if (any(map_lgl(Zu_df, ~ncol(.) > 2))) {
#       Zu_join <- reduce(.x = Zu_df[order(map_dbl(Zu_df, ncol), decreasing = TRUE)], .f = left_join)
#       
#     } else {
#       combn <- Zu_df[ranefs_no_int] %>% 
#         map(select, 1) %>% 
#         reduce(crossing)
#       Zu_join <- reduce(c(list(combn), Zu_df), left_join)
#       
#     }
#     
#     ## Sum for predictions
#     rand_pred <- Zu_join %>%
#       mutate(prediction = rowMeans(select(., contains("pred")))) %>%
#       select(ranefs_no_int, prediction)
#     
#     
#     
#     #### Fixed effects ####
#     beta <- as.matrix(column_to_rownames(fit$Beta[,-1], "Effect"))
#     Xb <- X %*% beta
#     
#     
#     
#     
# 
#     ## Extract predictions
#     fit_pred <- randef(fit) %>%
#       map("value") %>%
#       map(~tibble(term = names(.x), pred = .x), .vars = vars(term)) %>%
#       imap(~`names<-`(.x, c(.y, "pred"))) %>%
#       unname() %>%
#       imap(~`names<-`(.x, c(names(.x)[1], paste0("pred", .y)))) %>%
#       map(~rename_at(.x, vars(1), ~str_remove(., "u:"))) %>%
#       modify_if(.x = ., .p = ~str_detect(names(.)[1], ":"), ~separate(data = .x, col = 1, into = c("line_name", "environment"), sep = ":")) %>%
#       # Include the testing df, but only the columns (ranefs no int). This makes merging easier
#       c(list( distinct(select(test, ranefs_no_int)) ), .) %>%
#       reduce(.x = ., .f = left_join) %>%
#       ## Sum predictions
#       mutate(prediction = rowMeans(select(., contains("pred")))) %>%
#       select(ranefs_no_int, prediction)
# 
# 
#     ## Merge with test
#     test1 <- left_join(x = test, y = fit_pred, by = ranefs_no_int)
# 
#     # Add intercept?
#     if (add.intercept) {
#       test1$prediction <- test1$prediction + intercept
#     }
# 
#     to_return <- list(predictions = test1)
# 
#   }
# 
#   # Return
#   return(to_return)
# 
# }
  

## Functions to calculate heritability from lmer models
## Function for heritability using genotype difference basis
## 
## This code is from https://github.com/PaulSchmidtGit/Heritability
## 
herit2 <- function(object, ...) {
  UseMethod("herit2")
}


herit2.lmerModLmerTest <- function(object, gen.var = "line_name", type = c("cullis", "piepho")) {
  herit2.merMod(object = object, gen.var = gen.var, type = type)
}

herit2.merMod <- function(object, gen.var = "line_name", type = c("cullis", "piepho")) {
  
  ## Match type argument
  type <- match.arg(type)
  # Make sure gen.var is in the terms
  if (!any(grepl(pattern = gen.var, x = attr(terms(formula(object)), "term.labels")))) {
    stop(paste0(gen.var, " is not among the terms in the model object."))
  }
  
  ## Split based on type
  if (type == "cullis") {
    
    # extract estimated variance components (vc)
    vc   <- as.data.frame(VarCorr(object)) 
    
    # R = varcov-matrix for error term
    n  <- length(summary(object)$residuals) # numer of observations
    vcR <- subset(vc, grp == "Residual", vcov, drop = TRUE)     # error vc
    R <- diag(n) * vcR                   # R matrix = I * vcR
    
    # G = varcov-matrx for all random effects
    # varcov-matrix for genotypic effect
    n_g  <- summary(object)$ngrps[gen.var]    # number of genotypes
    vcG <- subset(vc, grp == gen.var, vcov, drop = TRUE)         # genotypic vc
    G  <- diag(n_g) * vcG               # gen part of G matrix = I * vc.g
    
    # Design Matrices
    X <- as.matrix(getME(object, "X")) # Design matrix fixed effects
    Ztlist <- (getME(object, "Ztlist")) # Design matrix random effects
    Zind <- which(! grepl(pattern = ":", x = names(Ztlist)))
    Z <- t(as.matrix(Ztlist[Zind][[which(grepl( pattern = gen.var, x = names(Ztlist)[Zind] ))]]))
    
    # Inverse of R
    Rinv <- solve(R)
    tryGinv <- try(Ginv <- solve(G), silent = TRUE)
    Ginv <- if (class(tryGinv) == "try-error") `dimnames<-`(MASS::ginv(G), dimnames(G)) else Ginv
    
    # Mixed Model Equation (HENDERSON 1986; SEARLE et al. 2006)
    C11 <- t(X) %*% Rinv %*% X
    C12 <- t(X) %*% Rinv %*% Z
    C21 <- t(Z) %*% Rinv %*% X
    C22 <- t(Z) %*% Rinv %*% Z + Ginv
    
    C <- as.matrix(rbind(cbind(C11, C12),  # Combine components into one matrix C
                         cbind(C21, C22)))
    
    ## Get the levels of gen.var ##
    gen.var.levels <- levels(model.frame(object)[[gen.var]])
    
    # Mixed Model Equation Solutions 
    tryInverse <- try(C.inv <- solve(C), silent = TRUE) # Inverse of C
    # If error, use generalized inverse
    C.inv <- if (class(tryInverse) == "try-error") `dimnames<-`(MASS::ginv(C), dimnames(C)) else C.inv
    
    C22.g <- C.inv[gen.var.levels, gen.var.levels] # subset of C.inv that refers to genotypic BLUPs
    
    # Mean variance of BLUP-difference from C22 matrix of genotypic BLUPs
    one <- matrix(1, nrow = n_g) # vector of 1s
    P.mu <- diag(n_g, n_g) - one %*% t(one)  # P.mu = matrix that centers for overall-mean
    vdBLUP.sum <- sum(diag((P.mu %*% C22.g)))        # sum of all variance of differences = trace of P.mu*C22.g
    vdBLUP.avg <- vdBLUP.sum * (2/(n_g*(n_g-1)))   # mean variance of BLUP-difference = divide sum by number of genotype pairs
    
    ## Return heritability
    H2 <- 1 - (vdBLUP.avg / 2 / vcG)
    # Set to NA if infinite
    H2 <- ifelse(is.infinite(H2), NA, H2)
    
  } else if (type == "piepho") {
    
    
    ## Refit the model with gen.var as fixed
    rhs <- as.character(formula(object))[3]
    new_formula <- reformulate(gsub(pattern = paste0("\\(1 \\| ", gen.var, "\\)"), replacement = gen.var, x = rhs), 
                               response = "value")
    
    ## Fit depending if mixed or fixed
    if (is.null(findbars(new_formula))) {
      object_f <- lm(formula = new_formula, data = model.frame(object))
      
    } else {
      object_f <- lmer(formula = new_formula, data = model.frame(object))
      
    }
    
    
    ## Get the variance components
    vc   <- as.data.frame(VarCorr(object)) 
    
    # Genotypic variance component
    vcG <- subset(vc, grp == gen.var, vcov, drop = TRUE)
    
    # Obtaining adjusted means based on genotypic BLUEs
    diffs_BLUE <- as.data.frame(emmeans(object = object_f, reformulate(gen.var, "pairwise"))$contrasts) 
    vdBLUE_avg <- mean(diffs_BLUE$SE^2) # mean variance of a difference = mean squared standard error of a difference
    
    # Return the heritabiliy
    H2 <- vcG/(vcG + vdBLUE_avg/2)
    
    
  }
  
  # Return heritability
  return(H2)
  
}


## A function to add or remove covariates based on AIC and vif
step_vif <- function(object, data, scope, direction, vif.threshold) {
  
  ## Get the terms to keep
  terms_keep <- terms_hold <- attr(terms(scope$lower), "term.labels")
  
  # Forward
  if (direction == "forward") {
    
    # Vector of all terms to add
    all_terms_add <- setdiff(attr(terms(scope$upper), "term.labels"), attr(terms(scope$lower), "term.labels"))
    terms_add <- all_terms_add
    
    ## Base formula AIC (no terms added)
    base_model <- object
    base_model_AIC <- AIC(base_model)
    
    # List of model formulas
    formula_list <- map(terms_add, ~reformulate(termlabels = c(terms_hold, .), response = "value")) %>%
      setNames(., terms_add)
    
    ## Fit models one-at-a-time
    fit_list <- map(formula_list, lm, data = data)
    
    # Get AIC and vif
    tidy_diag <- imap_dfr(fit_list, ~{
      vif_x <- vif(.x)
      colnames(vif_x) <- c("gvif", "df", "adj.gvif")
      as.data.frame(cbind(vif_x[.y,, drop = FALSE], AIC = AIC(.x)))
    }) %>% mutate(term = names(fit_list))
    
    ## Should the term be added?
    # Filter on suitable vif and AIC; then find lowest AIC
    tidy_diag1 <- subset(tidy_diag, AIC < base_model_AIC & gvif < vif.threshold)
    
    # Which term
    term_retain <- tidy_diag1$term[which.min(tidy_diag1$AIC)]
    # Which model to keep
    model_keep_which <- which(names(fit_list) == term_retain)
    
    
    ## While loop
    while (length(model_keep_which) > 0) {
      
      ## Set the new base model
      base_model <- fit_list[[model_keep_which]]
      base_model_AIC <- AIC(base_model)
      
      ## Remove the 'term_retain' from the list of terms to add
      terms_add <- setdiff(terms_add, term_retain)
      terms_hold <- c(terms_hold, term_retain)
      
      # List of model formulas
      formula_list <- map(terms_add, ~reformulate(termlabels = c(terms_hold, .), response = "value")) %>%
        setNames(., terms_add)
      
      ## Fit models one-at-a-time
      fit_list <- map(formula_list, lm, data = data)
      
      # Get AIC and vif
      tidy_diag <- imap_dfr(fit_list, ~{
        vif_x <- vif(.x)
        colnames(vif_x) <- c("gvif", "df", "adj.gvif")
        as.data.frame(cbind(vif_x[.y,, drop = FALSE], AIC = AIC(.x)))
      }) %>% mutate(term = names(fit_list))
      
      ## Should the term be added?
      # Filter on suitable vif and AIC; then find lowest AIC
      tidy_diag1 <- subset(tidy_diag, AIC < base_model_AIC & gvif < vif.threshold)
      
      # Which term
      term_retain <- tidy_diag1$term[which.min(tidy_diag1$AIC)]
      # Which model to keep
      model_keep_which <- which(names(fit_list) == term_retain)
      
    }
    
    # Model to return
    model_return <- base_model
    
  } # IF statement
  
  # Return the model
  return(model_return)
  
}





## Function to perform factorial regression analysis to identify covariates
## 
## This function will fit three models - a base model, base + main covariates, 
## and base + main covariates + interaction covariates
## 
## It uses stepwise elimination or aprior expectations to determine which
## covariates to keep
## 
## step.vif used variance inflation factors to remove covariates
## 
## step1 removes related covariates (i.e. tmean, tmin, tmax) based on 
## the first pass of mean squares
## 
## 
fact_reg <- function(base.formula = value ~ line_name, data, gen = "line_name", 
                     env = "environment", covariates, method = c("apriori", "step", "step1"),
                     criterion = c("AIC", "BIC")) {
  
  # Match args
  method <- match.arg(method)
  criterion <- match.arg(criterion)
  k <- ifelse(criterion == "AIC", 2, log(nrow(data)))
  
  # Fit the base models
  fit1 <- lm(formula = base.formula, data = data)
  fit1_alt <- update(fit1, formula = add_predictors(formula(fit1), reformulate(env)))
  
  # Split by apriori or step
  if (method == "apriori") {
  
    # Fit line name plus main covariates
    form2 <- add_predictors(formula(fit1), reformulate(covariates))
    fit2 <- update(object = fit1, formula = form2)
    
    # Fit line name + E-covariate + GE-covariates
    form3 <- add_predictors(form2, reformulate(paste0(paste0(gen, ":"), covariates)))
    fit3 <- update(object = fit2, formula = form3)

  } else if (method %in% c("step", "step1")) {
    
    # Determine the max degrees of freedom for covariates
    J <- length(unique(data[[env]]))
    H <- J - 1
    
    # Fit line name plus E-covariates
    # Add 1 covariate at-a-time and measure meansq (+ p.value)
    form_list <- map(covariates, ~add_predictors(formula(fit1), reformulate(.)))
    fit2_list <- fit_with(data, lm, form_list)
    fit2_anova_list <- map(fit2_list, anova)
    
    # Get a table of mean squares for each covariate
    fit2_covariate_meansq <- map_df(fit2_anova_list, tidy) %>% 
      filter(term %in% covariates) %>%
      # Sort by meansq
      arrange(desc(meansq)) %>%
      # Adjust p-values
      mutate(p.adj = p.adjust(p = p.value, method = "bonf")) %>%
      filter(p.adj < alpha)
    
    ## If step1, select covariates individually per growth stage
    if (method == "step1") {
      fit2_covariate_meansq_base <- fit2_covariate_meansq
      
      fit2_covariate_meansq_temp <- fit2_covariate_meansq_base %>% 
        filter(str_detect(term, "mint|maxt|tmean")) %>% 
        mutate(stage = str_extract(term, "flowering|grain_fill|early_vegetative|late_vegetative")) %>% 
        group_by(stage) %>% 
        top_n(x = ., n = 1, wt = meansq)
      
      fit2_covariate_meansq <- bind_rows(
        filter(fit2_covariate_meansq_base, str_detect(term, "mint|maxt|tmean", negate = TRUE)),
        fit2_covariate_meansq_temp
      ) 
      
    }
    
    
    
    # If there are no significant covariates, then assign fit2_fwd as fit1
    if (nrow(fit2_covariate_meansq) == 0) {
      fit2_fwd <- fit1
      
    } else {
      ## Significant covariates define the scope for forward regression
      scope <- list(lower = formula(fit1), upper = add_predictors(formula(fit1), reformulate(fit2_covariate_meansq$term)))
      fit2_fwd <- step(object = fit1, scope = scope, direction = "forward", trace = 0)
      # fit2_fwd_anova <- Anova(fit2_fwd, type = "II")
      
    }
    
    
    ## Did the model exhaust the number of possible covariates?
    ECs <- setdiff(attr(terms(fit2_fwd), "term.labels"), "line_name")
    # If so, drop the last one
    if (length(ECs) > H) {
      
      # Determine the interaction terms in this model
      drop_scope <- reformulate(ECs)
      
      fit2_drop <- drop1(fit2_fwd, scope = drop_scope)
      # Find the term to drop
      term_drop <- tidy(fit2_drop) %>%
        filter(is.finite(AIC)) %>%
        arrange(AIC) %>%
        slice(1) %>%
        pull(term)
      
      fit2_drop_form <- terms(fit2_fwd) %>%
        drop.terms(termobj = ., dropx = which(attr(., "term.labels") == term_drop), keep.response = TRUE) %>% 
        formula()
      
      # Refit
      fit2_fwd <- update(object = fit2_fwd, formula = fit2_drop_form)
      
    } 
    
    
    
    # Fit line name + E-covariate + GE-covariates
    # Perform the same procedure
    form_list <- map(covariates, ~add_predictors(formula(fit2_fwd), reformulate(paste0(paste0(gen, ":"), .))))
    fit3_list <- fit_with(data, lm, form_list)
    # fit3_anova_list <- map(fit3_list, ~Anova(., type = "II"))
    fit3_anova_list <- map(fit3_list, anova)
    
    
    # Get a table of mean squares for each covariate
    fit3_covariate_meansq <- map_df(fit3_anova_list, tidy) %>% 
      filter(str_detect(term, ":")) %>%
      mutate(meansq = sumsq / df,
             # Adjust p-values
             p.adj = p.adjust(p = p.value, method = "bonf")) %>%
      filter(p.adj < alpha) %>%
      # Sort by meansq
      arrange(desc(meansq))
    
    
    
    ## If step1, select covariates individually per growth stage
    if (method == "step1") {
      fit3_covariate_meansq_base <- fit3_covariate_meansq
      
      fit3_covariate_meansq_temp <- fit3_covariate_meansq_base %>% 
        filter(str_detect(term, "mint|maxt|tmean")) %>% 
        mutate(stage = str_extract(term, "flowering|grain_fill|early_vegetative|late_vegetative")) %>% 
        group_by(stage) %>% 
        top_n(x = ., n = 1, wt = meansq)
      
      fit3_covariate_meansq <- bind_rows(
        filter(fit3_covariate_meansq_base, str_detect(term, "mint|maxt|tmean", negate = TRUE)),
        fit3_covariate_meansq_temp
      ) 
      
    }
    
    
    
    # If there are no significant covariates, then assign fit2_fwd as fit1
    if (nrow(fit3_covariate_meansq) == 0) {
      fit3_fwd <- fit2_fwd
      
    } else {
      
      ## Significant covariates define the scope for forward regression
      scope <- list(lower = formula(fit2_fwd), upper = add_predictors(formula(fit2_fwd), reformulate(fit3_covariate_meansq$term)))
      fit3_fwd <- step(object = fit2_fwd, scope = scope, direction = "forward", trace = 0) 
      # fit3_fwd_anova <- Anova(fit3_fwd, type = "II")
      
    }
    
    
    ## Did the model exhaust the number of possible covariates?
    ECs <- attr(terms(fit3_fwd), "term.labels") %>%
      str_subset(., ":")
    # If so, drop the last one
    if (length(ECs) > H) {
      
      # Determine the interaction terms in this model
      drop_scope <- reformulate(ECs)
      
      fit3_drop <- drop1(fit3_fwd, scope = drop_scope)
      # Find the term to drop
      term_drop <- tidy(fit3_drop) %>%
        filter(is.finite(AIC)) %>%
        arrange(AIC) %>%
        slice(1) %>%
        pull(term)
      
      fit3_drop_form <- terms(fit3_fwd) %>%
        drop.terms(termobj = ., dropx = which(attr(., "term.labels") == term_drop), keep.response = TRUE) %>% 
        formula()
      
      # Refit
      fit3_fwd <- update(object = fit3_fwd, formula = fit3_drop_form)
      
    } 
    
    
    
    # Rename
    fit2 <- fit2_fwd
    fit3 <- fit3_fwd
    
  } 

  # Return a list of model fits
  list(fit1 = fit1, fit1_alt = fit1_alt, fit2 = fit2, fit3 = fit3)
    
}





# rapid_cv <- function(object, index, rapid = TRUE) {
# 
#   # Get the model frame
#   mf <- model.frame(object)
#   # Model matrix and crossprod inverse
#   X <- model.matrix(object)
#   XtX_inv <- solve(crossprod(X))
#   # Response
#   y <- model.response(mf)
#   # Coefficients
#   beta_hat <- coef(object)
#   
#   ## Turn index into the holdout set
#   all_index <- seq_along(y)
#   holdout_index <- lapply(index, setdiff, x = all_index)
# 
#   if (rapid) {
# 
#     # Iterate over indices
#     cv_out <- lapply(X = index, FUN = function(modelInd) {
#       # Subset X for the d testing datapoints
#       X_d <- X[-modelInd,,drop = FALSE]
#       d <- nrow(X_d)
#       # H matrix
#       H_d <- tcrossprod(X_d %*% XtX_inv, X_d)
#       H_d_inv <- ginv(diag(d) - H_d)
#       # Residuals
#       e_d <- y[-modelInd] - X_d %*% beta_hat
# 
#       # New betas
#       beta_hat_holdout <- beta_hat - (tcrossprod(XtX_inv, X_d) %*% H_d_inv %*% e_d)
# 
#       # yhat
#       y_hat <- X_d %*% beta_hat_holdout
#       # return
#       y_hat
# 
#     })
# 
#   } else {
# 
#     # Iterate over indices
#     cv_out <- lapply(X = index, FUN = function(modelInd) {
#       # Subset X for the training and testing
#       Xtrain <- X[modelInd,,drop = FALSE]
#       Xtest <-X[-modelInd,,drop = FALSE]
#       ytrain <- y[modelInd]
# 
#       # Refit
#       beta_hat <- lm.fit(x = Xtrain, y = ytrain)$coefficients
#       isNAbeta <- which(!is.na(beta_hat))
#       # Predict
#       ytest <- Xtest[,isNAbeta,drop = FALSE] %*% beta_hat[isNAbeta]
#       # Return
#       ytest
# 
#     })
# 
# 
#   }
# 
#   ## Pull out y
#   y_hat_all <- unlist(cv_out)
# 
#   # Calculate metrics and return
#   R2 <- cor(y_hat_all, y)^2
#   mse <- mean((y - y_hat_all)^2)
# 
#   c(R2 = R2, MSE = mse, RMSE = sqrt(mse))
# 
# }



## A function to implement a recursive feature addition algorithm that includes
## featurs on the basis of the RMSE of LOO cross-validation
## 
## params:
## x - a matrix of features
## y - a response
## 
## 
## 
rfa_loo <- function(object, data, scope, metric = c("RMSE", "R2"), index, env.col = "environment") {
  
  # Match arguments
  metric <- match.arg(metric)
  maximize <- ifelse(metric == "RMSE", FALSE, TRUE)
  
  # Assuming environments are here, find the maximum number of covariates that could be used
  J <- n_distinct(data[[env.col]])
  H <- J - 2
  
  ## Define the model frame for the upper scope
  mf <- model.frame(scope$upper, data = data)
  
  ## Perform the base cv
  base_fit <- object
  base_resample_pred <- rapid_cv(object = object, index = index)
  
  ## Start the algorithm
  # Define the available covariates and the used covariates
  cov_available <- setdiff(attr(terms(scope$upper), "term.labels"), attr(terms(scope$lower), "term.labels"))
  cov_used <- intersect(attr(terms(scope$upper), "term.labels"), attr(terms(scope$lower), "term.labels"))
  
  # Define a flag to stop the algorithm
  optMetric <- FALSE
  allCovUsed <- length(cov_available) == 0
  
  # A list to store output
  step_out <- list()
  i = 1
  
  # While loop
  while (!optMetric & !allCovUsed) {
    
    # Iterate over cov_available and build models
    cov_addition_objects <- map(cov_available, ~update(object = base_fit, formula = reformulate(c(cov_used, .), response = "value")))
    # Remove NAs
    object_is_na <- map_lgl(cov_addition_objects, ~any(is.na(coef(.x))))
    cov_addition_objects1 <- setNames(object = cov_addition_objects[!object_is_na], cov_available[!object_is_na])
    
    # Run cv - rapid if no interactions
    has_interactions <- map_lgl(cov_addition_objects1, ~any(str_detect(attr(terms(.x), "term.labels"), ":")))
    # cov_addition_pred <- map2(cov_addition_objects1, has_interactions, ~rapid_cv(object = .x, index = index, rapid = !.y))
    cov_addition_pred <- map(cov_addition_objects1, ~rapid_cv(object = .x, index = index))
    

    # Analyze the metrics
    cov_addition_metrics <- t(do.call("rbind", cov_addition_pred))
    cov_addition_metrics <- t(cbind(cov_addition_metrics, `(none)` = base_resample_pred))
    cov_addition_metrics <- cov_addition_metrics[order(cov_addition_metrics[,metric], decreasing = maximize),, drop = FALSE]
    
    # Select the covariate to add
    cov_retain <- row.names(cov_addition_metrics)[1]
    
    # If "none", stop
    if (cov_retain == "(none)") {
      optMetric <- TRUE
      
    } else {
      # Edit cov_available and cov_used
      cov_used <- union(cov_used, cov_retain)
      cov_available <- setdiff(cov_available, cov_used)
      
      allCovUsed <- length(cov_available) == 0 | length(cov_used) - 1 == H
      
      # Edit the metrics to compare
      base_resample_pred <- cov_addition_metrics[1,]
    }
    
    step_out[[i]] <- cov_addition_metrics
    i <- i + 1
    
  } # End the while loop
  
  # Return a list
  list(optVariables = cov_used, finalResults = step_out[[length(step_out)]][1,], stepTestResults = step_out)
  
}








## Function for recursive feature elimination using PLS
rfs_loo <- function(X, Y, comp = 3, metric = c("RMSE", "R2"), direction = c("forward", "backward"), index) {
  
  require(pls)
  
  # Match arguments
  metric <- match.arg(metric)
  direction <- match.arg(direction)
  
  stopifnot(is.matrix(Y))
  stopifnot(is.matrix(X))
  
  yall <- scale(Y)
  # Impute missing with 0
  yall[is.na(yall)] <- 0
  Xall <- scale(X)
  
  # Determine if Y is multivariate
  is.multivariate <- ncol(yall) > 1
  
  # Split on direction
  if (direction == "forward") {
    
    # Split based on multivariate-ness
    if (!is.multivariate) {
      
      X1 <- matrix(1, ncol = 1, nrow = nrow(Xall))
      ## Find the RMSE of the model with all variables
      base_fit <- plsr(yall ~ X1, validation = "LOO", scale = FALSE, center = FALSE)
      rmse <- rmseBase <- RMSEP(object = base_fit)$val[,,1]["CV"]
      r2 <- r2Base <- R2(base_fit)$val[,,1]
      
    } else {
      
      # Perform the resampling
      base_resample_pred <- lapply(X = index, FUN = function(modelInd) {
        X1 <- Xall[modelInd,,drop = FALSE]
        y1 <- yall[modelInd,,drop = FALSE]
        xTest <- Xall[-modelInd,,drop = FALSE]
        yTest <- yall[-modelInd,,drop = FALSE]
        fit <- plsr(y1 ~ X1)
        # Predict and return
        t(rbind(yTest, predict(object = fit, newdata = xTest, comp = seq_len(comp))))
      })
      
      # Row bind
      out <- do.call("rbind", base_resample_pred)
      # Calculate the metrics
      res <- out[,1] - out[,2]
      rmse <- rmseBase <- sqrt(mean(res^2, na.rm = TRUE))
      r2 <- r2Base <- cor(out)[1,2]^2
      
    }
    
    
    ## Start the algorithm
    
    # Set Xbase
    Xbase <- Xall[,NULL]
    # Define the available covariates and the used covariates
    cov_available <- colnames(Xall)
    cov_used <- character()
    
    # Define a flag to stop the algorithm
    optMetric <- FALSE
    allCovUsed <- length(cov_available) == 0
    
    # A list to store output
    step_out <- list()
    i = 1
    
    # While loop
    while (!optMetric & !allCovUsed) {
      
      # Iterate over cov_todrop and build models
      cov_pred <- sapply(X = cov_available, FUN = function(cov_add) {
        
        Xcov <- cbind(Xbase, Xall[,cov_add, drop = FALSE])
        
        # Split based on multivariate-ness
        if (!is.multivariate) {
          
          ## Find the RMSE of the model with all variables
          cov_fit <- plsr(yall ~ Xcov, validation = "LOO", scale = FALSE, center = FALSE)
          rmseCov <- RMSEP(object = cov_fit)$val[,,min(comp+1, ncol(Xcov)+1)]["CV"]
          r2Cov <- R2(cov_fit)$val[,,min(comp+1, ncol(Xcov)+1)]
          
        } else {
          
          # Resampling
          cov_resample_pred <- lapply(X = index, FUN = function(modelInd) {
            X1 <- Xcov[modelInd,,drop = FALSE]
            y1 <- yall[modelInd,,drop = FALSE]
            xTest <- Xcov[-modelInd,,drop = FALSE]
            yTest <- yall[-modelInd,,drop = FALSE]
            fit <- plsr(y1 ~ X1)
            # Predict and return
            t(rbind(yTest, predict(object = fit, newdata = xTest, comp = seq_len(min(comp, ncol(X1))))))
          })
          
          # Row bind
          out <- do.call("rbind", cov_resample_pred)
          # Calculate the metrics
          resCov <- out[,1] - out[,2]
          rmseCov <- sqrt(mean(resCov^2, na.rm = TRUE))
          r2Cov <- cor(out)[1,2]^2
          
        }
        
        # Return a matrix
        c(RMSE = unname(rmseCov), R2 = r2Cov)
        
      })
      
      # Analyze the metrics
      cov_metrics <- t(cbind(cov_pred, `(none)` = c(rmseBase, r2Base)))
      cov_metrics <- cov_metrics[order(cov_metrics[,metric], decreasing = FALSE),, drop = FALSE]
      
      # Select the covariate to add
      cov_retain <- row.names(cov_metrics)[1]
      
      # If "none", stop
      if (cov_retain == "(none)") {
        optMetric <- TRUE
        
      } else {
        # Edit cov_available and cov_used
        cov_used <- union(cov_used, cov_retain)
        cov_available <- setdiff(cov_available, cov_used)
        
        allCovUsed <- length(cov_available) == 0
        
        # Edit the metrics to compare
        rmseBase <- cov_metrics[1, "RMSE"]
        r2Base <- cov_metrics[1, "R2"]
        
        # Edit the model matrix
        Xbase <- cbind(Xbase, Xall[,cov_retain, drop = FALSE])
        
      }
      
      step_out[[i]] <- cov_metrics
      i <- i + 1
      
    } # End the while loop
    
    # Return a list
    out <- list(optVariables = cov_used, finalResults = step_out[[length(step_out)]][1,], stepTestResults = step_out)
    
    
  } else if (direction == "backward") {
    
    
    # Split based on multivariate-ness
    if (!is.multivariate) {
      
      ## Find the RMSE of the model with all variables
      base_fit <- plsr(yall ~ Xall, validation = "LOO", scale = FALSE, center = FALSE)
      rmse <- rmseBase <- RMSEP(object = base_fit)$val[,,comp+1]["CV"]
      r2 <- r2Base <- R2(base_fit)$val[,,comp+1]
      
    } else {
      
      # Perform the resampling
      base_resample_pred <- lapply(X = index, FUN = function(modelInd) {
        X1 <- Xall[modelInd,,drop = FALSE]
        y1 <- yall[modelInd,,drop = FALSE]
        xTest <- Xall[-modelInd,,drop = FALSE]
        yTest <- yall[-modelInd,,drop = FALSE]
        fit <- plsr(y1 ~ X1)
        # Predict and return
        t(rbind(yTest, predict(object = fit, newdata = xTest, comp = seq_len(comp))))
      })
      
      # Row bind
      out <- do.call("rbind", base_resample_pred)
      # Calculate the metrics
      res <- out[,1] - out[,2]
      rmse <- rmseBase <- sqrt(mean(res^2, na.rm = TRUE))
      r2 <- r2Base <- cor(out)[1,2]^2
      
    }
    
    
    ## Start the algorithm
    
    # Set Xbase
    Xbase <- Xall
    # Define the covariates that may be dropped
    cov_todrop <- colnames(Xbase)
    # Define the covariates that have been dropped
    cov_dropped <- character()
    # convert to named vector of column numbers
    col_todrop <- setNames(object = seq_len(ncol(Xbase)), cov_todrop)
    
    # Define a flag to stop the algorithm
    optMetric <- FALSE
    allCovUsed <- length(cov_todrop) == 0
    
    
    # A list to store output
    step_out <- list()
    i = 1
    
    # While loop
    while (!optMetric & !allCovUsed) {
      
      # Iterate over cov_todrop and build models
      cov_removal_pred <- sapply(X = col_todrop, FUN = function(col_drop) {
        
        Xcov <- Xbase[,-col_drop, drop = FALSE]
        
        # Split based on multivariate-ness
        if (!is.multivariate) {
          
          ## Find the RMSE of the model with all variables
          cov_fit <- plsr(yall ~ Xcov, validation = "LOO", scale = FALSE, center = FALSE)
          rmseCov <- RMSEP(object = cov_fit)$val[,,comp+1]["CV"]
          r2Cov <- R2(cov_fit)$val[,,comp+1]
          
        } else {
          
          # Resampling
          cov_resample_pred <- lapply(X = index, FUN = function(modelInd) {
            X1 <- Xcov[modelInd,,drop = FALSE]
            y1 <- yall[modelInd,,drop = FALSE]
            xTest <- Xcov[-modelInd,,drop = FALSE]
            yTest <- yall[-modelInd,,drop = FALSE]
            fit <- plsr(y1 ~ X1)
            # Predict and return
            t(rbind(yTest, predict(object = fit, newdata = xTest, comp = seq_len(min(comp, ncol(X1))))))
          })
          
          # Row bind
          out <- do.call("rbind", cov_resample_pred)
          # Calculate the metrics
          resCov <- out[,1] - out[,2]
          rmseCov <- sqrt(mean(resCov^2, na.rm = TRUE))
          r2Cov <- cor(out)[1,2]^2
          
        }
        
        # Return a matrix
        c(RMSE = unname(rmseCov), R2 = r2Cov)
        
      })
      
      # Analyze the metrics
      cov_removal_metrics <- t(cbind(cov_removal_pred, `(none)` = c(rmseBase, r2Base)))
      cov_removal_metrics <- cov_removal_metrics[order(cov_removal_metrics[,metric], decreasing = FALSE),, drop = FALSE]
      
      # Select the covariate to add
      cov_remove <- row.names(cov_removal_metrics)[1]
      
      # If "none", stop
      if (cov_remove == "(none)") {
        optMetric <- TRUE
        
      } else {
        # Edit cov_available and cov_used
        cov_dropped <- union(cov_dropped, cov_remove)
        cov_todrop <- setdiff(cov_todrop, cov_dropped)
        
        allCovUsed <- length(cov_todrop) == 0
        
        # Edit the metrics to compare
        rmseBase <- cov_removal_metrics[1, "RMSE"]
        r2Base <- cov_removal_metrics[1, "R2"]
        
        # Edit the model matrix
        Xbase <- Xbase[,colnames(Xbase) != cov_remove, drop = FALSE]
        col_todrop <- setNames(object = seq_len(ncol(Xbase)), cov_todrop)
        
      }
      
      step_out[[i]] <- cov_removal_metrics
      i <- i + 1
      
    } # End the while loop
    
    out <- list(optVariables = cov_todrop, finalResults = step_out[[length(step_out)]][1,], stepTestResults = step_out)
    
  } # End if statement
  
  # Return results
  return(out)
  
}



## A function for the recursive feature addition procedure
rfa_proc <- function(env.col = "environment") {
  
  covariates_use <- subset(trait_covariate_df, trait == unique(df$trait), covariate, drop = TRUE)
  
  ## Recursive feature addition ##
  ## Main effect
  # 1. Fit a base model
  base_fit <- lm(value ~ 1 + line_name, data = df1)
  # 2. Define the scope
  scope <- list(lower = formula(base_fit), upper = reformulate(c("line_name", covariates_use), response = "value"))
  # Run rfa
  rfa_out <- rfa_loo(object = base_fit, data = df1, scope = scope, metric = "RMSE", index = loo_indices, env.col = env.col)
  
  ## Interactions
  # 1. Fit a base model
  base_fit_int <- update(base_fit, formula = reformulate(rfa_out$optVariables, response = "value"))
  # 2. Define the scope
  scope <- list(lower = formula(base_fit_int), 
                upper = reformulate(c(rfa_out$optVariables, paste0("line_name:", covariates_use)), response = "value"))
  # Run rfa
  rfa_out_int <- rfa_loo(object = base_fit_int, data = df1, scope = scope, metric = "RMSE", index = loo_indices, env.col = env.col)
  
  
  ##  no soil
  covariates_use <- subset(trait_covariate_df, trait == unique(df$trait) & str_detect(covariate, "soil", negate = T), 
                           covariate, drop = TRUE) 
  ## Main effect
  # 1. Fit a base model
  base_fit <- lm(value ~ 1 + line_name, data = df1)
  # 2. Define the scope
  scope <- list(lower = formula(base_fit), upper = reformulate(c("line_name", covariates_use), response = "value"))
  # Run rfa
  rfa_out_nosoil <- rfa_loo(object = base_fit, data = df1, scope = scope, metric = "RMSE", index = loo_indices, env.col = env.col)
  
  ## Interactions
  # 1. Fit a base model
  base_fit_int <- update(base_fit, formula = reformulate(rfa_out_nosoil$optVariables, response = "value"))
  # 2. Define the scope
  scope <- list(lower = formula(base_fit_int), 
                upper = reformulate(c(rfa_out_nosoil$optVariables, paste0("line_name:", covariates_use)), response = "value"))
  # Run rfa
  rfa_out_int_nosoil <- rfa_loo(object = base_fit_int, data = df1, scope = scope, metric = "RMSE", index = loo_indices, env.col = env.col)
  
  
  ## Create a tibble
  tibble(
    feat_sel_type = "rfa_cv",
    model = c("model2", "model3"),
    adhoc = list(rfa_out, rfa_out_int),
    adhoc_nosoil = list(rfa_out_nosoil, rfa_out_int_nosoil)
  )
  
}



## A function for recursive feature elimination
rfe_proc <- function(main_fit, env.col = "environment") {
  
  df_env_eff <- coef(main_fit) %>% 
    tibble(term = names(.), effect = .) %>% 
    filter(str_detect(term, env.col)) %>% 
    mutate(term = str_remove(term, env.col)) %>% 
    add_row(term = last(levels(df1[[env.col]])), effect = -sum(.$effect))
  
  # LOO indices
  loo_indices <- df_env_eff %>%
    group_by(term) %>%
    crossv_loo_grouped() %>%
    pull(train) %>%
    map("idx")
  
  # run RFE
  Y <- as.data.frame(df_env_eff) %>% 
    column_to_rownames("term") %>%
    as.matrix()
  
  ## RFS - adhoc
  covariates_use <- subset(trait_covariate_df, trait == unique(df$trait), covariate, drop = TRUE)
  X <- distinct(select(df1, env.col, covariates_use)) %>% 
    as.data.frame(df_env_eff) %>% 
    column_to_rownames(env.col) %>% 
    as.matrix()
  # Execute
  rfs_main_adhoc <- rfs_loo(X = X, Y = Y, metric = "RMSE", direction = "forward", index = loo_indices)
  rfs_main_adhoc_back <- rfs_loo(X = X, Y = Y, metric = "RMSE", direction = "backward", index = loo_indices)
  
  ## RFS - adhoc no soil
  covariates_use <- subset(trait_covariate_df, trait == unique(df$trait) & str_detect(covariate, "soil", negate = T), 
                           covariate, drop = TRUE)    
  X <- distinct(select(df1, env.col, covariates_use)) %>% 
    as.data.frame(df_env_eff) %>% 
    column_to_rownames(env.col) %>% 
    as.matrix()
  
  # Run RFS
  rfs_main_adhoc_nosoil <- rfs_loo(X = X, Y = Y, metric = "RMSE", direction = "forward", index = loo_indices)
  rfs_main_adhoc_nosoil_back <- rfs_loo(X = X, Y = Y, metric = "RMSE", direction = "backward", index = loo_indices)
  
  
  ## Model interactions ##
  
  ## Return a df of residuals (GxE)
  df_gxe_eff <- add_residuals(df1, main_fit, var = "int_effect") %>%
    select(line_name, term = env.col, covariates_use, int_effect)
  
  # LOO indices
  loo_indices <- df_env_eff %>%
    group_by(term) %>%
    crossv_loo_grouped() %>%
    pull(train) %>%
    map("idx")
  
  # run RFE
  Y <- df_gxe_eff %>% select(line_name, term, int_effect) %>% 
    spread(term, int_effect) %>% 
    as.data.frame() %>% 
    column_to_rownames("line_name") %>%
    as.matrix() %>%
    t()
  
  # Adhoc
  covariates_use <- subset(trait_covariate_df, trait == unique(df$trait), covariate, drop = TRUE)
  X <- distinct(select(df1, env.col, covariates_use)) %>% 
    as.data.frame(df_env_eff) %>% 
    column_to_rownames(env.col) %>% 
    as.matrix()
  # Execute
  rfs_int_adhoc <- rfs_loo(X = X, Y = Y, metric = "RMSE", direction = "forward", index = loo_indices)
  rfs_int_adhoc_back <- rfs_loo(X = X, Y = Y, metric = "RMSE", direction = "back", index = loo_indices)
  
  # Adhoc - no soil
  covariates_use <- subset(trait_covariate_df, trait == unique(df$trait) & str_detect(covariate, "soil", negate = T), 
                           covariate, drop = TRUE)  
  X <- distinct(select(df1, env.col, covariates_use)) %>% 
    as.data.frame(df_env_eff) %>% 
    column_to_rownames(env.col) %>% 
    as.matrix()
  # Execute
  rfs_int_adhoc_nosoil <- rfs_loo(X = X, Y = Y, metric = "RMSE", direction = "forward", index = loo_indices)
  rfs_int_adhoc_nosoil_back <- rfs_loo(X = X, Y = Y, metric = "RMSE", direction = "backward", index = loo_indices)
  
  
  ## Return results
  tibble(
    feat_sel_type = "rfs_pls",
    direction = rep(c("forward", "backward"), each = 2),
    model = rep(c("model2", "model3"), 2),
    adhoc = list(rfs_main_adhoc, rfs_int_adhoc, rfs_main_adhoc_back, rfs_int_adhoc_back),
    adhoc_nosoil = list(rfs_main_adhoc_nosoil, rfs_int_adhoc_nosoil, rfs_main_adhoc_nosoil_back, rfs_int_adhoc_nosoil_back))
  
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
  


## Function for calculating a distance matrix based on two-way matrix
make_dist_mat <- function(x) {
  d1 <- as.matrix(dist(x))
  1 - (d1 / max(d1))
}



    

## The main prediction function for CV ##
# 
# x is a row with training and testing df
# 
# Other objects should be in the global environment:
# ec_model_building
# ec_model_scaled
# K
# model_fixed_forms
# model_rand_forms
# 
# 
# 
# 
genomewide_prediction <- function(x) {
  
  # Get training and test data
  row <- x
  tr <- unique(row$trait)
  
  # train <- droplevels(subset(data_to_model, id %in% row$train[[1]]$id)) %>%
  #   mutate(line_name = factor(line_name, levels = c(tp_geno, vp_geno)))
  # test <- subset(data_to_model, id %in% row$test[[1]]$id) %>%
  #   select(environment, line_name, value)
  
  ## Different subsetting method
  train <- row$train[[1]]
  test <- row$test[[1]]
  
  # Record the number of environment and observations used for training
  train_n <- summarize(train, nEnv = n_distinct(environment), nObs = n())
  
  
  # ## Get the covariate model and extract the formula
  # ec_model_form <- ec_model_building %>% 
  #   unnest(final_model) %>%
  #   subset(trait == row$trait & model == "model3_ammi", object, drop = TRUE) %>%
  #   # It's a list, so take the first element
  #   first() %>%
  #   formula()
  # 
  # # main_environment_covariates
  # main_environment_covariates <- all.vars(ec_model_form) %>% 
  #   subset(., map_lgl(., ~any(str_detect(string = str_subset(string = attr(terms(ec_model_form), "term.labels"), pattern = "\\|", negate = T), 
  #                                        pattern = .))))
  # 
  # # interaction_environment_covariates
  # interaction_environment_covariates <- all.vars(ec_model_form) %>% 
  #   subset(., map_lgl(., ~any(str_detect(string = str_subset(string = attr(terms(ec_model_form), "term.labels"), pattern = "\\|"), 
  #                                        pattern = .)))) %>% setdiff(., "line_name")
  # 
  # ## Subset covariates into a df
  # ec_df <- ec_tomodel_scaled %>%
  #   # Select environment and the relevant covariates
  #   select(env = environment, union(main_environment_covariates, interaction_environment_covariates)) %>%
  #   filter(env %in% levels(train$env)) %>%
  #   as.data.frame() %>%
  #   column_to_rownames("env")
  # 
  
  
  # ## Create a list of model formulas
  # model_fixed_forms <- formulas(
  #   .response = ~ value, 
  #   model1 = ~ 1,
  #   model2 = model1,
  #   # model2a = add_predictors(model2, ~ env),
  #   # model2b = reformulate(c("1", main_environment_covariates)),
  #   model3 = model2,
  #   # model3a = model2a,
  #   # model3b = model2b
  #   model4 = model1,
  #   model5 = model4,
  #   model2_P = model2,
  #   model3_P = model3,
  #   model4_P = model4,
  #   model5_P = model5,
  #   model3_GxE = model3,
  #   model5_GxL = model5
  # )
  
  
  ## Models for de novo fitting 
  ## Create a list of model formulas
  model_rand_forms <- formulas(
    .response = ~ value,
    model1 = ~ vs(line_name, Gu = K),
    model2 = add_predictors(model1, ~ vs(env, Gu = E)),
    # model2a = model1,
    # model2b = model1,
    model3 = add_predictors(model2, ~ vs(line_name:env1, Gu = GE)),
    # model3a = add_predictors(model1, ~ vs(line_name:env1, Gu = GE)),
    # model3b = model3a
    model4 = add_predictors(model1, ~vs(loc, Gu = L)),
    model5 = add_predictors(model4, ~vs(line_name:loc, Gu = GL)),
    # Modified models
    model2_P = add_predictors(model1, ~ vs(env, Gu = E_IPC)),
    model3_P = add_predictors(model2, ~ vs(line_name:env1, Gu = GE_IPC)),
    model4_P = add_predictors(model1, ~vs(loc, Gu = L_IPC)),
    model5_P = add_predictors(model4, ~vs(line_name:loc, Gu = GL_IPC)),
    model3_GxE = add_predictors(model2, ~ vs(line_name:env1, Gu = GxE)),
    model5_GxL = add_predictors(model4, ~vs(line_name:loc, Gu = GxL))
  ) %>% map(~ formula(delete.response(terms(.)))) # Remove response
  
  # Residual formula
  resid_form <- ~ vs(units)
  
  
  
  ## Create relationship matrices
  K <- K # Genomic
  # E <- Env_mat(x = ec_df[,main_environment_covariates, drop = FALSE], method = "Rincent2019")
  E <- subset(environmental_relmat_df, trait == tr & term == "main", EC, drop = TRUE)[[1]]
  L <- subset(location_relmat_df, trait == tr & time_frame == time_frame_use & term == "main", EC, drop = TRUE)[[1]]
  
  GE <- subset(environmental_relmat_df, trait == tr & term == "int", EC, drop = TRUE)[[1]] %>%
    kronecker(X = K, Y = ., make.dimnames = TRUE)
  GL <- subset(location_relmat_df, trait == tr & time_frame == time_frame_use & term == "int", EC, drop = TRUE)[[1]] %>%
    kronecker(X = K, Y = ., make.dimnames = TRUE)
  
  GxE <- kronecker(X = K, Y = E, make.dimnames = TRUE)
  GxL <- kronecker(X = K, Y = L, make.dimnames = TRUE)
  
  
  ## AMMI covariance matrices
  E_IPC <- subset(environmental_relmat_df, trait == tr & term == "main", K, drop = TRUE)[[1]]
  L_IPC <- subset(location_relmat_df, trait == tr & time_frame == time_frame_use & term == "int", K, drop = TRUE)[[1]]
  
  GE_IPC <- subset(environmental_relmat_df, trait == tr & term == "int", K, drop = TRUE)[[1]] %>%
    kronecker(X = K, Y = ., make.dimnames = TRUE)
  GL_IPC <- subset(location_relmat_df, trait == tr & time_frame == time_frame_use & term == "int", K, drop = TRUE)[[1]] %>%
    kronecker(X = K, Y = ., make.dimnames = TRUE)
  
  ## Define an expression that fits the model
  fit_mmer_exp <- expression({
    model_fit <- mmer(fixed = fixed_form, random = rand_form, rcov = resid_form,
                      data = train, date.warning = FALSE, verbose = TRUE)
  })
  
  
  
  #################
  ## Fit models and extract predictions
  #################
  
  prediction_out <- tibble(trait = tr, model = names(model_rand_forms)) %>%
    mutate(prediction = list(NULL))
  
  # Test df to merge
  test_merge <- select(test, line_name, env, loc = location, year, value) %>%
    mutate_if(is.factor, as.character)
  
  # Iterate over models
  for (m in seq(nrow(prediction_out))) {
    
    # Model name and formulas
    mod <- prediction_out$model[m]
    # fixed_form <- model_fixed_forms[[mod]]
    fixed_form <- value ~ 1
    rand_form <- model_rand_forms[[mod]]
    
    # ## Get the variance components from the full model
    # full_varcomp <- subset(prediction_out, model == mod, sigma, drop = T)[[1]]
    # # Assign values to separate objects
    # varG <- full_varcomp$`u:line_name`
    # varE <- full_varcomp$`u:env`
    # varGE <- full_varcomp$`u:;line_name:env`
    # varR <- full_varcomp$units
    
    # Use rrBLUP to fit model1
    if (mod == "model1") {
      
      model_fit <- kin.blup(data = as.data.frame(train), geno = "line_name", pheno = "value", K = K)
      
      # Fixed effects
      fixed_eff <- matrix(data = (model_fit$pred - model_fit$g)[1], nrow = 1, ncol = 1,
                          dimnames = list("(Intercept)", "estimate"))
      
      ## Random effects
      rand_eff <- list(`u:line_name` = c(model_fit$g))
      
      
    
    } else {
      
      ## Try to fit the model; capture the output
      model_stdout <- capture.output( eval(fit_mmer_exp) )
      
      # If model fit is empty, try using a smaller number of iterations; for instance find
      # the maximum logLik and use those iterations
      itry <- 1
      while (is_empty(model_fit) & itry == 1) {
        
        # Find the number of iterations that maximized the logLik
        best_iter <- model_stdout %>% 
          subset(., str_detect(., "singular", negate = T)) %>% 
          read_table(.) %>%
          subset(., LogLik == max(LogLik), iteration, drop = TRUE)
        
        # Refit
        eval(fit_mmer_exp)
        
        # Increase the counter
        itry = itry + 1
        
      }
      
      # If the model is still empty, create empty fixed and random effects
      if (is_empty(model_fit)) {
        fixed_eff <- matrix(as.numeric(NA), nrow = 1, ncol = 1, dimnames = list("(Intercept)", "estimate"))
        
        rand_eff <- list("u:line_name" = set_names(x = rep(NA, nlevels(test$line_name)), nm = levels(test$line_name)))
        
      } else {
        
        ## Fixed effects
        fixed_eff <- coef(model_fit) %>%
          select(term = Effect, estimate = Estimate) %>%
          column_to_rownames("term") %>%
          as.matrix()
        
        ## Random effects
        rand_eff <- map(randef(model_fit), "value")
        
      }

    }
      
      
    ## Vector of new column names for separation, if necessary
    separation_col_names <- str_extract(string = attr(terms(rand_form), "term.labels"), pattern = "[a-z_]{1,}:[a-z]{1,}") %>% 
      str_subset(., ":") %>% 
      str_split(., ":") %>% 
      unlist()
      
    ## Create an X matrix for test
    Xtest <- model.matrix(fixed_form, test)[,row.names(fixed_eff), drop = FALSE] # This seems to work
    # Calculate fixed effects by the formula Xb
    fixed_pred <- Xtest %*% fixed_eff
    
    
    ## Convert to a complete data.frame
    rand_eff_df <- rand_eff %>% 
      map(~tibble(term = names(.), pred = .x)) %>% 
      imap(~`names<-`(.x, c(str_remove_all(.y, "u:|[0-9]"), paste0("pred", .y)))) %>% 
      modify_if(~str_detect(names(.x)[1], ":"), 
                ~separate(.x, col = 1, into = separation_col_names, sep = ":")) %>% 
      .[order(map_dbl(., ncol), decreasing = T)] %>% 
      map(~left_join(test_merge, .x)) %>%
      reduce(full_join) %>% 
      mutate(pred_incomplete = rowSums(select(., contains("pred")))) %>% 
      select(-contains(":"))
    
    # Add predictions to the list
    prediction_out$prediction[[m]] <- mutate(rand_eff_df, pred_complete = pred_incomplete + c(fixed_pred))
    
  }
  
  # Return predictions
  return(list(prediction_out = prediction_out, train_n = train_n))
  
  
}


## Different version of prediction function
genomewide_prediction2 <- function(x, model.list, K, E, KE) {
  
  # Get training and test data
  row <- x
  tr <- unique(row$trait)
  
  
  # Edit train and test by adding covariates - centered only
  train <- row$train[[1]]
  test <- row$test[[1]]
  
  # Extract fixed and random formulas from the list
  model_fixed_forms <- model.list$fixed
  model_rand_forms <- model.list$random
  

  # Residual formula
  resid_form <- ~ vs(units)
  

  
  
  ## Create relationship matrices
  K <- K # Genomic
  E <- E
  KE <- KE 
  # Identity matrix for environments
  I <- diag(ncol(E)); dimnames(I) <- dimnames(E)
  GE <- kronecker(X = K, Y = KE, make.dimnames = TRUE)
  GxE <- kronecker(X = K, Y = E, make.dimnames = TRUE) # This is crossproduct of K and E covariates
  GI <- kronecker(X = K, Y = I, make.dimnames = TRUE)
  
  ## Define an expression that fits the model
  fit_mmer_exp <- expression({
    model_fit <- mmer(fixed = fixed_form, random = rand_form, rcov = resid_form,
                      data = train, date.warning = FALSE, verbose = TRUE, getPEV = TRUE)
  })
  
  
  
  #################
  ## Fit models and extract predictions
  #################
  
  prediction_out <- tibble(trait = tr, model = names(model_rand_forms)) %>%
    mutate(prediction = list(NULL))
  
  # Test df to merge
  test_merge <- test %>%
    select(which(names(.) %in% c("line_name", "env", "loc", "site", "year", "value"))) %>%
    mutate_if(is.factor, as.character)
  
  # Iterate over models
  for (m in seq(nrow(prediction_out))) {
    
    # Model name and formulas
    mod <- prediction_out$model[m]
    # fixed_form <- model_fixed_forms[[mod]]
    fixed_form <- model_fixed_forms[[mod]]
    rand_form <- model_rand_forms[[mod]]
    
    # ## Get the variance components from the full model
    # full_varcomp <- subset(prediction_out, model == mod, sigma, drop = T)[[1]]
    # # Assign values to separate objects
    # varG <- full_varcomp$`u:line_name`
    # varE <- full_varcomp$`u:env`
    # varGE <- full_varcomp$`u:;line_name:env`
    # varR <- full_varcomp$units
    
    # Use rrBLUP to fit model1 or model2_fr
    if (mod == "model1") {
      
      mf <- model.frame(value ~ line_name, train)
      y <- model.response(mf)
      Z <- model.matrix(~ -1 + line_name, mf)
      
      model_fit <- mixed.solve(y = y, Z = Z, K = K, SE = TRUE)
      
      # Fixed effects
      fixed_eff <- cbind(estimate = model_fit$beta, std_error = model_fit$beta.SE)
      row.names(fixed_eff) <- "(Intercept)"
      
      ## Random effects
      rand_eff <- list(`u:line_name` = cbind(u = model_fit$u, pev = model_fit$u.SE^2))
      
    } else if (mod == "model2_fr") {
      
      mf <- model.frame(formula = reformulate(c("line_name", covariate_list$main), response = "value"), data = train)
      y <- model.response(mf)
      # Covariates
      X <- model.matrix(reformulate(termlabels = covariate_list$main), data = mf)
      # If X is not full rank, drop the last term
      while (qr(crossprod(X))$rank < ncol(X)) X <- X[,-ncol(X)]
      
      Z <- model.matrix(~ -1 + line_name, data = mf)
      
      model_fit <- mixed.solve(y = y, Z = Z, K = K, X = X, SE = TRUE)
      
      ## Matrices of fixed and random effects
      fixed_eff <- as.matrix(model_fit$beta)

      ## Random effects
      rand_eff <- list(`u:line_name` = c(model_fit$u))
      
      
      # All other models use SOMMER
    } else {
      
      ## Try to fit the model; capture the output
      model_stdout <- capture.output( eval(fit_mmer_exp) )
      
      
      # # If model fit is empty, try using a smaller number of iterations; for instance find
      # # the maximum logLik and use those iterations
      # itry <- 1
      # while (is_empty(model_fit) & itry == 1) {
      #   
      #   # Find the number of iterations that maximized the logLik
      #   best_iter <- model_stdout %>% 
      #     subset(., str_detect(., "singular", negate = T)) %>% 
      #     read_table(.) %>%
      #     subset(., LogLik == max(LogLik), iteration, drop = TRUE)
      #   
      #   # Refit
      #   eval(fit_mmer_exp)
      #   
      #   # Increase the counter
      #   itry = itry + 1
      #   
      # }
      
      # If the model is still empty, create empty fixed and random effects
      if (is_empty(model_fit)) {
        fixed_eff <- matrix(as.numeric(NA), nrow = 1, ncol = 2, dimnames = list("(Intercept)", c("estimate", "std_error")))
        
        rand_eff <- list("u:line_name" = cbind(u = set_names(x = rep(NA, nlevels(test$line_name)), nm = levels(test$line_name)),
                                               pev = NA))
        
      } else {
        
        # Check to make sure the loglikihood was maximized
        maxed_LL <- which.max(model_fit$monitor[1,]) == ncol(model_fit$monitor)
        
        if (!maxed_LL) {
          # Refit the model using the iterations that maximized the LL
          model_fit <- mmer(fixed = fixed_form, random = rand_form, rcov = resid_form,
                            data = train, date.warning = FALSE, verbose = TRUE, getPEV = TRUE,
                            iters = which.max(model_fit$monitor[1,]))
        }
        
        ## Fixed effects
        fixed_eff <- coef(model_fit) %>%
          select(term = Effect, estimate = Estimate) %>%
          mutate(std_error = sqrt(c(model_fit$VarBeta))) %>%
          column_to_rownames("term") %>%
          as.matrix()
        
        ## Random effects
        rand_eff <- map(randef(model_fit), "value")
        # PEV
        pev <- map(model_fit$PevU, "value") %>% map(diag)
        # combine
        rand_eff <- map2(rand_eff, pev, ~cbind(u = .x, pev = .y))
        
      }
      
    }
    
    
    
    
    ## Vector of new column names for separation, if necessary
    separation_col_names <- str_extract(string = attr(terms(rand_form), "term.labels"), pattern = "[a-z_]{1,}:[a-z]{1,}") %>% 
      str_subset(., ":") %>% 
      str_split(., ":") %>% 
      unlist()
    
    ## Create an X matrix for test
    Xtest <- model.matrix(fixed_form, test)[,row.names(fixed_eff), drop = FALSE] # This seems to work
    # Calculate fixed effects by the formula Xb
    fixed_pred <- Xtest %*% fixed_eff[,1]
    fixed_pred_se <-  Xtest %*% fixed_eff[,2] ## This will not work with quantitative elements of X
    
    ## If model3_fr, calculate random effects
    if (mod == "model3_fr") {
      
      ## Create an X matrix for test
      Ztest <- test %>%
        mutate(intercept = 1) %>%
        select(environment, intercept, covariate_list$interaction) %>%
        distinct() %>%
        as.data.frame() %>%
        column_to_rownames("environment") %>%
        as.matrix() %>%
        t()
      
      ## Calculate random effects - main effect and environmental interaction
      rand_eff_hat <- do.call("cbind", rand_eff) %*% Ztest
      
      # Gather
      rand_eff_df <- rand_eff_hat %>%
        as.data.frame() %>% 
        rownames_to_column("line_name") %>% 
        gather(env, pred_incomplete, -line_name) %>%
        left_join(test_merge, .)
        
      
    
    } else {
      
      # ## Convert to a complete data.frame
      # rand_eff_df <- rand_eff %>% 
      #   map(~tibble(term = names(.), pred = .x)) %>% 
      #   imap(~`names<-`(.x, c(str_remove_all(.y, "u:|[0-9]"), paste0("pred", .y)))) %>% 
      #   modify_if(~str_detect(names(.x)[1], ":"), 
      #             ~separate(.x, col = 1, into = separation_col_names, sep = ":")) %>% 
      #   .[order(map_dbl(., ncol), decreasing = T)] %>% 
      #   map(~left_join(test_merge, .x)) %>%
      #   reduce(full_join) %>% 
      #   mutate(pred_incomplete = rowSums(select(., contains("pred")))) %>% 
      #   select(-contains(":"))
      
      ## Convert to a complete data.frame
      rand_eff_df <- rand_eff %>% 
        map(~rownames_to_column(as.data.frame(.), "term")) %>%
        map(~rename(., pred = u)) %>%
        imap(~`names<-`(.x, c(str_remove_all(.y, "u:|[0-9]"), paste0("pred", .y), paste0("pev", .y)))) %>% 
        modify_if(~str_detect(names(.x)[1], ":"), 
                  ~separate(.x, col = 1, into = separation_col_names, sep = ":")) %>% 
        .[order(map_dbl(., ncol), decreasing = T)] %>% 
        map(~left_join(test_merge, .x)) %>%
        reduce(full_join) %>% 
        mutate(pred_incomplete = rowSums(select(., contains("pred"))), pred_incomplete_pev = rowSums(select(., contains("pev")))) %>% 
        select(-contains(":"))
      
    }
    
    

    
    # Add predictions to the list
    prediction_out$prediction[[m]] <- mutate(rand_eff_df, pred_complete = pred_incomplete + c(fixed_pred),
                                             pred_complete_pev = pred_incomplete_pev + c(fixed_pred_se^2))
    
  }
  
  # Return predictions
  return(list(prediction_out = prediction_out))
  
  
}










# function to impute missing data using mean
impute <- function(x) {
  x1 <- x
  x1[is.na(x1)] <- mean(x1, na.rm = T)
  return(x1)
}

# Function to calculate the coefficient of variation
cv <- function(x, na.rm = FALSE) sd(x = x, na.rm = na.rm) / mean(x = x, na.rm = na.rm)


# First define a function to calculate bias
bias <- function(obs, pred) mean((pred - obs) / obs)


