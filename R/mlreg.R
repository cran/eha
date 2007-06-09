mlreg <-
function (formula = formula(data),
          data = parent.frame(), 
          na.action = getOption("na.action"),
          init = NULL,
          method = c("ML", "MPPL"),
          control = list(eps = 1e-8,
          maxiter = 10, n.points = 12, trace = FALSE),
          singular.ok = TRUE,
          model = FALSE, 
          x = FALSE,
          y = TRUE,
          boot = FALSE,
          geometric = FALSE,
          rs = NULL,
          frailty = NULL,
          max.survs = NULL) 
{
    warning("'mlreg' is deprecated; use 'coxreg' instead (see 'methods')")
    if (method[1] == "ML") method <- "ml"
    else if (method[1] == "MPPL") method <- "mppl"
    else stop(paste("Unknown method", as.character(method[1])))
    
    coxreg(formula,
          data, 
          na.action,
          init,
          method,
          control,
          singular.ok,
          model, 
          x,
          y,
          boot,
          geometric,
          rs,
          frailty,
          max.survs)
}
