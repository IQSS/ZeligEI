#' Quinn's Dynamic Ecological Inference Model
#'
#' Vignette: \url{http://docs.zeligproject.org/articles/zeligei_eidynamic.html}
#' @import methods
#' @export Zelig-eidynamic
#' @exportClass Zelig-eidynamic
#'
#' @include model-ei.R

zeidynamic <- setRefClass("Zelig-eidynamic",
                          contains = "Zelig-ei")


zeidynamic$methods(
  initialize = function() {
    callSuper()
    .self$name <- "eidynamic"
    .self$description <- "Quinn's Dynamic Ecological Inference Model"
    .self$fn <- quote(MCMCpack::MCMCdynamicEI)
    .self$packageauthors <- "Andrew D. Martin, Kevin M. Quinn, Jong Hee Park"
    .self$wrapper <- "eidynamic"
    .self$vignette.url <- "http://docs.zeligproject.org/articles/zeligei_eidynamic.html"
    ref1<-bibentry(
            bibtype="InCollection",
            title = "Ecological Inference in the Presence of Temporal Dependence.",
            booktitle = "Ecological Inference: New Methodological Strategies",
            author = person("Kevin", "Quinn"),
            year = 2004,
            publisher = "Cambridge University Press",
            organization = "Cambridge University Press",
            address = "New York",
            editor = c(person("Gary", "King"), person("Ori", "Rosen"), person("Martin", "Tanner")) 
            )
    .self$refs<-c(.self$refs,ref1)
  }
)

zeidynamic$methods(
  zelig = function(formula, data, N=NULL, ..., weights = NULL, by = NULL, bootstrap = FALSE, na.action="na.omit") {
    na.action <- checkZeligEIna.action(na.action)

    if(!identical(bootstrap,FALSE)){
      stop("Error: The bootstrap is not available for Markov chain Monte Carlo (MCMC) models.")
    }
    if(!is.null(weights)){
      stop("Error: This model is dynamic over time and currently Zelig does not have a weighting approach to work in this model.
        Check if you intended to use the W argument to adjust the temporal dependence among elements in the Quinn model.")
    }
    if(!is.null(by)){
      stop("Error: The `by' argument is not implemented for Quinn's Dynamic EI model.  Try the eiml model if this is required,
        or subset the data and run multiple models.")
    }

    cnvt <- convertEIformula(formula=formula, N=N, data=data, na.action=na.action)

    .self$zelig.call <- match.call(expand.dots = TRUE)
    .self$model.call <- match.call(expand.dots = TRUE)

    .self$model.call$r0 <- cnvt$r0
    .self$model.call$r1 <- cnvt$r1
    .self$model.call$c0 <- cnvt$c0
    .self$model.call$c1 <- cnvt$c1

    .self$model.call$N <- NULL       
    .self$model.call$formula <- NULL
    .self$model.call$data <- NULL
    .self$model.call$na.action <- NULL

    # Note, formula and data pass through the Zelig internals, but are ignored by the wrapped model
    callSuper(formula = formula, data = data, ..., weights = NULL, by = NULL, bootstrap = FALSE)
  }
)

zeidynamic$methods(
  param = function(z.out) {
    return(z.out)
  }
)

zeidynamic$methods(
  getcoef = function() {
    "Get estimated model coefficients"
    return(.self$zelig.out$z.out[[1]])
  } 
)
