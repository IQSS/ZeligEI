#' Wakefield's Hierarchical Ecological Inference Model
#'
#' Vignette: \url{http://docs.zeligproject.org/articles/zeligei_eihier.html}
#' @import methods
#' @export Zelig-eihier
#' @exportClass Zelig-eihier
#'
#' @include model-ei.R

zeihier <- setRefClass("Zelig-eihier",
                          contains = "Zelig-ei")


zeihier$methods(
  initialize = function() {
    callSuper()
    .self$name <- "eihier"
    .self$description <- "Wakefield's Hierarchical Ecological Inference Model"
    .self$fn <- quote(MCMCpack::MCMChierEI)
    .self$packageauthors <- "Andrew D. Martin, Kevin M. Quinn, Jong Hee Park"
    .self$wrapper <- "eihier"
    .self$vignette.url <- "http://docs.zeligproject.org/articles/zeligei_eihier.html"
    ref1<-bibentry(
            bibtype="Article",
            title = "Ecological Inference for 2 x 2 Tables.",
            author = person("Jonathan C.", "Wakefield"),
            journal = "Journal of the Royal Statistical Society, Series A.",
            volume = 167,
            number = 3,
            year = 2004,
            pages = "385--445")
    .self$refs<-c(.self$refs,ref1)
  }
)

zeihier$methods(
  zelig = function(formula, data, N=NULL, ..., weights = NULL, by = NULL, bootstrap = FALSE, na.action="na.omit") {
    na.action <- checkZeligEIna.action(na.action)

    if(!identical(bootstrap,FALSE)){
      stop("Error: The bootstrap is not available for Markov chain Monte Carlo (MCMC) models.")
    }
    if(!is.null(weights)){
      stop("Error: Weights are not implemented for the Wakefield Hierarchical EI model.  Try the eiml model if weights are required.")
    }
    if(!is.null(by)){
      stop("Error: The `by' argument is not implemented for the Wakefield Hierarchical EI model.  Try the eiml model if this is required,
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

zeihier$methods(
  param = function(z.out) {
    return(z.out)
  }
)

zeihier$methods(
  getcoef = function() {
    "Get estimated model coefficients"
    return(.self$zelig.out$z.out[[1]])
  } 
)

