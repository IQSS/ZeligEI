#' Multinomial Dirichlet model for Ecological Inference in RxC tables
#'
#' Vignette: \url{http://docs.zeligproject.org/articles/zeligei_eirxc.html}
#' @import methods
#' @export Zelig-eirxc
#' @exportClass Zelig-eirxc
#'
#' @include model-ei.R

zeirxc <- setRefClass("Zelig-eirxc",
                          contains = "Zelig-ei")

zeirxc$methods(
  initialize = function() {
    callSuper()
    .self$name <- "eirxc"
    .self$description <- "Multinomial Dirichlet model for Ecological Inference in RxC tables"
    .self$fn <- quote(eiPack::ei.MD.bayes)
    .self$packageauthors <- "Michael Kellerman, Olivia Lau"
    .self$wrapper <- "eirxc"
    .self$vignette.url <- "http://docs.zeligproject.org/articles/zeligei_eirxc.html"
    ref1<-bibentry(
            bibtype="Article",
            title = "Bayesian and Frequentist Inference for Ecological Inference: The R x C case.",
            author = c(
                person("Ori", "Rosen"),
                person("Wenxin", "Jiang"),
                person("Gary", "King"),
                person("Martin A.", "Tanner")
                ),
            journal = "Statistica Neerlandia",
            volume = 167,
            year = 2001,
            pages = "134--156")
    .self$refs<-c(.self$refs,ref1)
  }
)

zeirxc$methods(
  zelig = function(formula, data, N = NULL, ..., weights = NULL, by = NULL, bootstrap = FALSE, na.action="na.omit") {
    na.action <- checkZeligEIna.action(na.action)

    if(!identical(bootstrap,FALSE)){
      stop("Error: The bootstrap is not available for Markov chain Monte Carlo (MCMC) models.")
    }
    
    cnvt <- convertEIformula2(formula=formula, data=data, N=N, na.action=na.action, rxc=TRUE)
    localformula <- cnvt$formula
    localdata <- cnvt$data

    .self$zelig.call <- match.call(expand.dots = TRUE)
    .self$model.call <- match.call(expand.dots = TRUE)
    
    .self$model.call$N <- NULL
    .self$model.call$na.action <- NULL
    .self$model.call$total <- cnvt$totalName

    callSuper(formula = localformula, data = localdata, N=NULL, ..., weights = weights, by = by, bootstrap = bootstrap)
  }
)

zeirxc$methods(
  getcoef = function() {
    "Get estimated model coefficients"
    coeflist <- list()
    for(i in 1:length(.self$zelig.out$z.out)){
      coeflist[[i]]<-.self$zelig.out$z.out[[i]]$draw
    }
    return(coeflist) 
  } 
)
