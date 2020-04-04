
#' Sobol sequence
#'
#' I don't like the randtoolbox because you can't have two
#' separate concurrent Sobol sequences.
sobol_design <- R6::R6Class(
  "sobol_design",
  public = list(
    X = NULL,
    D = NULL,
    L = NULL,
    b = NULL,
    seed = NULL,
    use_lhs = NULL,
    initialize = function(D, L, seed=numeric(0)) {
      self$D <- D
      self$L <- L
      self$seed <- seed
      # D=1 sobol gives vector, as.matrix makes it matrix with 1 col
      self$X <- as.matrix(randtoolbox::sobol(n=L,dim=D, init=T, seed=seed, scrambling=1))
    },
    get.batch = function(L=self$L) {
      # if (length(self$seed) > 0) {
      #   set.seed(self$seed)
      #   self$seed <- self$seed + 1
      # }
      # D=1 sobol gives vector, as.matrix makes it matrix with 1 col
      newX <- as.matrix(randtoolbox::sobol(n=L, dim=self$D, init=F, scrambling=1))
      self$X <- rbind(self$X, newX)
      newX
    }
  )
)

