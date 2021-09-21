#' Immune response prediction
#'
#' \code{predict_immune_response} predicts immune response using two algorithms: multi-task elastic net
#' and bayesian efficient multi-kernel algorithm. While BEMKL can exploit information across different input and output datasets,
#' multi-task elastic net can only do so for response variables. Another advantage of BEMKL is missing data handling, which is not
#' the case for the other algorithm.
#'
#' These algorithms use model parameters learned during training on different types of data in
#' order to compute the immune response.
#'
#' @importFrom utils combn
#' @importFrom stats na.omit
#'
#' @export
#'
#' @param pathways numeric matrix with data
#' @param immunecells numeric matrix with data
#' @param tfs numeric matrix with data
#' @param lrpairs numeric matrix with data
#' @param ccpairs numeric matrix with data
#' @param cancertype string character
#'
#' @return Predictions for each model building.
#'
#------------------------------------------------------------------------------------------------

predict_immune_response <- function(pathways=NULL, immunecells=NULL, tfs=NULL, lrpairs=NULL, ccpairs=NULL, cancertype){

  if(missing(cancertype)) stop("cancer type needs to be specified")
  if(all(is.null(pathways),is.null(immunecells), is.null(tfs), is.null(lrpairs), is.null(ccpairs))) stop("none signature specified")

  # Initialize variables
  views <- c(
    pathways = "gaussian",
    immunecells = "gaussian",
    tfs = "gaussian",
    lrpairs = "gaussian",
    ccpairs = "gaussian"
  )

  algorithm <-  c("Multi_Task_EN") #,"BEMKL")

  # Check which views are missing
  miss_views <- c(ifelse(missing(pathways), NA, 1),
                  ifelse(missing(immunecells), NA, 2),
                  ifelse(missing(tfs), NA, 3),
                  ifelse(missing(lrpairs), NA, 4),
                  ifelse(missing(ccpairs), NA, 5))

  # Single views
  view_simples <- lapply(miss_views[!is.na(miss_views)], function(X) {
    tmp <- views[X]
    return(tmp)
  })

  # All corresponding views
  view_combinations <- view_simples

  all_predictions <- lapply(1:length(view_combinations), function(X){

    view_info <- view_combinations[[X]]
    view_name <- paste(names(view_info), collapse="_")
    view_data <- lapply(tolower(names(view_info)), function(x) as.data.frame(get(x)))
    names(view_data) <- names(view_info)
    message(X,".view source: ", view_name, "\n")

    # Predict immune response using model parameters
    summary_alg <- lapply(algorithm, function(alg){

      if (alg %in% c("BEMKL")){

        pred_alg <- predict_with_bemkl(view_name = view_name,
                                       view_info = view_info,
                                       view_data = view_data,
                                       learned_model = trained_models[[cancertype]][[view_name]])

      }else if (alg %in% c("Multi_Task_EN")){

        pred_alg <- predict_with_multitaskelasticnet(view_name = view_name,
                                                     view_info = view_info,
                                                     view_data = view_data,
                                                     learned_model = trained_models[[cancertype]][[view_name]])
      }
      return(pred_alg)
    })
    names(summary_alg) <- algorithm
    return(summary_alg)
  })
  names(all_predictions) <- names(unlist(view_combinations))
  return(all_predictions)

}
