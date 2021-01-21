#' @title Hutcheson test comparisons between multiple communities.

#' @description This function computes the p-values of the Hutcheson (1970) to test the significance of the difference between community diversity indexes.

#' @details This function performs Hutcheson's t-tests for comparing multiple diversity indexes pairwise. This test is based on Shannon diversity indices' value computed using a logarithm base specified by the user. The alternative hypothesis is one-sided, chosen automatically according to the sign of the difference between each pair of communities tested.

#' @note missing values will be replaced with zero values

#' @param x Numeric dataframe or matrix of the abundance of species per community. Each column corresponds to the communities and each row to a species.

#' @param shannon.base A numeric indicating the logarithm base for the Shanon indices. Defaults to \emph{e}.

#' @return A matrix whose entries are the p-values of the test. The names of the rows and columns are the names of the communities and its Shannon diversity index.

#' @seealso See \code{\link[ecolTest]{Hutcheson.T.test}} in \pkg{ecolTest} package.

#' @author David Ramirez Delgado \email{tucorreo@correo.com}.

#' @author Hugo Salinas \email{hugosal@comunidad.unam.mx}.

#' @examples

#' # poner algun ejemlpo aqui tambien

#' @export

multiple.Hutcheson.T.test<-function(x,shannon.base=exp(1)){
  x <- drop(as.matrix(x))
  if (!is.numeric(x))
    stop("input data must be numeric")
  if (any(x < 0, na.rm = TRUE))
    stop("input data must be non-negative")
  if (any(dim(x) < 2))
    stop("input data must have at least two columns and rows")
  nx<-ny<-1:ncol(x)
  N<-apply(x,2,sum)
  H<-(N*log(N, shannon.base)-apply(x*log(x,shannon.base), 2,sum,na.rm = TRUE))/N
  names(ny)<-names(nx)<-paste(colnames(x),round(H,2),sep = " ")
  f<-Vectorize(function(X,Y){Hutcheson.T.test(x[,X],x[,Y],
                shannon.base = shannon.base, alternative = "auto")$p.value})
  p_values<-outer(nx, ny,FUN = f)
  return(list(p.values = p_values))}
