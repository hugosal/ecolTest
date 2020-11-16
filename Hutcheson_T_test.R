#' @title Hutcheson T-test for two diversity indices.

#' 

#' @description This function performs Hutcheson (1970) to test the significance of the difference between two communities' diversity indices.

#' 

#' @details This function performs Hutcheson's t-test for comparing two diversity indices. This test is based on Shannon diversity indices' value computed using a logarithm base specified by the user. One-sided or two-sided tests are available.  

#' 

#' @note missing values will be replaced with zero values

#' 

#' @aliases Hutcheson.T.test.default

#' 

#' @param shanon.base A numeric indicating the logarithm base for the Shanon indices. Defaults to \emph{e}.

#' 

#' @param alternative A character indicating the alternative hypothesis. Can be "two.sided"(default), "less", "greater", or "auto"

#' 

#' @param difference  A numeric indicating the value hypothesized of the difference between the indices. Defaults to 0.

#' 

#' @return A list with class "htest" containing the following components:

#' 

#' \itemize{
#'   \item statistic: Value of the Hutcheson T statistic
#'   \item parameter: The degrees of freedom of the T statistic parameter
#'   \item p.value: The test's p-value.
#'   \item estimate: The Shannon diversity indices
#'   \item null.value:The  hypothesized value of the difference between the Shannon diversty indices
#'   \item method: Name of the test
#'   \item alternative: The alternative hypothesis
#'   \item data.name: Name of the data used in the test

#' }

#'

#' @seealso See \code{t.test}, \code{\link[t.test]{t.test}} in \pkg{stats}

#' 

#' @author Hugo Salinas \email{hugosal@comunidad.unam.mx}, this function is based on the Student's t-Test \code{\link[t.test]{t.test}} in \pkg{stats}.

#'

#' @references 

#' Zar, Jerrold H. 2010. Biostatistical Analysis. 5th ed. Pearson. pp. 174-176.

#' 

#' Hutcheson, Kermit. 1970. A Test for Comparing Diversities Based on the Shannon Formula. Journal of Theoretical Biology 29: 151-54.

#' 

#' @examples

#' # poner ejemplos, conseguir alguna base de datos que se publique con esto para usar aqui

#' # prueba una cola, etc

#' @export


Hutcheson.T.test<-function(x, y,
                           shanon.base = exp(1),
                            alternative = "two.sided", difference = 0){
  dname<-paste((deparse(substitute(x))),", ",(deparse(substitute(y))))
  x<-drop(as.matrix(x))
  y<-drop(as.matrix(y))
  
  if (!is.numeric(x)|!is.numeric(y)){
    stop("input data must be numeric")}
  
  if (any(c(x,y) < 0, na.rm = TRUE)){
    stop("input data must be non-negative")}
  
  if (any(is.na(c(x,y)))){
    x[is.na(x)]<-0
    y[is.na(y)]<-0
    warning("missing values replaced with zeroes")}
  
  alternative <- char.expand(alternative, c("two.sided", 
                                        "less", "greater","auto"))
  if (length(alternative) > 1L || is.na(alternative)){
    stop("alternative must be \"two.sided\", \"less\" or \"greater\"")}
  length_diff<-length(x)-length(y)
  if(length_diff>0){
    y<-c(y,rep(0,length_diff))
  }
  else if(length_diff<0){
    x<-c(x,rep(0,abs(length_diff)))
  }
  xy<-matrix(c(x,y),ncol=2)
  N<-apply(xy,2,sum)
  H<-(N*log(N, shanon.base)-apply(xy*log(xy,shanon.base), 2,sum,na.rm = TRUE))/N
  S<-(apply(xy*log(xy,shanon.base)**2, 2,sum,na.rm = TRUE) -
      ((apply(xy*log(xy,shanon.base), 2,sum,na.rm = TRUE)**2)/N))/(N**2)
  HutchesonTstat<- (diff(H[c(2,1)])-difference)/sqrt(sum(S))
  df<-(sum(S)**2)/(sum(S**2/N))
  estimate_dif<-diff(H[c(2,1)])
  if (alternative == "auto") {
    alternative <-if(estimate_dif<0){"less"}else{"greater"}}
  
  if (alternative == "less") {
    pval <- pt(HutchesonTstat, df)
  }
  else if (alternative == "greater") {
    pval <- pt(HutchesonTstat, df, lower.tail = FALSE)
  }
  else {
    pval <- 2 * pt(-abs(HutchesonTstat), df)
  }
  names(HutchesonTstat) <- "Hutcheson T"
  names(df) <- "df"
  names(H) <- c("Diversity of x","Diversity of y")
  mu<-difference
  names(mu)<-"difference in H'"
  rval <- list(statistic = HutchesonTstat, parameter = df, p.value = pval,
               estimate = H, null.value = mu,
               method="Hutcheson T-test for two communities",
              alternative = alternative,data.name=dname)
    class(rval) <- "htest"
    return(rval)
  }