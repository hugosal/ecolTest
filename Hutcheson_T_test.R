

#funcion y prueba de funcion t de hutchinson para diferencias
#en indice de Shanon entre dos comunidades

#hugo Salinas
#14/9/2020

Hutcheson.T.test<-function(x, y,
                            conf.level = 0.95, shanon.base = exp(1),
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
    warning("Missing values replaced with zeroes")}
  
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
  
  xy<-matrix(c(x,y),ncol=2) #juntar los dos vectores en una sola matriz
  #calcular indice shanon para las dos de una vez

  N<-apply(xy,2,sum)# numero de individuos por comunidad

  H<-(N*log(N, shanon.base)-apply(xy*log(xy,shanon.base), 2,sum,na.rm = TRUE))/N
  
  S<-(apply(xy*log(xy,shanon.base)**2, 2,sum,na.rm = TRUE) -
        
      ((apply(xy*log(xy,shanon.base), 2,sum,na.rm = TRUE)**2)/N))/(N**2)

  HutchesonTstat<- (diff(H[c(2,1)])-difference)/sqrt(sum(S))
  
  df<-(sum(S)**2)/(sum(S**2/N)) #en el texto creo que los redondea, aqui ya no lo hago

  
  #lo que sigue es nada mas para el print bonito de prueba de hip
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
  #names(estimate)<-"Difference in H'"
  names(H) <- c("Diversity of x","Diversity of y")
  mu<-difference
  names(mu)<-"difference in H'"
  rval <- list(statistic = HutchesonTstat, parameter = df, p.value = pval,
               estimate = H, null.value = mu,
               stderr = stderr, method="Hutcheson T-test for two communities",
              alternative = alternative,data.name=dname)
    class(rval) <- "htest"
    return(rval)
  }
# 
#library(vegan)
#data(BCI)
#comunidad<-BCI #datos de vegan


# #Ejemplo de Zar:
# tab<-matrix(c(47,35,7,5,3,2, 48,23,11,13,8,2),nrow = 2,byrow = T)
# comunidad<-data.frame(tab)
# colnames(comunidad)<-c("oak","corn","black","beech","cherry", "other")
# 
# 
# Hutcheson.T.test(x = comunidad[1,] ,y = comunidad[2,],shanon.base = 10,
#                   alternative = "two.sided") #mismo resultado de Zar
# 
# 
# #pruebas de otras direcciones, que creo que son consistentes
# Hutcheson.T.test(x = comunidad[1,] ,y = comunidad[2,],shanon.base = 10,
#                   alternative = "greater")
# 
# Hutcheson.T.test(x = comunidad[1,] ,y = comunidad[2,],shanon.base = 10,
#                   alternative = "less")
