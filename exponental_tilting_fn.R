exponential_tilting_fn <- function( par, ffq,  prec=0)     {
  Xbeta<-    exp(ffq%*%par) 
  sum( exp(ffq%*%par)   ) 
  
}
