exponental_tilting_grad <- function( par, ffq, prec=0 )     {
  Xbeta<-    exp(ffq%*%par) 
  t(ffq) %*% exp(ffq%*%par) 
}
