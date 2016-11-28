safe.exp <- function (x) 
{
  xx = sign(x) * pmin(abs(x), 500)
  return(exp(xx))
}
