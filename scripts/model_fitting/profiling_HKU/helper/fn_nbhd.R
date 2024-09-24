# define the function to get the neighbourhood
fn_nbhd <- function(object, time, unit) {
  nbhd_list = list()
  if(time > 1){
    nbhd_list = c(nbhd_list, list(c(unit, time-1)))
  }
  if(unit > 1){
    nbhd_list = c(nbhd_list, list(c(unit-1, time)))
  }
  return(nbhd_list)
}