toDate <- function(times){
    if (is.numeric(times)) return(times)
    times * 365.2425 + as.Date("0000-01-01")
}
