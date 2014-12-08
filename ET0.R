# Berechnung der Grasreferenzverdunstung

net_radiation <- function (date, latitude, Sunshine, Temp, U, albedo) {
  library("lubridate")
  phi <- latitude
  alpha <- albedo
  
  xi <-  0.0172 * yday(date) - 1.39        
  
  R_0 <- 2425 + 1735 * sin(xi) + 44 * (phi - 51.0) * (sin(xi)-1)  
  
  S_0 <- 12.3 + sin(xi) * (4.3 + (phi - 51.0)/6)
  
  S_r <- Sunshine/S_0 
  
  R_G <- R_0 * (0.19 + 0.55 * S_r)
  
  R_L <- 10.8 + 0.205 * Temp 
  
  RnL <- R_L * (1.64*R_G/R_0-0.22) * (0.34 - 0.0044 * sqrt( U * es(Temp) ))
  
  L <- 249.8 - 0.242 * Temp
  
  RnK <- (1-alpha)  * R_G/L
  
  Rn <- RnK - RnL
  return(Rn)
}

gammastern <- function(v_2, gamma) gamma *(1+0.34*v_2)

es <- function(Temp) 6.11*exp(17.62 * Temp / (243.12 + Temp))

s <- function(Temp) es(Temp)* 4284/(243.12 + Temp)^2

ET0 <- function(Temp, v_2, Datum, Breitengrad, Sonnenscheindauer, U, albedo, gamma){
	s(Temp)/(s(Temp) + gammastern(v_2, gamma)) * 
		net_radiation(Datum, Breitengrad, Sonnenscheindauer, Temp, U, albedo) +
		90*gamma / (s(Temp) + gammastern(v_2, gamma)) * 
		es(Temp) / (Temp + 273) *
		(1-U/100) * v_2
}



# Beispiel (wird nicht ausgeführt, wenn die Datei gesourced wird)
if(FALSE){
  Temp <- 30 # °C
  U <- 60 # %
  v_2 <- 3 # m/sec
  Sonnenscheindauer <- 10 # h
  Breitengrad <- 51 # °
  Datum <- as.Date("2014-06-01")
  
  gamma <- 0.65 # hPa/K
  albedo <- 0.23
  
  ET0(Temp, v_2, Datum, Breitengrad, Sonnenscheindauer, U, albedo, gamma)
}
