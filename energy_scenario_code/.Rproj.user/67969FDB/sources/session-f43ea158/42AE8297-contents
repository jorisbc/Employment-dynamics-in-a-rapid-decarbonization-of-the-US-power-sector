# color stuff
coolwarm_hcl <- colorspace::diverging_hcl(11,
                                          h = c(250, 10), c = 100, l = c(37, 88), power = c(0.7, 1.7))
brbg_hcl <- colorspace::diverging_hcl(11,
                                      h = c(180, 50), c = 80, l = c(20, 95), power = c(0.7, 1.3))
viridis_hcl <- colorspace::sequential_hcl(11,
                                          h = c(300, 75), c = c(35, 95), l = c(15, 90), power = c(0.8, 1.2))
plasma_hcl <- colorspace::sequential_hcl(11,
                                         h = c(-100, 100), c = c(60, 100), l = c(15, 95), power = c(2, 0.9))

pal <- function(col, border = "transparent") {
  n <- length(col)
  plot(0, 0, type="n", xlim = c(0, 1), ylim = c(0, 1),
       axes = FALSE, xlab = "", ylab = "")
  rect(0:(n-1)/n, 0, 1:n/n, 1, col = col, border = border)
}
par(mar = rep(0, 4), mfrow = c(2, 2))
pal(coolwarm_hcl)
pal(viridis_hcl)
pal(brbg_hcl)
pal(plasma_hcl)

my_colors_l <- list(coolwarm_hcl[c(5,1,3,4,2)], # A
                    coolwarm_hcl[c(11,7,9,10,8)], # B
                    viridis_hcl[c(1,3,2)], # C
                    viridis_hcl[c(9,4)], # D
                    c(coolwarm_hcl[c(6)],'black'), # E
                    brbg_hcl[c(4,1,3,5)], # F
                    brbg_hcl[c(11,8,10)], # G
                    plasma_hcl[c(11,7,9,10,8)]) # H
my_colors <- unlist(my_colors_l)


colfull <- data.frame(
tech_full = c("Batteries", "BECCS", "Biomass", "Coal CCS", "Coal", "CSP", "PV dist.", 
              "Gas CC CCS", "Gas CC", "Gas CT", "Geothermal", "Hydro", "Nuclear", "O-G-S", 
              "PHS", "PV utility", "Wind off", "Wind on" ),
tech_colors = c(viridis_hcl[5],brbg_hcl[9],brbg_hcl[11],'darkgrey','black',plasma_hcl[11],plasma_hcl[9],
                plasma_hcl[4],plasma_hcl[5],plasma_hcl[6],coolwarm_hcl[8],'blue','magenta',coolwarm_hcl[6],
                coolwarm_hcl[1],viridis_hcl[11],brbg_hcl[1],brbg_hcl[3])  
)
# tech_short = c("Batteries", "Biomass", "Biomass", "Coal", "Coal", "CSP", "PV", 
#                "Gas", "Gas", "Gas", "Geothermal", "Hydro", "Nuclear", "Gas", 
#                "Hydro", "PV", "Wind", "Wind" )