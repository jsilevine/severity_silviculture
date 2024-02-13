library(sf)

perims <- st_read("data/perimeters/all_perimiters.shp")
perims <- perims[,c("FIRE_NA", "YEAR_")]
names(perims)[1] <- "Fire_Name"
names(perims)[2] <- "Fire_Year"
perims$Start_Day <- 152
perims$End_Day <- 258
perims$Fire_ID <- c("DIXIE_2021", "NORTH_COMPLEX_2020", "SHEEP_2020", "WALKER_2019", "SUGAR_2021")

plot(perims[,1])

st_write(perims, "data/perimeters/gee_perimeters.shp", append = FALSE)
st_write(perims[perims$Fire_Name == "SUGAR",], "data/perimeters/gee_test.shp", append = FALSE)
