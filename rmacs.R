titertable <- read.titerTable("2010_2020_raw_titer.csv")

map <- acmap(table = titertable)
map <- optimizeMap(
  map = map,
  number_of_dimensions = 2,
  number_of_optimizations = 100,
  minimum_column_basis = "none"
)

view(map)
selectedOptimization(map)
plot(map)


antigentic_distance<-mapDistances(map)