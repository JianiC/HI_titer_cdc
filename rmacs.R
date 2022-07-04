titertable <- read.titerTable("2010_2020_raw_titer_filter.csv")

map <- acmap(table = titertable)
map <- optimizeMap(
  map = map,
  number_of_dimensions = 2,
  number_of_optimizations = 500,
  minimum_column_basis = "none"
)
selectedOptimization(map)

view(map)



antigentic_distance<-mapDistances(map)

## add color 
agCoords(map, optimization_number = 1)%>%
  mutate()
  merge(uga_titer_strain_filteryear,by=c("ag_strain_name","ag_strain_name"))
  ggplot(aes(),color=ag_colelction_year)+
  geom_point()