source("setup.R")


titertable <- read.titerTable("2010_2020_raw_titer_filter.csv")
titertable2 <- read.titerTable("2010_2020_raw_titer_filter2.csv")
map2 <- acmap(table = titertable2)
map2 <- optimizeMap(
  map = map2,
  number_of_dimensions = 2,
  number_of_optimizations = 5000,
  minimum_column_basis = "none"
)
selectedOptimization(map2)

view(map2)

## diagnostics

agCohesion(map)

srCohesion(map)

mapCohesion(map)


checkHemisphering(
  map2,
  optimization_number = 1,
  grid_spacing = 0.25,
  stress_lim = 0.1,
  options = list()
)

## dimension test 
dimensionTestMap(
  map2,
  dimensions_to_test = 1:5,
  test_proportion = 0.1,
  minimum_column_basis = "none",
  fixed_column_bases = rep(NA, numSera(map2)),
  number_of_optimizations = 100,
  replicates_per_dimension = 100,
  options = list()
)

antigentic_distance<-mapDistances(map2)

## add color 
as.data.frame(agCoords(map2, optimization_number = 1))%>%
  rownames_to_column()%>%
  setNames(c("ag_strain_name","Dim1","Dim2"))%>%
  merge(uga_titer_strain_filteryear,by=c("ag_strain_name","ag_strain_name"))%>%
  ggplot(aes(x=Dim1,y=Dim2,group=ag_collection_year, color=as.factor(ag_collection_year)))+
  geom_point()
  