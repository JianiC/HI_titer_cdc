source("setup.R")

set.seed(850909)
titertable <- read.titerTable("2010_2020_raw_titer_filter.csv")

map <- acmap(table = titertable)
map <- optimizeMap(
  map = map,
  number_of_dimensions = 2,
  number_of_optimizations = 1000,
  minimum_column_basis = "none",
  options = list(dim_annealing = TRUE)
)

view(map)
mapCohesion(map)

 

titertable2 <- read.titerTable("2010_2020_raw_titer_filter2.csv")
map2 <- acmap(table = titertable2)
map2 <- optimizeMap(
  map = map2,
  number_of_dimensions = 2,
  number_of_optimizations = 2000,
  minimum_column_basis = "none",
  #options = list(dim_annealing = TRUE)
)
selectedOptimization(map2)

view(map2)

## diagnostics

agCohesion(map)

srCohesion(map)


optimizeAgReactivity(
  map2,
  optimization_number = 1,
  reactivity_stress_weighting = 1,
  fixed_ag_reactivities = rep(NA, numAntigens(map2 )),
  start_pars = rep(0, numAntigens(map2)),
  reoptimize = FALSE,
  number_of_optimizations = 100,
  options = list()
)


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
  options = list(dim_annealing = TRUE)
)

moveTrappedPoints(
  map2,
  optimization_number = 1,
  grid_spacing = 0.25,
  max_iterations = 10,
  options = list()
)

antigentic_distance<-mapDistances(map2)

## add color 
as.data.frame(agCoords(map, optimization_number = 1))%>%
  rownames_to_column()%>%
  setNames(c("ag_epi_isolate_id","Dim1","Dim2"))%>%
  merge(uga_titer_strain_filteryear,by=c("ag_epi_isolate_id","ag_epi_isolate_id"))%>%
  ggplot(aes(x=Dim1,y=Dim2,group=ag_collection_year, color=as.factor(ag_collection_year)))+
  geom_point()+
  #scale_color_brewer(palette="Set3")+
  theme_bw()
  