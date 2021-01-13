#----------------------------------------------------#
# Plot results from validate theory pipeline
#----------------------------------------------------#

# Load results.

# Please specify the name of the directory containing the saved data here.
dir_name <- "20210110121248"
save_directory <- paste0("data/",dir_name)
for(file_name in list.files(save_directory))
{
  load(paste0(save_directory,"/",file_name))
}

# Data for plotting
distributions <- configs$distributions
seed_location <- configs$seed_location
rhos <- sapply(distributions,FUN = function(distribution){
  attributes(distribution$partition(seed_location))[["rho"]]
  }
)
avg_empirical_conductance <- rowMeans(empirical_conductances)

# Plots
plot(x = rhos,
     y = avg_empirical_conductance,
     type = "l")
lines(x = rhos,
      y = population_conductances,
      type = "l")
