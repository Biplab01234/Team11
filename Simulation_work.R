
# Install devtools 
install.packages("devtools")

# Install spatstat.core from GitHub
devtools::install_github("spatstat/spatstat.core")

install.packages("FNN")

# Install the INLA package from the INLA repository

options(timeout = 3000)  
install.packages("INLA", repos = "https://inla.r-inla-download.org/R/stable", dep = TRUE)

install.packages("fmesher")
install.packages("INLA", repos=c(getOption("repos"), INLA="https://inla.r-inla-download.org/R/stable"), dep=TRUE)
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install(c("graph", "Rgraphviz"), dep=TRUE)

install.packages("rgeos")



# Loading Necessary Libraries
library(sp)
library(raster)
library(spdep)       # For spatial dependencies
library(MASS)        # For statistical methods
library(INLA)        # For integrated nested Laplace approximations
library(fields)      # For tools for spatial data
library(sf)
library(spatstat)
library(dplyr)
library(osmdata)
library(ggplot2)
library(spatstat.core)
library(FNN)


# Set the seed for reproducibility
set.seed(123)

# Function to generate random points within a polygon boundary
generate_points <- function(n, poly){
  points <- spsample(poly, n, "random")
  return(coordinates(points))
}

# Function to create a SpatialPolygons object from an extent
create_polygon <- function(ext) {
  coords <- matrix(c(ext@xmin, ext@ymin,
                     ext@xmin, ext@ymax,
                     ext@xmax, ext@ymax,
                     ext@xmax, ext@ymin,
                     ext@xmin, ext@ymin), ncol = 2, byrow = TRUE)
  p <- Polygon(coords)
  ps <- Polygons(list(p), ID = "1")
  sp_poly <- SpatialPolygons(list(ps))
  return(sp_poly)
}

# Define polygons
north_poly <- create_polygon(raster::extent(494, 496, 2150, 2156))
south_poly <- create_polygon(raster::extent(496, 500, 2145, 2150))

# Generate theft locations
theft_locations_north <- generate_points(689, north_poly)
theft_locations_south <- generate_points(3327, south_poly)

# Simulate covariates for the 90 areal units
covariates <- data.frame(
  Pop15 = rnorm(90, mean=1000, sd=200),
  Apart = rnorm(90, mean=500, sd=100),
  Eco = rnorm(90, mean=700, sd=150),
  Employ = rnorm(90, mean=600, sd=120),
  inBorn = rnorm(90, mean=800, sd=160),
  Health = rnorm(90, mean=900, sd=180),
  Scholar = round(rnorm(90, mean=8, sd=1)),
  Extor = rpois(90, lambda=2),
  Murder = rpois(90, lambda=1),
  Burg = rpois(90, lambda=3),
  Shop = rpois(90, lambda=4),
  Public = rpois(90, lambda=5),
  Street = rpois(90, lambda=6),
  Kidnap = rpois(90, lambda=0.5)
)

# Simulate spatial coordinates for the 90 blocks (this is arbitrary)
set.seed(123)  # Ensure reproducibility
covariates$x <- runif(90, min = 494, max = 500)
covariates$y <- runif(90, min = 2145, max = 2156)

# Convert theft locations to data frames
north_theft_locations_df <- data.frame(x = theft_locations_north[, 1], y = theft_locations_north[, 2])
south_theft_locations_df <- data.frame(x = theft_locations_south[, 1], y = theft_locations_south[, 2])

# Find the nearest block for each theft location in the North region
nearest_indices_north <- get.knnx(cbind(covariates$x, covariates$y), cbind(north_theft_locations_df$x, north_theft_locations_df$y), k = 1)$nn.index
north_theft_locations_df <- cbind(north_theft_locations_df, covariates[nearest_indices_north, ])

# Repeat for the South region
nearest_indices_south <- get.knnx(cbind(covariates$x, covariates$y), cbind(south_theft_locations_df$x, south_theft_locations_df$y), k = 1)$nn.index
south_theft_locations_df <- cbind(south_theft_locations_df, covariates[nearest_indices_south, ])

# Convert to ppp objects with marks (covariates)
north_theft_ppp <- spatstat.geom::ppp(x = north_theft_locations_df$x, 
                                      y = north_theft_locations_df$y, 
                                      window = spatstat.geom::owin(xrange = range(north_theft_locations_df$x), 
                                                                   yrange = range(north_theft_locations_df$y)),
                                      marks = north_theft_locations_df[, -c(1, 2)])

south_theft_ppp <- spatstat.geom::ppp(x = south_theft_locations_df$x, 
                                      y = south_theft_locations_df$y, 
                                      window = spatstat.geom::owin(xrange = range(south_theft_locations_df$x), 
                                                                   yrange = range(south_theft_locations_df$y)),
                                      marks = south_theft_locations_df[, -c(1, 2)])





# Perform linear modeling for each covariate separately for the North region
model_north_pop15 <- lm(y ~ Pop15, data = north_theft_locations_df)
model_north_apart <- lm(y ~ Apart, data = north_theft_locations_df)
model_north_eco <- lm(y ~ Eco, data = north_theft_locations_df)
model_north_employ <- lm(y ~ Employ, data = north_theft_locations_df)
model_north_inBorn <- lm(y ~ inBorn, data = north_theft_locations_df)
model_north_health <- lm(y ~ Health, data = north_theft_locations_df)
model_north_scholar <- lm(y ~ Scholar, data = north_theft_locations_df)
model_north_extor <- lm(y ~ Extor, data = north_theft_locations_df)
model_north_murder <- lm(y ~ Murder, data = north_theft_locations_df)
model_north_burg <- lm(y ~ Burg, data = north_theft_locations_df)
model_north_shop <- lm(y ~ Shop, data = north_theft_locations_df)
model_north_public <- lm(y ~ Public, data = north_theft_locations_df)
model_north_street <- lm(y ~ Street, data = north_theft_locations_df)
model_north_kidnap <- lm(y ~ Kidnap, data = north_theft_locations_df)

# Summarize the model results for the North region
summary(model_north_pop15)
summary(model_north_apart)
summary(model_north_eco)
summary(model_north_employ)
summary(model_north_inBorn)
summary(model_north_health)
summary(model_north_scholar)
summary(model_north_extor)
summary(model_north_murder)
summary(model_north_burg)
summary(model_north_shop)
summary(model_north_public)
summary(model_north_street)
summary(model_north_kidnap)



# Perform linear modeling for each covariate separately for the South region
model_south_pop15 <- lm(y ~ Pop15, data = south_theft_locations_df)
model_south_apart <- lm(y ~ Apart, data = south_theft_locations_df)
model_south_eco <- lm(y ~ Eco, data = south_theft_locations_df)
model_south_employ <- lm(y ~ Employ, data = south_theft_locations_df)
model_south_inBorn <- lm(y ~ inBorn, data = south_theft_locations_df)
model_south_health <- lm(y ~ Health, data = south_theft_locations_df)
model_south_scholar <- lm(y ~ Scholar, data = south_theft_locations_df)
model_south_extor <- lm(y ~ Extor, data = south_theft_locations_df)
model_south_murder <- lm(y ~ Murder, data = south_theft_locations_df)
model_south_burg <- lm(y ~ Burg, data = south_theft_locations_df)
model_south_shop <- lm(y ~ Shop, data = south_theft_locations_df)
model_south_public <- lm(y ~ Public, data = south_theft_locations_df)
model_south_street <- lm(y ~ Street, data = south_theft_locations_df)
model_south_kidnap <- lm(y ~ Kidnap, data = south_theft_locations_df)

# Summarize the model results for the South region
summary(model_south_pop15)
summary(model_south_apart)
summary(model_south_eco)
summary(model_south_employ)
summary(model_south_inBorn)
summary(model_south_health)
summary(model_south_scholar)
summary(model_south_extor)
summary(model_south_murder)
summary(model_south_burg)
summary(model_south_shop)
summary(model_south_public)
summary(model_south_street)
summary(model_south_kidnap)



# Function to identify significant covariates from a list of linear models
identify_significant_covariates <- function(model_list, significance_level = 0.05) {
  significant_covariates <- list()
  
  for (model_name in names(model_list)) {
    model <- model_list[[model_name]]
    p_values <- summary(model)$coefficients[, "Pr(>|t|)"]
    
    # Exclude the intercept
    p_values <- p_values[names(p_values) != "(Intercept)"]
    
    # Find significant covariates
    significant_covariates[[model_name]] <- names(p_values[p_values < significance_level])
  }
  
  return(significant_covariates)
}

# Create a list of models for the North region
models_north <- list(
  model_north_pop15 = model_north_pop15,
  model_north_apart = model_north_apart,
  model_north_eco = model_north_eco,
  model_north_employ = model_north_employ,
  model_north_inBorn = model_north_inBorn,
  model_north_health = model_north_health,
  model_north_scholar = model_north_scholar,
  model_north_extor = model_north_extor,
  model_north_murder = model_north_murder,
  model_north_burg = model_north_burg,
  model_north_shop = model_north_shop,
  model_north_public = model_north_public,
  model_north_street = model_north_street,
  model_north_kidnap = model_north_kidnap
)

# Create a list of models for the South region
models_south <- list(
  model_south_pop15 = model_south_pop15,
  model_south_apart = model_south_apart,
  model_south_eco = model_south_eco,
  model_south_employ = model_south_employ,
  model_south_inBorn = model_south_inBorn,
  model_south_health = model_south_health,
  model_south_scholar = model_south_scholar,
  model_south_extor = model_south_extor,
  model_south_murder = model_south_murder,
  model_south_burg = model_south_burg,
  model_south_shop = model_south_shop,
  model_south_public = model_south_public,
  model_south_street = model_south_street,
  model_south_kidnap = model_south_kidnap
)


# Identify significant covariates for the North region
significant_north <- identify_significant_covariates(models_north)

# Identify significant covariates for the South region
significant_south <- identify_significant_covariates(models_south)

# Print significant covariates
print("Significant Covariates - North:")
print(significant_north)

print("Significant Covariates - South:")
print(significant_south)





# Function to perform stepwise selection based on BIC
perform_stepwise_selection <- function(initial_model) {
  step(initial_model, direction = "both", k = log(nrow(initial_model$model)))
}

# Applying stepwise selection for the North region
north_models_stepwise <- lapply(models_north, perform_stepwise_selection)

# Applying stepwise selection for the South region
south_models_stepwise <- lapply(models_south, perform_stepwise_selection)

# Function to extract model formula
extract_model_formula <- function(model_list) {
  sapply(model_list, function(model) as.character(formula(model)))
}

# Extracting final model formulas after stepwise selection
final_formulas_north <- extract_model_formula(north_models_stepwise)
final_formulas_south <- extract_model_formula(south_models_stepwise)

# Print final model formulas
print("Final Model Formulas - North:")
print(final_formulas_north)

print("Final Model Formulas - South:")
print(final_formulas_south)



#LGCP

# Assuming 'significant_north' contains the significant covariates for the North region

# Creating the formula for LGCP model for the North region
covariates_north_lgcp <- unlist(significant_north)
formula_lgcp_north <- as.formula(paste("y ~", paste(covariates_north_lgcp, collapse = " + ")))

# Fit LGCP model for the North region
model_lgcp_north <- inla(formula_lgcp_north, family = "lgcp", data = north_theft_locations_df)

# Repeat the process for the South region
covariates_south_lgcp <- unlist(significant_south)
formula_lgcp_south <- as.formula(paste("y ~", paste(covariates_south_lgcp, collapse = " + ")))

# Fit LGCP model for the South region
model_lgcp_south <- inla(formula_lgcp_south, family = "lgcp", data = south_theft_locations_df)




# Convert your data to a format suitable for INLA
# Assuming north_theft_locations_df and south_theft_locations_df have coordinates and covariates
north_points <- ppp(north_theft_locations_df$x, north_theft_locations_df$y, window = owin(c(min(north_theft_locations_df$x), max(north_theft_locations_df$x)), c(min(north_theft_locations_df$y), max(north_theft_locations_df$y))))
south_points <- ppp(south_theft_locations_df$x, south_theft_locations_df$y, window = owin(c(min(south_theft_locations_df$x), max(south_theft_locations_df$x)), c(min(south_theft_locations_df$y), max(south_theft_locations_df$y))))

# Specify the LGCP model formula
# Note: 'f' represents a spatial random field and is a key part of specifying an LGCP
formula_lgcp_north <- y ~ f(spatial.field, model = inla.spde2.pcmatern) + covariate1 + covariate2 # Add your significant covariates

# Fit the LGCP model for the North region
model_lgcp_north <- inla(formula_lgcp_north, data = as.data.frame(north_points), family = "poisson", control.compute = list(config = TRUE))

# Repeat the process for the South region with its respective formula and data






#NHPP
# Creating the formula for NHPP model for the North region
covariates_north_nhpp <- unlist(significant_north)
formula_nhpp_north <- as.formula(paste("y ~", paste(covariates_north_nhpp, collapse = " + ")))

# Fit NHPP model for the North region
model_nhpp_north <- inla(formula_nhpp_north, family = "poisson", data = north_theft_locations_df)

# Repeat the process for the South region
covariates_south_nhpp <- unlist(significant_south)
formula_nhpp_south <- as.formula(paste("y ~", paste(covariates_south_nhpp, collapse = " + ")))

# Fit NHPP model for the South region
model_nhpp_south <- inla(formula_nhpp_south, family = "poisson", data = south_theft_locations_df)






# Bayesian Inference for LGCP and NHPP models

# Load the INLA library if not already loaded
library(INLA)

# Inference for LGCP model for the North region
summary_lgcp_north <- summary(model_lgcp_north)
print(summary_lgcp_north)

# Inference for LGCP model for the South region
summary_lgcp_south <- summary(model_lgcp_south)
print(summary_lgcp_south)

# Inference for NHPP model for the North region
summary_nhpp_north <- summary(model_nhpp_north)
print(summary_nhpp_north)

# Inference for NHPP model for the South region
summary_nhpp_south <- summary(model_nhpp_south)
print(summary_nhpp_south)






# Calculate Euclidean distances between theft and recovery points
calculate_distances <- function(theft_points, recovery_points){
  distances <- numeric(length = nrow(recovery_points))
  for (i in 1:nrow(recovery_points)) {
    distances[i] <- sqrt(sum((theft_points[i, ] - recovery_points[i, ])^2))
  }
  return(distances)
}

# Calculate the distances for the subset of theft locations
distances <- calculate_distances(all_theft_locations[sampled_indices, ], recovery_points)



# Plotting Figure 1: Theft locations in the North and South regions
plot_theft_locations <- function(north_poly, south_poly, theft_locations_north, theft_locations_south) {
  # North region plot
  plot(north_poly, col='grey', main='North')
  points(theft_locations_north, col='blue', pch=20)
  
  # South region plot
  plot(south_poly, col='grey', main='South')
  points(theft_locations_south, col='blue', pch=20)
}

plot_theft_locations(north_poly, south_poly, theft_locations_north, theft_locations_south)

# Plotting Figure 2: Recovery locations and Histogram of distances
plot_recovery_and_histogram <- function(recovery_points, distances) {
  # Recovery locations plot
  plot(north_poly, col='grey', main='Recovered') # assuming recovery points are in the north for simplicity
  points(recovery_points, col='red', pch=20)
  
  # Histogram of distances
  hist(distances, breaks=50, main='Histogram of Distance', xlab='Distance (km)', ylab='Frequency')
}

plot_recovery_and_histogram(recovery_points, distances)








# Function to generate random points within a polygon boundary for Belo Horizonte
generate_points_bh <- function(n, poly){
  points <- spsample(poly, n, "random")
  return(coordinates(points))
}

# Assuming Belo Horizonte as a rectangular area for simplicity
bh_poly <- as(raster::extent(XXX, YYY, ZZZ, WWW), 'SpatialPolygons') # Replace with actual extent

# Generate theft locations for recovered cars in Belo Horizonte
theft_locations_bh <- generate_points_bh(5250, bh_poly)

# Function to generate recovery points close to theft points in Belo Horizonte
generate_recovery_points_bh <- function(theft_points, n){
  recoveries <- matrix(nrow = n, ncol = 2)
  for (i in 1:n){
    # Distribution to simulate closeness
    if(runif(1) < 0.15){ # Approx. 15% within 200m
      distance <- runif(1, min = 0, max = 0.2) # Max 200m
    } else {
      distance <- runif(1, min = 0.2, max = 5) # Arbitrarily chosen upper limit
    }
    angle <- runif(1, 0, 2 * pi)
    recoveries[i, ] <- theft_points[i, ] + c(cos(angle) * distance, sin(angle) * distance)
  }
  return(recoveries)
}

# Generate recovery locations in Belo Horizonte
recovery_points_bh <- generate_recovery_points_bh(theft_locations_bh, 5250)





spatial.r
Displaying spatial.r.