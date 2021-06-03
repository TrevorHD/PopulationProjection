##### Create transition matrix ----------------------------------------------------------------------------

# Function that creates a transition matrix given number and names of stages
A.create <- function(values, stage.names = NULL){
  
  # Let values be a vector of numerics representing the values of the transition matrix 
  # (from left to right and top to bottom)

  # Let stage.names be a vector of characters representing the names of stages or classes
  
  # Create matrix
  A <- matrix(data = values, nrow = sqrt(length(values)), byrow = TRUE)
  
  # Use stage/class names for rows and columns
  colnames(A) <- c(stage.names)
  rownames(A) <- c(stage.names)
  
  # Return the transition matrix
  return(A)}





##### Initial population size -----------------------------------------------------------------------------

# Function that creates a column vector of initial counts
A0.create <- function(values, stage.names = NULL){
  
  # Let values be a vector of numerics representing the values of the transition matrix 
  # (from left to right and top to bottom)
  
  # Let stage.names be a vector of characters representing the names of stages or classes
  
  # Create matrix
  A0 <- matrix(data = values, ncol = 1)
  
  # Use stage/class names for rows and columns
  rownames(A0) <- rownames(A)
  colnames(A0) <- "Abundance"
  
  # Return the transition matrix
  return(A0)}





##### Project population ----------------------------------------------------------------------------------

# Function to project population
pop.project <- function(steps, A, A0, output.raw = FALSE, output.integer = FALSE,
                        graph = "none"){
  
  # steps: number of time steps to project over
  # A: transition matrix with dimensions n x n
  # A0: initial abundance matrix with dimensions n x 1
  # output.raw: should raw outputs instead of a summary be given for each time step?
  # output.integer: should abundances be rounded to the nearest whole number?
  # graph: "abundance" for population counts, "proportions" for proportional abundance
  
  # Create empty matrix of counts
  counts <- matrix(0, nrow = 4, ncol = steps + 1)
  
  # Rename rows and columns
  rownames(counts) <- rownames(A0)
  colnames(counts) <- seq(0, steps)
  
  # Set first column as initial population
  counts[, 1] <- A0
  
  # Populate matrix of counts by calculating projecting population for each time step
  for(i in 2:(steps + 1)){
    counts[, i] <-  A %*% counts[, i - 1]}
  
  # Same as above, but calculates proportional abundance instead of absolute abundance
  counts.prop <- counts
  for(i in 0:(steps + 1)){
    for(j in 1:nrow(counts)){
      counts.prop[j, i] <- counts[j, i]/sum(counts[, i])}}
  
  # Calculate stable age distribution (when steps -> Inf)
  stable.prop <- Re(eigen(A)$vectors[, 1])/sum(Re(eigen(A)$vectors[, 1]))
  
  # Calculate discrete growth rate lambda for all time steps
  counts.lambda <- c()
  for(i in 0:steps){
    counts.lambda.val <- sum(counts[, i + 1])/sum(counts[, i])
    counts.lambda <- append(counts.lambda, counts.lambda.val)}
  
  # Growth rate at stable age distribution
  stable.lambda <- Re(eigen(A)$values[1])
  
  # Plot absolute abundances
  if(graph == "abundance"){
    par(mar = c(5, 5, 4, 2))
    plot(0, 0, pch = "", ylim = c(0, 1.1*max(counts)), xlim = c(0, steps + 1), 
         xlab = "Time",  ylab = "Abundance", xaxt = "n", cex.lab = 1.8, cex.axis = 1.3)
    for(i in 1:nrow(counts)){
      points(counts[i, ], col = rainbow(nrow(counts))[i], type = "l", lwd = 2)}
    axis(1, at = seq(1, steps + 1), labels = seq(0, steps), cex.axis = 1.3)
    legend("topleft", col = rainbow(nrow(counts)), lwd = rep(2, nrow(counts)), 
           legend = if(is.null(rownames(counts))){1:nrow(counts)} else {rownames(counts)}, 
           bty = "n")}
  
  # Plot relative abundances
  if(graph == "proportions"){
    par(mar = c(5, 5, 4, 2))
    plot(0, 0, pch = "", ylim = c(0, 1), xlim = c(0, steps + 1), 
         xlab = "Time",  ylab = "Proportional Abundance", xaxt = "n", cex.lab = 1.8, cex.axis = 1.3)
    for(i in 1:nrow(counts)){
      points(counts.prop[i, ], col = rainbow(nrow(counts.prop))[i], type = "l", lwd = 2)}
    axis(1, at = seq(1, steps + 1), labels = seq(0, steps), cex.axis = 1.3)
    legend("topleft", col = rainbow(nrow(counts.prop)), lwd = rep(2, nrow(counts.prop)), 
           legend = if(is.null(rownames(counts.prop))){1:nrow(counts.prop)} else {rownames(counts.prop)}, 
           bty = "n")}
  
  # Compile final and stable-state outputs
  out <- list(FinalAbundance = if(output.integer == FALSE){
                                 counts[, ncol(counts)]} else {round(counts[, ncol(counts)], 0)},
              FinalPropAbundance = round(counts.prop[, ncol(counts.prop)], 4),
              StablePropAbundance = round(stable.prop, 4),
              FinalLambda = round(counts.lambda[length(counts.lambda)], 4),
              StableLambda = round(stable.lambda, 4))
  
  # Compile raw outputs
  out.raw <- list(Abundance = if(output.integer == FALSE){counts} else {round(counts, 0)},
                  PropAbundance = round(counts.prop, 4),
                  Lambda = round(counts.lambda, 4))
  
  # Output values to console
  if(output.raw == FALSE){
    return(out)} else {return(out.raw)}}





##### Example run -----------------------------------------------------------------------------------------

# Create transition matrix with stage names
A <- A.create(c(0,    0,    3,    2,
                0.67, 0.61, 0,    0,
                0,    0.29, 0.72, 0,
                0,    0,    0.24, 0.95),
              stage.names = c("Juvenile", "Sub-adult", "Adult (young)", "Adult (old)"))

# Create matrix of initial abundances
A0 <- A0.create(c(50, 50, 50, 50),
                stage.names = c("Juvenile", "Sub-adult", "Adult (young)", "Adult (old)"))

# Run population projection
pop.project(steps = 10, A = A, A0 = A0, graph = "abundance",
            output.raw = FALSE, output.integer = TRUE)

