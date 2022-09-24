#' Simple implementation of particle swarm optimization
#' 
#' @seealso https://en.wikipedia.org/wiki/Particle_swarm_optimization
#' 
#' @examples
#' 
#' FUN <- function(x) mean(abs(x - 0.5))
#' 
#' pso(
#'   n_dim = 20,
#'   n_particles = 10,
#'   max_iter = 1000,
#'   inertia = 0.1,
#'   step_size_local = 0.1,
#'   step_size_global = 0.1,
#'   min = 0,
#'   max = 1,
#'   FUN = FUN
#' )
#' 
#' result <- pso(
#'   n_dim = 2,
#'   n_particles = 5,
#'   max_iter = 1000,
#'   inertia = 0.5,
#'   step_size_local = 0.1,
#'   step_size_global = 0.9,
#'   min = 0,
#'   max = 1,
#'   FUN = FUN
#' )
#' result$global_best_value
#' lapply(result$positions, function(x) any(x < 0))
#' result$positions[[999]]
#' result$global_best_position
#' 
#' ls_lines <- vector(mode = "list", length = 5)
#' for (i in 1:5) {
#'   ls_lines[[i]] <- matrix(
#'     unlist(lapply(result$positions[1:20], function(x) x[i,])),
#'     ncol = 2, byrow = TRUE
#'   )
#' }
#' 
#' plot(c(0,1), c(1, 0), col = "white")
#' abline(h = 0.5, lty = 3, col = "grey")
#' abline(v = 0.5, lty = 3, col = "grey")
#' 
#' for (i in 1:5) {
#'   lines(ls_lines[[i]], col = i)
#'   points(ls_lines[[i]], col = i, pch = 19)
#' }
pso <- function(n_dim,
                n_particles,
                max_iter,
                inertia = 0.1,
                step_size_local = 0.1,
                step_size_global = 0.1,
                min = 0,
                max = 1,
                FUN,
                ...) {
  
  starting_positions <- runif(n = n_particles * n_dim, min = min, max = max)
  particle_positions <- matrix(
    data = starting_positions, ncol = n_dim, nrow = n_particles
  )
  starting_positions <- particle_positions
  
  particle_positions_best <- particle_positions
  
  particle_position_value <- apply(X = particle_positions, 1, FUN = FUN)
  particle_position_value_best <- particle_position_value
  
  global_best_value <- min(particle_position_value)
  global_best_position_idx <- which.min(particle_position_value)
  global_best_position <- particle_positions[global_best_position_idx, ]
  global_best_position_matrix <- matrix(
    data = global_best_position, byrow = TRUE, ncol = n_dim, nrow = n_particles
  )
  
  particle_velocity <- runif(
    n = n_particles,
    min = -abs(max - min),
    max = abs(max - min)
  )
  particle_dimension_velocity <- matrix(
    data = particle_velocity, ncol = n_dim, nrow = n_particles, byrow = FALSE
  )
  
  iter <- 1
  positions_list <- vector(mode = "list", length = max_iter)
  values_list <- vector(mode = "list", length = max_iter)
  velocity_list <- vector(mode = "list", length = max_iter)
  
  while (iter <= max_iter && global_best_value > 0.01) {
    if (iter %% 10 == 0) {
      message(sprintf("Iteration %s of %s", iter, max_iter))
    }
    
    local_shift <- matrix(
      runif(n = n_particles * n_dim, min = 0, max = 1), 
      ncol = n_dim, nrow = n_particles
    )
    global_shift <- matrix(
      runif(n = n_particles * n_dim, min = 0, max = 1), 
      ncol = n_dim, nrow = n_particles
    )
    
    particle_dimension_velocity <- inertia * particle_dimension_velocity +
      step_size_local * local_shift * (starting_positions - particle_positions) +
      step_size_global * global_shift * (global_best_position_matrix - particle_positions)
    
    # limit the velocity to not move outside borders of search space
    particle_dimension_velocity <- pmin(particle_dimension_velocity, max - particle_positions)
    particle_dimension_velocity <- pmax(particle_dimension_velocity, min - particle_positions)
    
    particle_positions <- particle_positions + particle_dimension_velocity
    particle_position_value <- apply(X = particle_positions, 1, FUN = FUN)
    
    improvement_idx <- which(particle_position_value < particle_position_value_best)
    
    particle_positions_best[improvement_idx, ] <- particle_positions[improvement_idx, ]
    
    if (min(particle_position_value) < global_best_value) {
      global_best_value <- min(particle_position_value)
      global_best_position_idx <- which.min(particle_position_value)
      global_best_position <- particle_positions[global_best_position_idx, ]
      global_best_position_matrix <- matrix(
        data = global_best_position, byrow = TRUE, ncol = n_dim, nrow = n_particles
      )
    }
    
    positions_list[[iter]] <- particle_positions
    values_list[[iter]] <- particle_position_value
    velocity_list[[iter]] <- particle_dimension_velocity
    iter <- iter + 1
  }
  
  return(
    list(
      global_best_value = global_best_value,
      global_best_position = global_best_position,
      global_best_position_matrix = global_best_position_matrix,
      starting_positions = starting_positions,
      positions = positions_list,
      values = values_list,
      velocities = velocity_list
    )
  )
}
