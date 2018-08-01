
pr <- data_idn$pr$positive / data_idn$pr$examined

ex <- raster::extract(data_idn$pop_raster, data_idn$pr[, 3:4], cellnumbers = TRUE)

cells_with_mulitple_points <- unique(ex[, 'cells'][duplicated(ex[, 'cells'])])


cell_range <- rep(NA, length(cells_with_mulitple_points))


big_enough <- data_idn$pr$examined >= 500

for(i in seq_along(cell_range)){
  subset_bool <- ex[, 'cells'] == cells_with_mulitple_points[i]
  pr_sub <- pr[subset_bool & big_enough]
  cell_range[i] <- max(dist(pr_sub))
}

cell_range <- cell_range[is.finite(cell_range)]

max(cell_range)

