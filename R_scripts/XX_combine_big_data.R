suppressMessages(suppressWarnings({
  library(stringr)
  library(magrittr)
  library(dplyr)
  library(arrow)
  library(oceanmap)
  library(ggplot2)
  library(conflicted)
  library(rasterVis)
  library(qs)
  library(qs2)
}))
conflicts_prefer(dplyr::filter(), dplyr::select(), terra::extract(), terra::adjacent())

fixnum <- function(n, digits = 4) {str_flatten(c(rep("0", digits-nchar(as.character(n))), as.character(n)))}
clean_memory <- function() {
  gc()
  terra::terraOptions(memfrac=0.5)
}
extract_adj <- function(r, cells) {
  adj_cells <- sapply(X = r, FUN = adjacent, cells = cells, directions = "queen")[,1] %>% unique()
  adj_vals <- matrix(nrow = length(adj_cells), ncol = 365, data = NA)
  for (ac in 1:length(adj_cells)) {
    adj_vals[ac,] <- extract(r, y = adj_cells[ac]) %>% unlist() %>% unname()
  }
  colMeans(adj_vals, na.rm = T)
}

# BARRA-C2 land fraction dist mask
load(file.path("D:", "FRDC-Seaweed-Raw-Data", "BARRA-C2", "BARRA_C2_cell_mask.Rdata")) # distance from coast cell_mask
cells <- ncell(cell_mask) # total number of cells
# mask_rast <- terra::rast(cell_mask)

# Pre-load which cells are valid based on the cell mask
cell_list <-  as.data.frame(cell_mask)
cell_list$ncell <- 1:cells
cell_list <- cell_list %>% 
  filter(!is.na(layer) & layer < 97.5) %>% 
  # This will become the metadata of how the adjacentcy searches went
  mutate(salt_adjs = 0, salt_success = NA,
         u_adjs = 0, u_success = NA,
         v_adjs = 0, v_success = NA)
rm(cell_mask, cells)
clean_memory()

cell_mat <- matrix(nrow = nrow(cell_list), ncol = 365, data = NA)

# Get BARRA data
BARRA_C2_rsdsdir_data <- terra::rast(file.path("data_raw", "BARRA-C2", "BARRA_C2_rsdsdir_data.tif"))
BARRA_C2_ts_data <- terra::rast(file.path("data_raw", "BARRA-C2", "BARRA_C2_ts_data.tif"))

# Get BRAN data - has already been projected (not masked)
BRAN_salt_data <- terra::rast(file.path("data_raw", "BRAN2023", "BRAN_salt_data.tif"))
BRAN_u_data <- terra::rast(file.path("data_raw", "BRAN2023", "BRAN_u_data.tif"))
BRAN_v_data <- terra::rast(file.path("data_raw", "BRAN2023", "BRAN_v_data.tif"))

overwrite <- F
cell_list$fname <- NA
for (i in 1:nrow(cell_list)) {
  cell_list$fname[i] <- file.path("data", "cell_inputs", paste0("cell_inputs_", fixnum(cell_list$ncell[i], 8), ".parquet"))
}

# First try -------------------------------------------------------------------------------------------------------
# On the first try, cell_inputs are only saved if there is data available for all variables (T, I, S, U and V)
progress_list <- cell_list %>% 
  filter(is.na(salt_success))
all_cells <- progress_list$ncell
prog_track <- as.integer(seq(0, length(all_cells), length.out = 101))

for (c in 1:length(all_cells)) {
  cell_inputs <- data.frame(
    rsdsdir = extract(BARRA_C2_rsdsdir_data, all_cells[c]) %>% unlist() %>% unname(),
    ts = extract(BARRA_C2_ts_data, all_cells[c]) %>% unlist() %>% unname()
  )
  
  # For each value, if no valid data is found the script searches up to three cells away (~16 km) in all directions before giving up 
  salt <- extract(BRAN_salt_data, all_cells[c]) %>% unlist() %>% unname()
  u <- extract(BRAN_u_data, all_cells[c]) %>% unlist() %>% unname()
  v <- extract(BRAN_v_data, all_cells[c]) %>% unlist() %>% unname()

  # Final success checks
  cell_list$salt_success[cell_list$ncell == all_cells[c]] <- ifelse(all(is.na(salt)), F, T)
  cell_list$u_success[cell_list$ncell == all_cells[c]] <- ifelse(all(is.na(u)), F, T)
  cell_list$v_success[cell_list$ncell == all_cells[c]] <- ifelse(all(is.na(v)), F, T)
  
  if(all(!is.na(salt)) & all(!is.na(u)) & all(!is.na(v))) {
    cell_inputs$salt <- salt
    cell_inputs$u <- u
    cell_inputs$v <- v
    cell_inputs$cell_ID <- all_cells[c]
    arrow::write_parquet(cell_inputs, progress_list$fname[progress_list$ncell == all_cells[c]])
  }
  if (c %in% prog_track){
    print(paste0(c, " of ", length(all_cells), " cells done, ", round(100*c/length(all_cells),1), "% finished at ", Sys.time()))
  }
}

# Second try ------------------------------------------------------------------------------------------------------
# On the second try, only cells without saved data are tried again. cell_inputs are saved regardless of what data is available
# progress_list <- cell_list %>% 
#   filter(!salt_success & )
all_cells <- progress_list$ncell
prog_track <- as.integer(seq(0, length(all_cells), length.out = 101))

for (c in 1:length(all_cells)) {
  fname <- file.path("data", "cell_inputs", paste0("cell_inputs_", fixnum(all_cells[c], 8), ".parquet"))
  
  if (!file.exists(fname)) {
    cell_inputs <- data.frame(
      rsdsdir = extract(BARRA_C2_rsdsdir_data, all_cells[c]) %>% unlist() %>% unname(),
      ts = extract(BARRA_C2_ts_data, all_cells[c]) %>% unlist() %>% unname()
    )
    
    # For each value, if no valid data is found the script searches up to three cells away (~16 km) in all directions before giving up 
    salt <- extract(BRAN_salt_data, all_cells[c]) %>% unlist() %>% unname()
    if (all(is.na(salt))) {
      # First try - get cells adjacent to original cell
      salt <- extract_adj(r = BRAN_salt_data, cells = all_cells[c])
      if (all(is.na(salt))) {
        # Second try - get cells adjacent to the first try cells
        adj_cells <- sapply(X = BRAN_salt_data, FUN = adjacent, cells = all_cells[c], directions = "queen")[,1] %>% unique()
        salt <- extract_adj(r = BRAN_salt_data, cells = adj_cells)
        if (all(is.na(salt))) {
          # Third try - get cells adjacent to the second try cells
          adj_cells <- sapply(X = BRAN_salt_data, FUN = adjacent, cells = adj_cells, directions = "queen")[,1] %>% unique()
          salt <- extract_adj(r = BRAN_salt_data, cells = adj_cells)
          cell_list$salt_adjs[cell_list$ncell == all_cells[c]] <- 3
        } else {cell_list$salt_adjs[cell_list$ncell == all_cells[c]] <- 2}
      } else {cell_list$salt_adjs[cell_list$ncell == all_cells[c]] <- 1}
    }
    
    u <- extract(BRAN_u_data, all_cells[c]) %>% unlist() %>% unname()
    if (all(is.na(u))) {
      # First try - get cells adjacent to original cell
      u <- extract_adj(r = BRAN_u_data, cells = all_cells[c])
      if (all(is.na(u))) {
        # Second try - get cells adjacent to the first try cells
        adj_cells <- sapply(X = BRAN_u_data, FUN = adjacent, cells = all_cells[c], directions = "queen")[,1] %>% unique()
        u <- extract_adj(r = BRAN_u_data, cells = adj_cells)
        if (all(is.na(u))) {
          # Third try - get cells adjacent to the second try cells
          adj_cells <- sapply(X = BRAN_u_data, FUN = adjacent, cells = adj_cells, directions = "queen")[,1] %>% unique()
          u <- extract_adj(r = BRAN_u_data, cells = adj_cells)
          cell_list$u_adjs[cell_list$ncell == all_cells[c]] <- 3
        } else {cell_list$u_adjs[cell_list$ncell == all_cells[c]] <- 2}
      } else {cell_list$u_adjs[cell_list$ncell == all_cells[c]] <- 1}
    }
    
    v <- extract(BRAN_v_data, all_cells[c]) %>% unlist() %>% unname()
    if (all(is.na(v))) {
      # First try - get cells adjacent to original cell
      v <- extract_adj(r = BRAN_v_data, cells = all_cells[c])
      if (all(is.na(v))) {
        # Second try - get cells adjacent to the first try cells
        adj_cells <- sapply(X = BRAN_v_data, FUN = adjacent, cells = all_cells[c], directions = "queen")[,1] %>% unique()
        v <- extract_adj(r = BRAN_v_data, cells = adj_cells)
        if (all(is.na(v))) {
          # Third try - get cells adjacent to the second try cells
          adj_cells <- sapply(X = BRAN_v_data, FUN = adjacent, cells = adj_cells, directions = "queen")[,1] %>% unique()
          v <- extract_adj(r = BRAN_v_data, cells = adj_cells)
          cell_list$v_adjs[cell_list$ncell == all_cells[c]] <- 3
        } else {cell_list$v_adjs[cell_list$ncell == all_cells[c]] <- 2}
      } else {cell_list$v_adjs[cell_list$ncell == all_cells[c]] <- 1}
    }
    
    # Final success checks
    cell_list$salt_success[cell_list$ncell == all_cells[c]] <- ifelse(all(is.na(salt)), F, T)
    cell_list$u_success[cell_list$ncell == all_cells[c]] <- ifelse(all(is.na(u)), F, T)
    cell_list$v_success[cell_list$ncell == all_cells[c]] <- ifelse(all(is.na(v)), F, T)
    
    cell_inputs$salt <- salt
    cell_inputs$u <- u
    cell_inputs$v <- v
    cell_inputs$cell_ID <- all_cells[c]
    
    arrow::write_parquet(cell_inputs, fname)
  }
  if (c %in% as.integer(seq(0, length(all_cells), length.out = 101))){
    print(paste0(c, " of ", length(all_cells), " cells done, ", round(100*c/length(all_cells),1), "% finished at ", Sys.time()))
  }
}

arrow::write_parquet(cell_list, "data/cell_list_metadata.parquet")

# This is how to judge whether a cell is empty based on its contents extracted from the rasterlayer (possibly helpful for the BRAN data)
# !is.na(terra::extract(mask_rast, 1000)$layer) & terra::extract(mask_rast, 1000)$layer < 97.5



