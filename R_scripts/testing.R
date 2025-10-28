library(rbenchmark)

benchmark(
  with_pivot = {
    cols <- colnames(gapfilling_inputs_chunked)
    cols <- cols[!cols %in% c("yday", "cell_no", "state")]
    gapfilling_inputs_chunked %>% 
      pivot_longer(
        cols = all_of(cols), 
        names_to = "variable", 
        values_to = "value", 
        names_transform = list(variable = as.factor)
      ) %>% 
      group_by(state, cell_no, variable) %>% 
      reframe(missing_vals = sum(is.na(value))/400) %>% 
      pivot_wider(names_from = variable, values_from = missing_vals)
    },
    no_pivot = {
      gapfilling_inputs_chunked %>% 
        group_by(state, cell_no) %>% 
        reframe(
          hz = sum(is.na(hz))/400,
          I_input = sum(is.na(I_input))/400,
          Kd_490 = sum(is.na(Kd_490))/400,
          S_input = sum(is.na(S_input))/400,
          T_input = sum(is.na(T_input))/400,
          UV_input = sum(is.na(UV_input))/400
        )
    }
)

