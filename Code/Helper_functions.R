# P_INFL_OO6                    ################################################
# Script: Helper Functions      ################################################
# Author: Trang Tu              ################################################

if (!require("pacman", quietly = TRUE)) {
  install.packages("pacman")
}

library(pacman)
pacman::p_load(
  here,
  rio,
  dplyr, 
  ggplot2,
  tidyr,
  patchwork
)

################################################################################
# Function to clean character missingness --------------------------------------
clean_missing <- function(df) {
  char_cols <- sapply(df, is.character)
  
  df[char_cols] <- lapply(df[char_cols], function(x) {
    x <- trimws(x) # Remove leading and/or trailing whitespace from character strings
    x[x %in% c("", ".", "NA")] <- NA
    x
  })
  df
}


################################################################################
# Function to summarize ICA analysis results -----------------------------------

# Helper function for x-axis limit:
make_x_scale <- function(x, ref = NULL, min_width = 0.3) {
  
  vals <- c(x, ref)
  
  rng <- range(vals, na.rm = TRUE)
  
  xmin <- rng[1]
  xmax <- rng[2]
  
  current_width <- xmax - xmin
  
  if (current_width < min_width) {
    
    center <- mean(rng)
    
    xmin <- center - min_width / 2
    xmax <- center + min_width / 2
  }
  
  scale_x_continuous(
    limits = c(xmin, xmax),
    expand = expansion(mult = c(0.02, 0.05))
  )
}

make_panel_limits <- function(df, ref = NULL, min_width = 0.2) {
  
  out <- df %>%
    group_by(State) %>%
    summarise(
      xmin = min(Value, na.rm = TRUE),
      xmax = max(Value, na.rm = TRUE),
      .groups = "drop"
    )
  
  if (!is.null(ref)) {
    out <- out %>%
      left_join(ref, by = "State") %>%
      mutate(
        xmin = pmin(xmin, Truth),
        xmax = pmax(xmax, Truth)
      )
  }
  
  out %>%
    rowwise() %>%
    mutate(
      width = xmax - xmin,
      center = (xmin + xmax) / 2,
      xmin = ifelse(
        width < min_width,
        center - min_width / 2,
        xmin
      ),
      xmax = ifelse(
        width < min_width,
        center + min_width / 2,
        xmax
      )
    ) %>%
    ungroup()
}

# Main function
summarise_ICA <- function(
    dat,
    metric = "R2_H",
    bins = 50,
    alpha_levels = c(0.05, 0.10, 0.20),
    output_dir = here("Plot"),
    save_plot = TRUE,
    dpi = 300,
    width = 10,
    height = 7.3,
    print_summary = TRUE,
    print_plot = TRUE
) {
  
  # Summary statistics -------------------
  
  x <- dat[[metric]]
  
  # density-based mode
  dens <- density(x, na.rm = TRUE)
  mode_est <- dens$x[which.max(dens$y)]
  
  summary_tbl <- data.frame(
    Statistic = c(
      "Mean",
      "Median",
      "Mode",
      "Min",
      "Max"
    ),
    Value = c(
      round(mean(x, na.rm = TRUE), 3),
      round(median(x, na.rm = TRUE),3),
      round(mode_est, 3),
      round( min(x, na.rm = TRUE), 3),
      round(max(x, na.rm = TRUE), 3)
    )
  )
  
  # Symmetric density intervals
  for(a in alpha_levels){
    
    q <- quantile(
      x,
      probs = c(a/2, 1-a/2),
      na.rm = TRUE
    )
    
    summary_tbl <- rbind(
      summary_tbl,
      data.frame(
        Statistic =
          sprintf("SDI_%d%%", round((1-a)*100)),
        Value =
          sprintf("[%.3f, %.3f]", q[1], q[2])
      )
    )
  }
  
  # Histogram -------------------
  
  p_hist <- ggplot(dat,
                   aes(x = .data[[metric]])) +
    geom_histogram(
      bins = bins,
      fill = "grey80",
      color = "black"
    ) +
    theme_bw() +
    labs(
      x = if (metric == "R2_H") expression(R[H]^2) else metric,
      y = "Frequency"
    ) + 
    theme(
      plot.title = element_blank()) + 
    geom_vline(xintercept = mean(x, na.rm = TRUE), 
               color = "#bf3729", linewidth = 1.5) +
    geom_vline(
      xintercept = quantile(
        x,
        probs = c(0.1/2, 1-0.1/2),
        na.rm = TRUE
      )[1], 
      color = "black",
      linetype = "dashed",
      linewidth = 1
    ) +
    geom_vline(
      xintercept = quantile(
        x,
        probs = c(0.1/2, 1-0.1/2),
        na.rm = TRUE
      )[2], 
      color = "black",
      linetype = "dashed",
      linewidth = 1
    ) +
    scale_x_continuous(breaks = seq(0, 1, by =0.1))
  
  
  # Helper for annotation positions -------------------
  
  make_label_pos <- function(df) {
    df %>%
      group_by(State) %>%
      summarise(
        xmax =  max(Value), #quantile(Value, 0.95, na.rm = TRUE),
        ymax = max(density(Value)$y),
        .groups = "drop"
      )
  }
  
  # Mean distributions -------------------
  
  mean_dat <- dat %>%
    transmute(
      mean_Y_S0,
      mean_Y_S1
    ) %>%
    pivot_longer(
      everything(),
      names_to = "Variable",
      values_to = "Value"
    ) %>%
    mutate(
      State = ifelse(grepl("S0", Variable), "S0", "S1")
    )
  
  mean_ref <- data.frame(
    State = c("S0", "S1"),
    Truth = c(
      unique(dat$mean.S0)[1],
      unique(dat$mean.S1)[1]
    )
  )
  
  mean_est <- mean_dat %>%
    group_by(State) %>%
    summarise(
      Mean = mean(Value, na.rm = TRUE),
      .groups = "drop"
    )
  
  mean_labels <- mean_est %>%
    left_join(mean_ref, by = "State") %>%
    left_join(make_label_pos(mean_dat), by = "State") %>%
    mutate(
      label = paste0(
        "Est mean: ", round(Mean, 2),
        "\nEmp mean: ", round(Truth, 2)
      )
    )
  mean_limits <- make_panel_limits(
    mean_dat,
    mean_ref,
    min_width = 0.2
  )
  
  
  p_mean <- ggplot(mean_dat, aes(Value)) +
    geom_density(
      fill = "grey80",
      color = "black"
    ) +
    geom_vline(
      data = mean_est,
      aes(xintercept = Mean),
      color = "#2f357c",
      linewidth = 1
    ) +
    geom_vline(
      data = mean_ref,
      aes(xintercept = Truth),
      color = "#bf3729",
      linetype = "dashed",
      linewidth = 1
    ) +
    geom_label(
      data = mean_labels,
      aes(label = label),
      x = Inf,
      y = Inf,
      inherit.aes = FALSE,
      hjust = 1.05,
      vjust = 1.1,
      size = 3.5,
      alpha = 0.6
    ) +
    geom_blank(
      data = mean_limits,
      aes(x = xmin, y = 0)
    ) +
    geom_blank(
      data = mean_limits,
      aes(x = xmax, y = 0)
    ) +
    facet_wrap(~State, scales = "free") +
    labs(
      x = "Values of estimated sample mean",
      y = "Density"
    ) +
    coord_cartesian(clip = "off") +
    theme_bw()
  
  # Variance distributions -------------------
  
  var_dat <- dat %>%
    transmute(
      var_Y_S0,
      var_Y_S1
    ) %>%
    pivot_longer(
      everything(),
      names_to = "Variable",
      values_to = "Value"
    ) %>%
    mutate(
      State = ifelse(grepl("S0", Variable), "S0", "S1")
    )
  
  var_ref <- data.frame(
    State = c("S0", "S1"),
    Truth = c(
      unique(dat$var.S0)[1],
      unique(dat$var.S1)[1]
    )
  )
  
  var_est <- var_dat %>%
    group_by(State) %>%
    summarise(
      Mean = mean(Value, na.rm = TRUE),
      .groups = "drop"
    )
  
  var_labels <- var_est %>%
    left_join(var_ref, by = "State") %>%
    left_join(make_label_pos(var_dat), by = "State") %>%
    mutate(
      label = paste0(
        "Est var: ", round(Mean, 2),
        "\nEmp var: ", round(Truth, 2)
      )
    )
  var_limits <- make_panel_limits(
    var_dat,
    var_ref,
    min_width = 0.2
  )
  
  p_var <- ggplot(var_dat, aes(Value)) +
    geom_density(
      fill = "grey80",
      color = "black"
    ) +
    geom_vline(
      data = var_est,
      aes(xintercept = Mean),
      color = "#2f357c",
      linewidth = 1
    ) +
    geom_vline(
      data = var_ref,
      aes(xintercept = Truth),
      color = "#bf3729",
      linetype = "dashed",
      linewidth = 1
    ) +
    geom_label(
      data = var_labels,
      aes(label = label),
      x = Inf,
      y = Inf,
      inherit.aes = FALSE,
      hjust = 1.05,
      vjust = 1.1,
      size = 3.5,
      alpha = 0.6
    ) +
    geom_blank(
      data = var_limits,
      aes(x = xmin, y = 0)
    ) +
    geom_blank(
      data = var_limits,
      aes(x = xmax, y = 0)
    ) +
    facet_wrap(~State, scales = "free") +
    labs(
      x = "Values of estimated sample variance",
      y = "Density"
    ) +
    coord_cartesian(clip = "off") +
    theme_bw()
  
  # Deviance distributions -------------------
  
  dev_dat <- dat %>%
    select(dev_S0, dev_S1) %>%
    pivot_longer(
      everything(),
      names_to = "Variable",
      values_to = "Value"
    ) %>%
    mutate(
      State = ifelse(grepl("S0", Variable), "S0", "S1")
    )
  
  dev_est <- dev_dat %>%
    group_by(State) %>%
    summarise(
      Mean = mean(Value, na.rm = TRUE),
      .groups = "drop"
    )
  
  dev_labels <- dev_est %>%
    left_join(make_label_pos(dev_dat), by = "State") %>%
    mutate(
      label = paste0(
        "Mean deviance: ",
        round(Mean, 3)
      )
    )
  dev_limits <- make_panel_limits(
    dev_dat,
    min_width = 0.2
  )

  
  p_dev <- ggplot(dev_dat, aes(Value)) +
    geom_density(
      fill = "grey80",
      color = "black"
    ) +
    geom_vline(
      data = dev_est,
      aes(xintercept = Mean),
      color = "#2f357c",
      linewidth = 1
    ) +
    geom_label(
      data = dev_labels,
      aes(label = label),
      x = Inf,
      y = Inf,
      inherit.aes = FALSE,
      hjust = 1.05,
      vjust = 1.1,
      size = 3.5
    ) +
    geom_blank(
      data = dev_limits,
      aes(x = xmin, y = 0)
    ) +
    geom_blank(
      data = dev_limits,
      aes(x = xmax, y = 0)
    ) +
    facet_wrap(
      ~State,
      scales = "free"
    ) +
    labs(
      x = "Deviance value",
      y = "Density"
    ) +
    coord_cartesian(clip = "off") +
    theme_bw()
  
  
  # AIC distributions -------------------
  
  AIC_dat <- dat %>%
    select(AIC_S0, AIC_S1) %>%
    pivot_longer(
      everything(),
      names_to = "Variable",
      values_to = "Value"
    ) %>%
    mutate(
      State = ifelse(grepl("S0", Variable), "S0", "S1")
    )
  
  AIC_est <- AIC_dat %>%
    group_by(State) %>%
    summarise(
      Mean = mean(Value, na.rm = TRUE),
      .groups = "drop"
    )
  
  AIC_labels <- AIC_est %>%
    left_join(make_label_pos(AIC_dat), by = "State") %>%
    mutate(
      label = paste0(
        "Mean AIC: ",
        round(Mean, 0)
      )
    )
  AIC_limits <- make_panel_limits(
    AIC_dat,
    min_width = 0.2
  )
  
  
  p_AIC <- ggplot(AIC_dat, aes(Value)) +
    geom_density(
      fill = "grey80",
      color = "black"
    ) +
    geom_vline(
      data = AIC_est,
      aes(xintercept = Mean),
      color = "#2f357c",
      linewidth = 1
    ) +
    geom_label(
      data = AIC_labels,
      aes(label = label),
      x = Inf,
      y = Inf,
      inherit.aes = FALSE,
      hjust = 1.05,
      vjust = 1.1,
      size = 3.5,
      alpha = 0.6
    ) +
    geom_blank(
      data = AIC_limits,
      aes(x = xmin, y = 0)
    ) +
    geom_blank(
      data = AIC_limits,
      aes(x = xmax, y = 0)
    ) +
    facet_wrap(
      ~State,
      scales = "free"
    ) +
    labs(
      x = "AIC value",
      y = "Density"
    ) +
    coord_cartesian(clip = "off") +
    theme_bw()
  
  
  # Export the plots --------------
  if (isTRUE(save_plot)){
    data_name <- deparse(substitute(dat))
    
    if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
    
    save_plot_f <- function(plot, suffix, path = output_dir) {
      if (is.null(output_dir)) {
        return(invisible(NULL))
      }
      file_name <- file.path(
        output_dir,
        paste0(data_name, "_", suffix, ".png")
      )
      ggsave(
        filename = file_name,
        plot = plot,
        width = width,
        height = height,
        dpi = dpi,
        units = "cm",
        bg = "white"
      )
    }
    
    save_plot_f(
      p_hist,
      paste0(metric, "_histogram")
    )
    
    save_plot_f(
      p_mean,
      "mean_distribution"
    )
    
    save_plot_f(
      p_var,
      "variance_distribution"
    )
    
    save_plot_f(
      p_dev,
      "deviance_distribution"
    )
    
    save_plot_f(
      p_AIC,
      "AIC_distribution"
    )
    
  }
  
  
  # Display -------------------
  if( print_summary){
    print(summary_tbl)
    
  }
  if (print_plot){
    print(p_hist)
    print(p_dev)
    print(p_AIC)
    print(p_mean)
    print(p_var)
  }
  
  invisible(
    list(
      summary = summary_tbl,
      histogram = p_hist,
      dev_plot = p_dev,
      mean_plot = p_mean,
      var_plot = p_var
    )
  )
}
