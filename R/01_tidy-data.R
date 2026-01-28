# 01_tidy-data.R ####
# Author: Simon Dedman
# Date: 2026-01-28
# Purpose: Explore and clean shark data for analysis
#
# Description:
#   - Load raw shark data
#   - Generate histograms for numeric columns to identify outliers
#   - Generate bar plots for character columns to reveal patterns
#   - Filter data based on exploratory findings (gear type, country)
#
# Inputs:
#   - data/sharkdata.rda
#
# Outputs:
#   - outputs/01_explore-clean/histograms/*.png
#   - outputs/01_explore-clean/barplots/*.png


# Functions ####
## plot_numeric_histograms ####
#' Plot histograms for all numeric columns
#'
#' @param df Data frame to plot
#' @param output_dir Output directory path (will be created if needed)
#'
plot_numeric_histograms <- function(df, output_dir) {
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  df_name <- deparse(substitute(df))

  numeric_cols <- df |>
    dplyr::select(tidyselect::where(is.numeric)) |>
    names()

  purrr::walk(numeric_cols, \(col_name) {
    p <- ggplot2::ggplot(df, ggplot2::aes(x = .data[[col_name]])) +
      ggplot2::geom_histogram(bins = 30) +
      ggplot2::labs(
        title = paste("Histogram of", col_name),
        x = col_name,
        y = "Frequency"
      ) +
      ggplot2::theme_minimal()

    ggplot2::ggsave(
      filename = file.path(output_dir, paste0(df_name, "_", col_name, "_histogram.png")),
      plot = p,
      width = 6,
      height = 4
    )
  })
}

## plot_character_barplots ####
#' Plot bar plots for all character columns
#'
#' @param df Data frame to plot
#' @param output_dir Output directory path (will be created if needed)
#'
plot_character_barplots <- function(df, output_dir) {
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  df_name <- deparse(substitute(df))

  char_cols <- df |>
    dplyr::select(tidyselect::where(is.character)) |>
    names()

  purrr::walk(char_cols, \(col_name) {
    p <- ggplot2::ggplot(
      df,
      ggplot2::aes(x = forcats::fct_infreq(.data[[col_name]]))
    ) +
      ggplot2::geom_bar() +
      ggplot2::labs(
        title = paste("Bar plot of", col_name),
        x = col_name,
        y = "Count"
      ) +
      ggplot2::theme_minimal() +
      ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))

    ggplot2::ggsave(
      filename = file.path(output_dir, paste0(df_name, "_", col_name, "_barplot.png")),
      plot = p,
      width = 8,
      height = 5
    )
  })
}

## plot_spatial_map ####
#' Plot spatial distribution map
#'
#' @param df Data frame with lon/lat columns
#' @param output_path Full output file path
#' @param lon_col Name of longitude column (default "longitude")
#' @param lat_col Name of latitude column (default "latitude")
#' @param title Plot title
#'
plot_spatial_map <- function(df,
                             output_path,
                             lon_col = "longitude",
                             lat_col = "latitude",
                             title = "Spatial Distribution") {
  dir.create(dirname(output_path), recursive = TRUE, showWarnings = FALSE)

  world <- rnaturalearth::ne_countries(scale = "medium", returnclass = "sf")

  lon_vals <- df[[lon_col]]
  lat_vals <- df[[lat_col]]

  p <- ggplot2::ggplot() +
    ggplot2::geom_sf(data = world, fill = "antiquewhite") +
    ggplot2::geom_point(
      data = df,
      ggplot2::aes(x = .data[[lon_col]], y = .data[[lat_col]]),
      color = "blue",
      alpha = 0.5,
      size = 1
    ) +
    ggplot2::coord_sf(
      xlim = c(
        min(lon_vals) - 0.01 * abs(min(lon_vals)),
        max(lon_vals) + 0.01 * abs(max(lon_vals))
      ),
      ylim = c(
        min(lat_vals) - 0.01 * abs(min(lat_vals)),
        max(lat_vals) + 0.01 * abs(max(lat_vals))
      ),
      expand = FALSE
    ) +
    ggplot2::labs(
      title = title,
      x = "Longitude",
      y = "Latitude"
    ) +
    ggplot2::theme_minimal()

  ggplot2::ggsave(
    filename = output_path,
    plot = p,
    width = 8,
    height = 5
  )
}

usethis::use_package("rnaturalearth")


# Shark data (sharkdata) ####
## Load data ####
load(here::here("data", "sharkdata.rda"))


## Explore numeric columns ####
plot_numeric_histograms(
  sharkdata,
  here::here("outputs", "01_explore-clean", "histograms")
)

# Based on histograms, remove outliers from specific columns: check:
# lat ~27
# lon ~-80: both gear==Polyball, USA
# hook_buoy_number >60: gear==Gillnet


## Explore character columns ####
plot_character_barplots(
  sharkdata,
  here::here("outputs", "01_explore-clean", "barplots")
)

# Based on bar plots:
# country==USA is minority, drop
# gear != Drumline is minority, drop
# bait: could subset to bonito only?


## Filter data ####
sharkdata <- sharkdata |> # 1461 to 1364
  dplyr::filter(
    gear %in% c("Drumline-top", "Drumline-bottom"),
    country != "USA"
  )


## Spatial distribution ####
plot_spatial_map(
  sharkdata,
  here::here("outputs", "01_explore-clean", "sharkdata_spatial_distribution_cleaned.png"),
  title = "Shark Data Spatial Distribution After Cleaning"
)




# Static Explanatory Data (expvarstat) ####
## Load data ####
load(here::here("data", "expvarstat.rda"))

## Explore numeric columns ####
plot_numeric_histograms(
  expvarstat,
  here::here("outputs", "01_explore-clean", "histograms")
)

## Explore character columns ####
plot_character_barplots(
  expvarstat,
  here::here("outputs", "01_explore-clean", "barplots")
)

## Spatial distribution ####
plot_spatial_map(
  expvarstat,
  here::here("outputs", "01_explore-clean", "expvarstat_spatial_distribution.png"),
  title = "Static Explanatory Data Spatial Distribution"
)


# Dynamic Explanatory Data (expvardy) ####
## Load data ####
load(here::here("data", "expvardy.rda"))

## Explore numeric columns ####
plot_numeric_histograms(
  expvardy,
  here::here("outputs", "01_explore-clean", "histograms")
)

## Explore character columns ####
plot_character_barplots(
  expvardy,
  here::here("outputs", "01_explore-clean", "barplots")
)

## Spatial distribution ####
plot_spatial_map(
  expvardy,
  here::here("outputs", "01_explore-clean", "expvardy_spatial_distribution.png"),
  title = "Dynamic Explanatory Data Spatial Distribution"
)




# Join tables ####
## summarise to remove duplicates ####
expvarstatunique <- expvarstat |> # remove dupes, nrow 31473 before; 30055 after; 1418 removed. 159 vars before, 140 after, 19 removed
  dplyr::group_by(latitude, longitude) |>
  dplyr::summarise(
    dplyr::across(tidyselect::where(is.numeric), mean, na.rm = TRUE),
    dplyr::across(tidyselect::where(~ is.character(.x) | lubridate::is.POSIXt(.x)), dplyr::first)
  ) |> # library(lubridate)
  dplyr::mutate(
    dplyr::across(tidyselect::where(is.numeric), ~ ifelse(is.nan(.x), NA, .x)), # convert NaN to NA. POSIX needs lubridate
    dplyr::across(tidyselect::where(~ is.character(.x)), ~ ifelse(is.nan(.x), NA, .x))
  ) |> # https://community.rstudio.com/t/why-does-tidyrs-fill-work-with-nas-but-not-nans/25506/5
  dplyr::ungroup()

expvardyunique <- expvardy |> # remove dupes, nrow 31473 before; 30055 after; 1418 removed. 159 vars before, 140 after, 19 removed
  dplyr::group_by(latitude, longitude, event_dt, event_ts) |>
  dplyr::summarise(
    dplyr::across(tidyselect::where(is.numeric), mean, na.rm = TRUE),
    dplyr::across(tidyselect::where(~ is.character(.x) | lubridate::is.POSIXt(.x)), dplyr::first)
  ) |> # library(lubridate)
  dplyr::mutate(
    dplyr::across(tidyselect::where(is.numeric), ~ ifelse(is.nan(.x), NA, .x)), # convert NaN to NA. POSIX needs lubridate
    dplyr::across(tidyselect::where(~ is.character(.x)), ~ ifelse(is.nan(.x), NA, .x))
  ) |> # https://community.rstudio.com/t/why-does-tidyrs-fill-work-with-nas-but-not-nans/25506/5
  dplyr::ungroup()

sharkdatajoined <- sharkdata |> # 1364*26 -> 1364*33
  dplyr::left_join(
    expvarstatunique,
    by = c("latitude", "longitude")
  ) |>
  dplyr::left_join(
    expvardyunique,
    by = c("event_dt", "event_ts", "latitude", "longitude")
  )


## Explore numeric columns ####
plot_numeric_histograms(
  sharkdatajoined,
  here::here("outputs", "01_explore-clean", "histograms")
)

## Explore character columns ####
plot_character_barplots(
  sharkdatajoined,
  here::here("outputs", "01_explore-clean", "barplots")
)

## Spatial distribution ####
plot_spatial_map(
  sharkdatajoined,
  here::here("outputs", "01_explore-clean", "expvardy_spatial_distribution.png"),
  title = "Dynamic Explanatory Data Spatial Distribution"
)

# Save cleaned and joined data ####
usethis::use_data(sharkdatajoined, overwrite = TRUE)
