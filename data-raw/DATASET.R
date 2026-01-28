## code to prepare `DATASET` dataset goes here
sharkdata <- readr::read_csv("data-raw/mydata.csv")
expvarstat <- readr::read_csv("data-raw/expvar-static.csv")
expvardy <- readr::read_csv("data-raw/expvar-dynamic.csv")

usethis::use_data(sharkdata, overwrite = TRUE)
usethis::use_data(expvarstat, overwrite = TRUE)
usethis::use_data(expvardy, overwrite = TRUE)
