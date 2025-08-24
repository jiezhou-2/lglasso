## code to prepare `DATASET` dataset goes here

usethis::use_data(sample_data, overwrite = TRUE,compress = "xz")

control=list(tau0=c(0.8,1))

usethis::use_data(control)
