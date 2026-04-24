# Analysis of sports records

Supplementary material, data and code for the paper "Analysis of sport records: A predictive approach".

# Data

We collected the annual records (since 2021) and world records in some atheltics and aquatics disciplines by manually scraping the following two website:
* https://worldathletics.org/
* https://www.worldaquatics.com/

The data can be loaded with:
```R
data <- read.csv("records.csv")
```

# Data analysis

```R
require(tidyverse)
source("libARGumbel.R")
# for reproducibility in calibration methods
set.seed(1234)

df <- data |>
  filter(grepl("^\\d+ metres$", type)) |>
  as_tibble() |>
  group_by(gender, type, wr) |>
  reframe(
    # the column "est" will contain all the estimated models
    est = list(argumbel$new(y))
  ) |>
  rowwise() |>
  mutate(
    # the column "calib" will contain all the calibration objects
    calib = list(calibration$new(est, nboot = 300)),
    # the predicted probability of breaking the world record follows:
    # estimative method
    pbwr_est = 1-calib$pest(wr),
    # calibrated probability method
    pbwr_CP = 1-calib$pCP(wr),
    # calibrated quantiles method
    pbwr_CQ = 1-calib$pCQ(wr)
  ) |>
  ungroup() |>
  select(gender, type, wr, starts_with("pbwr_"))

df |>
  mutate(distance = as.numeric(stringr::str_extract(type, "\\d+")), .after=gender) |>
  select(-type) |>
  arrange(gender, distance) |>
  knitr::kable()
```

|gender | distance|        wr|  pbwr_est|   pbwr_CP|   pbwr_CQ|
|:------|--------:|---------:|---------:|---------:|---------:|
|men    |      100| 10.438413| 0.0267524| 0.0377474| 0.0422536|
|men    |      200| 10.422095| 0.0829305| 0.1072238| 0.1158028|
|men    |      400|  9.295840| 0.1438132| 0.1740427| 0.1801209|
|men    |      800|  7.927856| 0.1428723| 0.1723210| 0.1811742|
|men    |    10000|  6.365372| 0.0740983| 0.0978116| 0.1031630|
|women  |      100|  9.532889| 0.0226110| 0.0313772| 0.0362805|
|women  |      200|  9.372071| 0.0531532| 0.0700532| 0.0750856|
|women  |      400|  8.403361| 0.0288948| 0.0386277| 0.0463125|
|women  |      800|  7.062147| 0.0501746| 0.0661526| 0.0739420|
|women  |    10000|  5.766547| 0.0639978| 0.0769602| 0.0901594|
