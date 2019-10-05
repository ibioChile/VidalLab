suppressMessages({
  library("sleuth")
})

setwd("/Users/pamelacamejo/Documents/IBIO/Rodrigo_Gutierrez/Projects/Arabidopsis_nitrate/Data")
sample_id <- dir(file.path("canales2014.kalisto"))
kal_dirs <- file.path(sample_id, 'abundance.h5')

