## TF properties
load_tf_census = function(){
  message("Load tf census")
  tf_census = readRDS("inst/extdata/annotations/tf_annotation.rds") %>%
    distinct(tf) %>%
    pull(tf)

  return(tf_census)
}
