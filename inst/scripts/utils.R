## TF properties
load_tf_census = function(){
  message('Load TF census')
  tf_census = read_delim('inst/extdata/annotations/human_tf_census.txt',
                         delim = "\t", col_names = "tf") %>%
    distinct(tf) %>%
    pull(tf)
  
  return(tf_census)
}