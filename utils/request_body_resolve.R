# function for varsome variant post request body
exchange_fun <- function(param_array){
  variants_char <- character()
  for (i in 2:length(param_array)-1) {
    if( i != 1){
      character_ <- param_array[[i]]
      character_ <- gsub("\\["," ",character_)
      character_ <- gsub("\\]"," ",character_)
      character_ <- gsub(","," ",character_)
      variants_char <- append(variants_char, character_)
    }
  }
  variants_char <- paste(variants_char, collapse = ',')
  variants_char <- paste("{'variants':","[",variants_char,"]}")
  variants_char <- gsub(" ","",variants_char)
  variants_char <- gsub("'",'"',variants_char)
  
  return(variants_char)
}

request_body <- function (param_array) { # change list to json charactor for post request use
  # the maximum request number is 1000, if the number over 1000, we should split them to groups
  new_list <- split(param_array, ceiling(seq_along(param_array) / 500))
  if(length(new_list) > 1){
    params_list <- list()
    # exchange each group to request body format, and pust it to list as a return
    for (i in 1:length(new_list)) {
      # new_list[i] each group
      variants_char <- exchange_fun(new_list[[i]])
      params_list <- append(params_list, variants_char)
    }
    return (params_list)
  }else{
    variants_char <- list()
    variants_char <- exchange_fun(param_array)
    return(variants_char)
  }
}

# function for liftover get request body
liftover_body_fun <- function(array){
  liftover_char <- list()
  for (i in 2:47) {
    item_ <- array[[i]]
    item_ <- gsub("\\["," ",item_)
    item_ <- gsub("\\]"," ",item_)
    item_ <- gsub(","," ",item_)
    item_ <- gsub("\"","",item_)
    item_ <- gsub(" ","",item_)
    liftover_char <- append(liftover_char, item_)
  }
  return(liftover_char)
}