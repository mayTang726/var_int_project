# function for varsome variant post request body
request_body <- function (array) {
  # change list to json charactor for post request use
  variants_char <- character()
  for (i in 2:length(array)-1) {
    if( i != 1){
      character_ <- array[[i]]
      character_ <- gsub("\\["," ",character_)
      character_ <- gsub("\\]"," ",character_)
      character_ <- gsub(","," ",character_)
      variants_char <- append(variants_char, character_)
    }
  }
  print(variants_char)
  variants_char <- paste(variants_char, collapse = ',')
  variants_char <- paste("{'variants':","[",variants_char,"]}")
  variants_char <- gsub(" ","",variants_char)
  variants_char <- gsub("'",'"',variants_char)
  return (variants_char)
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