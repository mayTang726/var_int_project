# install.packages(c("DBI", "RMySQL"))
install.packages("sqldf")
install.packages("httr")
library(jsonlite)
# library(DBI)
# library(RMySQL)

connect_function <- function(param) {
  if (!param$API) {
    print("inter !API part") # nolint
    connect <- dbConnect(RMySQL::MySQL(), dbname = param$dbname , host = param$host, port = param$port, user = param$user, password = param$password) # nolint
    library(sqldf)
    filt_data <- sqldf("select * from cosmic_mutantcensus_v99_grch37 where COSMIC_SAMPLE_ID = COSS2947897") # nolint
    print(filt_data)
    search_result <- list(
      param_type = param$param_type,
      filt_data
    )
    return(search_result)
  }else {
    print("inter API part")
    # the parameter2:undefined.MutationTaster  should use GET request , if others should use other request , the juadge shoud be changed # nolint: line_length_linter.
    library(httr)
    # library(jsonlite) # nolint
    if (param$req_type == "GET") { # GET request
      if (length(param[["token"]]) != 0) {
        headers <- param$token
      }else {
        headers <- ""
      }
      # varsome get request just for liftover
      if (param$param_type == "varsome") {
        # overlift request, the response will be a chromosome variant id in hg38
        chr_list_liftover <- list()
        for (i in param$body) {
          url = paste(param$url, i)
          url <- gsub(" ","",url)
          response <- GET(url, content_type("application/json"), add_headers(Authorization=headers))
          status_code <- status_code(response)
          if(status_code == '200'){
            # print('request success!')
            filt_data <- content(response, as = "parsed")
            # chr_list_liftover[[length(chr_list_liftover) + 1]] <- filt_data
            chr_list_liftover <- append(chr_list_liftover, filt_data)
          }else{
            print('resuest fail!')
          }
        }
        print(chr_list_liftover)
        return(chr_list_liftover)
      }else{
        print('进来')
        if (param$param_type == "mutation") {
          response_result <- list()
          for (i in param$body) {
            url <- ""
            url <- paste(param$url, i)
            url <- gsub(" ", "", url)
            print(url)
            response <- GET(url, content_type("application/json")) # nolint
            status_code <- status_code(response)
            print(response)
            if (status_code == 200) {
              filt_data <- content(response, as = "text")
              lines <- readLines(textConnection(filt_data))
              keys <- strsplit(lines[1], "\t")[[1]]
              data_list <- list()
              for (line in lines[-1]) {
                values <- strsplit(line, "\t")[[1]]
                obj <- list()
                for (i in seq_along(keys)) {
                  obj[keys[i]] <- values[i] # nolint
                }
                data_list <- c(data_list, list(obj))
              }
              response_result <- append(response_result, data_list)
              print(response_result)
            }else {
              stop("HTTP get request error: ", status_code)
            }
          }
          search_result <- list(
            param_type = param$param_type,
            response_data = response_result
          )
          return(search_result)
        }else{
          if (status_code == 200) {
            filt_data <- content(response, as = "text")
            if (param$param_type == "oncokb") {
              library(jsonlite)
              filt_data <- fromJSON(filt_data)
              if (filt_data$oncogenic == "Oncogenic" | filt_data$oncogenic == "Likely Oncogenic") {
                # it means (likely) pathagenic , score +0.5
                search_result <- list(
                  param_type = param$param_type,
                  filt_data
                )
                return(search_result)
              }
            }else if (param$param_type == "civic") {
              # 需要使用post请求
              # print(12345)
              # print(response)
            }
            
          }else {
            stop("HTTP get request error: ", status_code)
          }
        }
      }
    }else{ # POST request
      # now only varsome use post request
      headers <- ""
      if (length(param[["token"]]) != 0) {
        headers <- c(
          Authorization=param$token,
          'Content-Type'="application/json"
        )
      }
      if (param$param_type == "SIFT") {
        # request limit numer : 300,  so we need to split the number of the request body
        split_list <- split(param$body, ceiling(seq_along(param$body) / 300))
        response_list <- list()
        for (i in 1:length(split_list)) {
          ##  change list to json format as post request body
          sift_request_body_json <- toJSON(list(hgvs_notations = split_list[[i]]), pretty = TRUE)
          response <- POST(param$url, body = sift_request_body_json, add_headers(headers),encode = "json")
          status_code <- status_code(response)
          if (status_code == 200) {
            response_content <- content(response, as = "text", encoding = "UTF-8")
            response_list[[i]] <- fromJSON(response_content)
          }else {
            # print(response)
            print(paste("HTTP post request error:", status_code))
          }
        }
        
        if(length(response_list) != 0){
          search_result <- list(
            param_type = param$param_type,
            write_status = 'succuss',
            response_data = response_list
          )
        }
      }else if (param$param_type == "varsome"){
        response <- POST(param$url, body = param$body, add_headers(headers),encode = "json") 
        status_code <- status_code(response)
        if (status_code == 200) {
          response_data <- ""
          response_data <- as.character(response)
          search_result <- list(
            param_type = param$param_type,
            write_status = 'succuss',
            response_data = response_data
          )
          return(search_result)
          # print(paste('http post request success:', content))
        }else{
          # print(response)
          print(paste("HTTP post request error:", status_code))
        }
      }
    }
  }
}
save.image()
