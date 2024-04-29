install.packages(c("DBI", "RMySQL"))
install.packages("sqldf")
install.packages("httr")
library(jsonlite)
library(DBI)
library(RMySQL)

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
      if (param$param_type == "varsome_liftover") {
        # overlift request, the response will be a chromosome variant id in hg38
        chr_list_liftover <- list()
        url <- param$url
        for (i in param$body) {
          url = paste(url, param$body[[i]])
          print(paste('url：', url) )
          response <- GET(param$url, content_type("application/json"), add_headers(Authorization=headers))
          status_code <- status_code(response)
          if(status_code == '200'){
            filt_data <- content(response, as = "text")
            chr_list_liftover[[length(chr_list_liftover) + 1]] <- filt_data
          }
        }
        return(chr_list_liftover)
      }else{
        response <- GET(param$url, content_type("application/json"), add_headers(Authorization=headers)) # nolint
        status_code <- status_code(response)
        print(response)
        if (status_code == 200) {
          filt_data <- content(response, as = "text")
          if (param$param_type == "mutation") {
            print('filt_data')
            print(filt_data)
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
            print(data_list)
            prediction_list <- list()
            obj <- list()
            for (i in seq_along(data_list)) {
              obj["transcript_stable"] <- data_list[[i]]$transcript_stable
              obj["prediction"] <- data_list[[i]]$prediction
              prediction_list <- c(prediction_list, list(obj))
            }
            # print(prediction_list)
            search_result <- list(
              param_type = param$param_type,
              prediction_list
            )
            return(search_result)
          }else if (param$param_type == "SIFT") {
            # SIFT
            # 用VarSome 调用
            library(jsonlite)
            library(xml2)
            stop_for_status(response)
          }else if (param$param_type == "oncokb") {
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
    }else{ # POST request
      # now only varsome use post request
      if (length(param[["token"]]) != 0) {
        headers <- c(
          Authorization=param$token,
          'Content-Type'="application/json"
        )
      }else {
        headers <- ""
      }
      response <- POST(param$url, body = param$body, add_headers(headers),encode = "json") 
      status_code <- status_code(response)
      if (status_code == 200) {
        # response_content <- content(response, "json")
        # response_data <- fromJSON(response)
        # print(paste('http post request success:', response))
        # cat(response)
        response_data <- as.character(response)
        folder_oath <- "data/"
        file_name <- c(paste(param$param_type, ".json"))
        file_path <- file.path(folder_oath, file_name)
        writeLines(response_data, file_path)
        # print(paste('response:',response_content))
        search_result <- list(
          param_type = param$param_type,
          write_status = 'succuss'
        )
        return(search_result)
        # print(paste('http post request success:', content))
      }else{
        print(response)
        print(paste("HTTP post request error:", status_code))
      }
    }
  }
}
save.image()
