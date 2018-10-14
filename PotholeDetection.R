library(jpeg)
library(data.table)

SegmentPothole <- function(img){
  dim_pothole <- dim(img)
  
  pothole_grey <- matrix(data=0, nrow = dim_pothole[1], 
                         ncol=dim_pothole[2])
  
  sapply(1:dim_pothole[1], function(i){
    sapply(1:dim_pothole[2], function(j){
      pothole_grey[i,j] <<- mean(img[i,j,])
    })
  })
  
  pothole_vec <- round(as.double(pothole_grey)*255)
  
  val_count <- as.data.frame(table(pothole_vec))
  setDT(val_count)
  val_count$pothole_vec <- as.numeric(as.character(val_count$pothole_vec))
  
  hist_maxima <- unlist(matrix(val_count[which.max(Freq)])[,1])
  
  val_count[, perp_dist := abs(hist_maxima[2] * pothole_vec - 
                                 hist_maxima[1] * Freq)/sqrt(sum(hist_maxima^2))]
  
  grey_thresh <- val_count[pothole_vec < hist_maxima[1]][which.max(perp_dist)]$pothole_vec/255
  
  pothole_segmented <- array(data=0, dim = c(dim_pothole[1], 
                                             dim_pothole[2], 3))
  
  pothole_segmented_mat <- matrix(data=0, nrow = dim_pothole[1], 
                                  ncol =dim_pothole[2])
  
  sapply(1:dim_pothole[1], function(i){
    sapply(1:dim_pothole[2], function(j){
      if(mean(img[i,j,]) <= grey_thresh){
        pothole_segmented[i,j,] <<- 1
        pothole_segmented_mat[i, j] <<- 1
      } 
    })
  })
  
  return(list(pothole_grey, pothole_segmented, pothole_segmented_mat))
  
}

ExtractPotholeRegions <- function(img){
  
  dim_pothole <- dim(img)
  
  segment_list <- SegmentPothole(img)
  
  pothole_grey <- segment_list[[1]]
  pothole_segmented <- segment_list[[2]]
  pothole_segmented_mat <- segment_list[[3]]
  
  row_seq <- seq(1, dim_pothole[1], 9)
  col_seq <- seq(1, dim_pothole[2], 9)
  
  pothole_mat <- matrix(data=0, nrow = dim_pothole[1], 
                        ncol=dim_pothole[2])
  
  pothole_sec <- array(data=0, dim = c(dim_pothole[1], 
                                       dim_pothole[2], 3))
  
  
  sapply(1:(length(row_seq)-1), function(i){
    sapply(1:(length(col_seq)-1), function(j){
      curr_section <- as.double(pothole_segmented_mat[row_seq[i]:min(row_seq[i+1]-1, dim_pothole[1]),
                                                      col_seq[j]:min(col_seq[j+1]-1, dim_pothole[2])])
      
      
      if(sum(curr_section)/length(curr_section) > 0.15){
        pothole_mat[row_seq[i]:min(row_seq[i+1]-1, dim_pothole[1]),
                    col_seq[j]:min(col_seq[j+1]-1, dim_pothole[2])] <<- 1
        pothole_sec[row_seq[i]:min(row_seq[i+1]-1, dim_pothole[1]),
                    col_seq[j]:min(col_seq[j+1]-1, dim_pothole[2]), ] <<- 1
      }
      
    })
  })
  
  return(list(pothole_grey, pothole_segmented, pothole_mat, pothole_sec))
  
}


GetPotholeFeatures <- function(img){
  
  dim_pothole <- dim(img)
  
  extract_list <- ExtractPotholeRegions(img)
  
  LOG_filter <- matrix(data= c(0,0,3,2,2,2,3,0,0,
                               0,2,3,5,5,5,3,2,0,
                               3,3,5,3,0,3,5,3,3,
                               2,5,3,-12,-23,-12,3,5,2,
                               2,5,0,-23,-40,-23,0,5,2,
                               2,5,3,-12,-23,-12,3,5,2,
                               3,3,5,3,0,3,5,3,3,
                               0,2,3,5,5,5,3,2,0,
                               0,0,3,2,2,2,3,0,0), nrow=9, ncol=9)
  
  pothole_LOG <- matrix(data=0, nrow = dim_pothole[1], 
                        ncol=dim_pothole[2])
  
  pothole_grey <- extract_list[[1]]
  pothole_mat <- extract_list[[3]]
  
  sapply(1:(dim_pothole[1]-8), function(i){
    sapply(1:(dim_pothole[2]-8), function(j){
      curr_LOG <- pothole_grey[i:(i+8), j:(j+8)] %*% LOG_filter
      
      pothole_LOG[i+4, j+4] <<- curr_LOG[5,5]
    })
  })
  
  pothole_count <- sum(pothole_mat)
  pothole_grey_vals <- matrix(data=0, nrow = pothole_count, ncol=1)
  nopothole_grey_vals <- matrix(data=0, nrow = dim_pothole[1] * dim_pothole[2] -pothole_count, ncol=1)
  pothole_LOG_vals <- matrix(data=0, nrow = pothole_count, ncol=1)
  nopothole_LOG_vals <- matrix(data=0, nrow = dim_pothole[1] * dim_pothole[2] -pothole_count, ncol=1)
  
  curr_ph_index <- 1
  curr_noph_index <- 1
  
  sapply(1:dim_pothole[1], function(i){
    sapply(1:dim_pothole[2], function(j){
      if(pothole_mat[i,j] == 1){
        pothole_grey_vals[curr_ph_index, 1] <<- pothole_grey[i, j]
        pothole_LOG_vals[curr_ph_index, 1] <<- pothole_LOG[i, j]
        curr_ph_index <<- curr_ph_index+1
      }else{
        nopothole_grey_vals[curr_noph_index, 1] <<- pothole_grey[i, j]
        nopothole_LOG_vals[curr_noph_index, 1] <<- pothole_LOG[i, j]
        curr_noph_index <<- curr_noph_index+1
      }
      
    })
  })
  
  sd_grey <- sd(pothole_grey_vals[, 1]) - sd(nopothole_grey_vals[, 1])
  sd_LOG <- sd(pothole_LOG_vals[, 1]) - sd(nopothole_LOG_vals[, 1])
  
  return(list(extract_list[[2]], extract_list[[4]], sd_grey, sd_LOG))
  
}


GetPotHoleFreaturesAll <- function(input_path){
  
  file_names <- list.files(input_path, patter='*.jpg')
  
  results <- NULL
  
  lapply(file_names, function(file_name){
    
    print(file_name)
    
    curr_image <- readJPEG(paste0(input_path, file_name))
    curr_features <- GetPotholeFeatures(curr_image)
    
    curr_table <- data.table(file_name = file_name,
                             sd_grey = curr_features[[3]],
                             sg_LOG = curr_features[[4]],
                             pred = ifelse(curr_features[[3]] > 0 && 
                                             curr_features[[4]] > 0, 
                                           'pothole',
                                           'nopothole'))
    
    results <<- rbind.data.frame(results, curr_table)
  })
  
  return(results)
  
}

pothole_results <- GetPotHoleFreaturesAll('/Users/apple/Desktop/PotholeDetection/pothole_images/')

nopothole_results <- GetPotHoleFreaturesAll('/Users/apple/Desktop/PotholeDetection/nopothole_images/')


#################################
pothole <- readJPEG('/Users/apple/Desktop/PotholeDetection/pothole_images/pothole8.jpg')

final_list <- GetPotholeFeatures(pothole)

writeJPEG(final_list[[1]], '/Users/apple/Desktop/PotholeDetection/pothole_segmented.jpg')

writeJPEG(final_list[[2]], '/Users/apple/Desktop/PotholeDetection/pothole_sec.jpg')


nopothole <- readJPEG('/Users/apple/Desktop/PotholeDetection/nopothole_images/nopothole1.jpg')

final_list_no <- GetPotholeFeatures(nopothole)

writeJPEG(final_list_no[[1]], '/Users/apple/Desktop/PotholeDetection/nopothole_segmented.jpg')

writeJPEG(final_list_no[[2]], '/Users/apple/Desktop/PotholeDetection/nopothole_sec.jpg')
