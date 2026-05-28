#Function for bootstrapping on the basis of an individual ID when there are multiple rows per individual
#Example: If a given ID is sampled twice then all rows for that person will be repeated twice in the bootstrap sample

bootstrap_cluster <- function(data, id_var) {
  
  #unique IDs in data
  ids <- unique(data[[id_var]])

  #sample with replacement, and also assign a new bootstrap ID number which is unique
  boot_samp <- data.frame(
    boot_id = sample(ids, length(ids), replace = TRUE),
    unique_boot_id=1:length(ids)
  )
  
  #create bootstrap sample, retaining all rows in data for each individual for each time they are in boot_samp
  boot_data <- left_join(
    boot_samp,
    data,
    by = setNames(id_var, "boot_id"),
    relationship = "many-to-many"
  )
  
  #make id_var the unique bootstrap id
  boot_data[[id_var]]<-boot_data$unique_boot_id
  boot_data$unique_boot_id<-NULL
  boot_data$boot_id<-NULL
  
  return(boot_data)
  
}

#This is a test of the function

# df <- data.frame(
#   id = c(1,1,2,2,3,3),
#   y = c(10,11,20,21,30,31)
# )
# 
# bootstrap_cluster(df,"id")

