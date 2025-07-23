# Preprocessing MovieLensToy data


## Libraries
library(magrittr)
library(dplyr)
library(reshape2)


## import data 
df <- read.table("dev/MovieLens100.data")
colnames(df) <- c("user_id", "item_id", "rating", "timestamp")


## filter to have at least 200 ratings per user and 300 ratings per item
nb_user_min <- 200
nb_item_min <- 300

users_keep <- as.vector(which(table(df$user_id) >= nb_user_min))
length(users_keep)
items_keep <- as.vector(which(table(df$item_id) >= nb_item_min))
length(items_keep)

df_sub <- df %>% 
  filter(user_id %in% users_keep, 
         item_id %in% items_keep)


# reshape into dataframe 
MovieLensToy <- as.matrix(dcast(df_sub, user_id ~ item_id, value.var = "rating")[,-1])
usethis::use_data(MovieLensToy, overwrite = TRUE)
