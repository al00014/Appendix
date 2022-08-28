library(tidyverse)
library(BAPC)
library(parallel)
#Read data
setwd("./original data/pop/")
population_list <- list.files('./original data/pop/')
population <- lapply(population_list[str_detect(population_list, "CSV")], 
                     function(x) read_csv(x))

setwd("./original data/DALY_num/")
deathnumber_list <- list.files('./original data/DALY_num/')
deathnumber <- lapply(deathnumber_list[str_detect(deathnumber_list, "csv")], 
                      function(x) read_csv(x))

population_combine <- population[[1]]
deathnumber_combine <- deathnumber[[1]]
for(i in 2:7){
  population_combine <- rbind(population_combine, population[[i]])
  deathnumber_combine <- rbind(deathnumber_combine, deathnumber[[i]])
}

setwd("./original data/")

#Data Preparation
base_table <- expand.grid(c(1990,1995,2000,2005,2010,2015,2019),
                          c(33:41,43:55,57:63,71:72,66:69,101:102,97:99,75:95,
                            121:123,105:109,111:119,385,422,125:133,135:136,160,139:157,522,
                            161:165,6:8,351,22:23,25:30,10:20,183,186,
                            168:173,175:182,184:185,187,189:191,435,193:198,200:218,
                            1,44575:44578),
                            c(294), 3, c(28,5:20,30,160))
names(base_table) <- c("Year","Area","Cause","Sex","Age")

dat <- left_join(base_table, population_combine, by=c("Year"="year_id","Area"="location_id",
                                                      "Sex"="sex_id","Age"="age_group_id"))[,-c(6:12,14:15)]
names(dat)[6] <- c("Population")

dat <- left_join(dat, deathnumber_combine, by=c("Year"="year","Area"="location_id","Cause"="cause_id",
                                                "Sex"="sex_id","Age"="age_id"))[,-c(7:14,16:17)]
names(dat)[7] <- c("DeathNumber")

dat <- dat %>%
  mutate(Population = ceiling(ifelse(Population <= 1, 1, Population)),
         DeathNumber = ceiling(ifelse(DeathNumber <= 1, 1, DeathNumber / 10000)))

dat_list <- dat %>% split(list(.$Year, .$Area, .$Cause, .$Sex))

age_combine <- function(dat){
  dat[2,6] <- dat[2,6] + dat[1,6]
  dat[1,6] <- NA
  dat[2,7] <- dat[2,7] + dat[1,7]
  dat[1,7] <- NA
  dat <- dat[-1,]
}
dat_list <- lapply(dat_list, age_combine)

Population_list <- list()
DeathNumber_list <- list()
for (areaNum in c(33:41,43:55,57:63,71:72,66:69,101:102,97:99,75:95,
                  121:123,105:109,111:119,385,422,125:133,135:136,160,139:157,522,
                  161:165,6:8,351,22:23,25:30,10:20,183,186,
                  168:173,175:182,184:185,187,189:191,435,193:198,200:218,
                  1,44575:44578)){
  for (causeNum in c(294)){
    for (sexNum in c(3)){
      data_PopulationAPC <- matrix(nrow=0, ncol=18)
      data_DeathNumberAPC <- matrix(nrow=0, ncol=18)
      for (yearNum in c(1990,1995,2000,2005,2010,2015,2019)){
        PopulationAPC_name <- paste(yearNum, areaNum, causeNum, sexNum, sep = ".")
        PopulationAPC_00 <- t(matrix(dat_list[[PopulationAPC_name]][,"Population"]))
        data_PopulationAPC <- rbind(data_PopulationAPC, PopulationAPC_00)
        
        DeathNumberAPC_name <- paste(yearNum, areaNum, causeNum, sexNum, sep = ".")
        DeathNumberAPC_00 <- t(matrix(dat_list[[DeathNumberAPC_name]][,"DeathNumber"]))
        data_DeathNumberAPC <- rbind(data_DeathNumberAPC, DeathNumberAPC_00)
      }
      PopulationAPC_name_1 <- paste(areaNum, causeNum, sexNum, sep = ".")
      Population_list[[PopulationAPC_name_1]]<-data_PopulationAPC
      
      DeathNumberAPC_name_1 <- paste(areaNum, causeNum, sexNum, sep = ".")
      DeathNumber_list[[DeathNumberAPC_name_1]] <- data_DeathNumberAPC
    }}}

Population_list_APC <- lapply(Population_list, function(x){
  library(tidyverse)
  NA_matrix <- matrix(0, nrow = 6, ncol = 18)
  rbind(x, NA_matrix) %>% as.data.frame(x)})

DeathNumber_list_APC <- lapply(DeathNumber_list, function(x){
  library(tidyverse)
  NA_matrix <- matrix(NA, nrow = 6, ncol = 18)
  rbind(x, NA_matrix) %>% as.data.frame(x)})

APC_object <- list()
for(i in c(names(Population_list_APC))){
  APC_object[[i]] <- APCList(DeathNumber_list_APC[[i]], Population_list_APC[[i]], gf = 5)
}

#The BAPC model
function_BAPC <- function(x){
  library(BAPC)
  Population_std <- c(sum(131580709.2,500535402.8),585162350.4,536689497.6,519603412.1,492675757.3,442843296.3,
                      385628936.3,352742337.2,286291792.8,232409454.8,212591110,185415509.9,160643040.8,
                      123488561,84513452,61309223,35219664.5,20504179.2) #Population estimates from GBD 1990
  BAPC(x, predict = list(npredict = 6, retro = FALSE),
       model = list(age = list(model = "rw2", prior = "loggamma", param = c(1, 0.00005), initial = 4, scale.model = FALSE),
                    period = list(include = TRUE, model = "rw2", prior = "loggamma", param = c(1, 0.00005), initial = 4, scale.model = FALSE),
                    cohort = list(include = TRUE, model = "rw2", prior = "loggamma", param = c(1, 0.00005), initial = 4, scale.model = FALSE),
                    overdis = list(include = TRUE, model = "iid", prior = "loggamma", param = c(1, 0.005), initial = 4)),
       secondDiff = FALSE, stdweight = Population_std, verbose = FALSE)
}

cl.cores <- detectCores()
cl <- makeCluster(cl.cores)

apc_model_list <- parLapply(cl, APC_object, function_BAPC)

stopCluster(cl)

apc_model_list <- lapply(apc_model_list, function(x) qapc(x, percentiles=c(0.025, 0.975)))

save(apc_model_list, file = "./intermediate result/BAPC_result.Rdata")
