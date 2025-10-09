## creating a table to describe the continuous family distribution
library(ggpubr)

famdist <- data.frame(
  family = c("NO", "NO2", "GU", "LO", "RG",
             "exGAUS","NOF", "PE", "PE2", "SN1", "SN2","TF", "TF2", 
             "EGB2","GT", "JSU", "JSUo", "NET", "SHASH", "SHASHo", "SHASHo2",
             "SEP1","SEP2","SEP3","SEP4", "SST","ST1", "ST2", "ST3", "ST4", "ST5"),
  param = c("2", "2", "2", "2", "2", 
            "3","3","3","3","3","3","3","3",
            "4","4","4","2","4","4","4","4","4","4","4","4","4","4","4","4", "4", "4"),
  skewness = c("symmetrical", "symmetrical", "negativly skewed", "symmetrical", "postively skewed", 
               "postively skewed","symmetrical", "symmetrical", "symmetrical", "both","both", "symmetrical", "symmetrical",
               "both","symmetrical", "both", "both", "symmetrical", "both", "both", "both",
               "both","both","both","both", "both", "both", "both", "both", "both", "both"),
  kurtosis = c("mesokurtic", "mesokurtic", "mesokurtic", "leptokurtic", "leptokurtic",
               "leptokurtic","mesokurtic", "both", "both", "mesokurtic","mesokurtic", "leptokurtic", "leptokurtic",
               "leptokurtic","both", "leptokurtic", "leptokurtic", "leptokurtic", "both", "both", "both",
               "both","both","both","both","leptokurtic" ,"leptokurtic", "leptokurtic", "leptokurtic", "leptokurtic", "leptokurtic")
)
# change column names 
colnames(famdist) <- c("Continous family distribution", "No. of parameters", "Skewness", "Kurtosis")

# make them table 
famdist_table <- ggtexttable(famdist, rows = NULL)


pdf("famdist_table.pdf", width = 6, height = 8.75)
famdist_table
dev.off()


## creating a table to show all sims and the best continuous family distribution that describes the sim

sim_results <- data.frame(
  sim = c("1", "2", "3", "4", "5", "6", "7", "8"),
  sim_skew = c("symmetrical", "positively skewed", "symmetrical", "positively skewed",
               "both", "both", "positively skewed", "both"),
  sim_kurt = c("mesokurtic", "mesokurtic", "leptokurtic", "leptokurtic",
               "mesokurtic", "leptokurtic", "both", "both"),
  best_model = c("NO", "SN2", "SST", "ST3",
                 "SHASHo", "JSUo", "JSUo", "SHASH"),
  model_skew = c("symmetrical", "both", "both", "both",
                 "both", "both", "both", "both"),
  model_kurt = c("mesokurtic", "mesokurtic", "leptokurtic", "leptokurtic",
                 "both", "leptokurtic", "leptokurtic", "both")
)
# change column names 
colnames(sim_results) <- c("simulations", "simulation skewness", "simulated kurtosis", "best family", 
                       "best family skewness", "best family kurtosis")

# make them table 
sim_results_table <- ggtexttable(sim_results, rows = NULL, theme = ttheme("light", base_size = 12))

pdf("sim_result.pdf", width = 10, height = 2)
sim_results_table
dev.off()

# another one 
sim_results2 <- data.frame(
  sim = c("1", "2", "3", "4"),
  sim_skew = c("symmetrical", "positively skewed", "symmetrical", "positively skewed"),
  sim_kurt = c("mesokurtic", "mesokurtic", "leptokurtic", "leptokurtic"),
  best_model = c("NO", "SN2", "SST", "ST3"),
  model_skew = c("symmetrical", "both", "both", "both"),
  model_kurt = c("mesokurtic", "mesokurtic", "leptokurtic", "leptokurtic")
)
# change column names 
colnames(sim_results) <- c("simulations", "simulated skewness", "simulated kurtosis", "best family", 
                           "best family skewness", "best family kurtosis")
# make them table 
sim_results_table2 <- ggtexttable(sim_results2, rows = NULL, theme = ttheme("light", base_size = 12))



# another one with changing parameters 
sim_change_results <- data.frame(
  sim = c("5", "6", "7"),
  sim_skew = c("negative to positive", "positively skewed", "negative to positive"),
  sim_kurt = c("leptokurtic", "mesokurtic to leptokurtic", "mesokurtic to leptokurtic"),
  best_model = c("SST", "SEP2", "ST3"),
  model_skew = c("both", "both", "both"),
  model_kurt = c("leptokurtic", "both", "leptokurtic")
)
# change column names 
colnames(sim_change_results) <- c("simulations", "simulated skewness", "simulated kurtosis", "best family", 
                           "best family skewness", "best family kurtosis")
# make them table 
sim_change_results_table <- ggtexttable(sim_change_results, rows = NULL, theme = ttheme("light", base_size = 12))

pdf("more_sim_result.pdf", width = 10, height = 2)
sim_change_results_table
dev.off()


