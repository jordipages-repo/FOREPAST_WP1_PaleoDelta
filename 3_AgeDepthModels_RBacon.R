# # # # # #
# Age depth models 
# Jordi F. Pagès
# 03-03-2020
# Universitat de Barcelona
# CEAB-CSIC
# # # # # # 

# Libraries
library(rbacon)


# # # # # # # # # # # # # # # # # #
# 1) 1st step after receiving report from DirectAMS: Calculate calibrated age BP from pMC -----
# # # # # # # # # # # # # # # # # #
# cores <- read.csv(file = "Data/RadiocarbonData/PMC_FinalReport_2020Feb27.csv")
# myAge <- pMC.age(mn = cores$pMC, sdev = cores$pMC_1sigma_error)
# cores$myAge <- myAge[1:15]
# cores$myAge_error <- myAge[16:30]
# # write_csv(cores, path = "Data/RadiocarbonData/FinalReport_2020Feb27_Ages_calculated_fromPMCbacon.csv")


# # # # # # # # # # # # # # # # # #
# 2) Run age depht models with rbacon -----
# # # # # # # # # # # # # # # # # #

# # # #
# ERMS TANCADA CORE
# # # #
Bacon(core = "ErmsTancada", coredir = "Data/RadiocarbonData/", acc.shape = 1.5, acc.mean = 50, postbomb = 1)

pdf(file = "Figures/RadiocarbonFigs/ErmsTancada_free_surface.pdf")
agedepth(model.only = T)
dev.off()


# # # #
# LLANADA CORE
# # # #
Bacon(core = "Llanada", coredir = "Data/RadiocarbonData/", acc.shape = 1.5, acc.mean = 50, postbomb = 1)

# pdf(file = "Figures/RadiocarbonFigs/Llanada_fixed_surface.pdf")
agedepth(model.only = T)
# dev.off()





