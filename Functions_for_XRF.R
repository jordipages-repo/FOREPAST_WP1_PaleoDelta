# # # # # #
# Functions to process XRF data
# Jordi F. Pagès
# 01-10-2019
# CEAB - CSIC 
# # # # # # 

# The coefficient of variation is just the ratio between sd/mean
cv <- function(x,...){
  sd(x)/mean(x)*100
}

# To be able to do "not %in%"
'%!in%' <- function(x,y)!('%in%'(x,y))
