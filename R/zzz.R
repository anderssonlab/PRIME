# Global Variables
# 
# This file is included to avoid R CMD check NOTES regarding visible binding of global variables.
# The variables listed here are used in dplyr operations and for foreach iterations.
utils::globalVariables(c(
  "score",  # used as a column name in multiple functions
  "i",      # used in foreach loops
  "."       # used in magrittr pipelines (e.g. . %>% ...)
))