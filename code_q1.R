# code.R
rm(list = ls())
source('variational_functions.R')

y = read.delim('../homework 3/hwk3-data.txt', header = FALSE)[,1]
pri = list(
  mu = -0.5,
  s2 =  4.0,
  a_phi = 5,
  b_phi = 4
  )





# EOF