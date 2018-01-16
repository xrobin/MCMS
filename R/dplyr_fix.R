## quiets concerns of R CMD check re non-standard evaluation
# Messages like:
# make_unique_modifications: no visible binding for global variable
# 'start'
# Undefined global functions or variables:
#   modifications pair ratio start
# See https://stackoverflow.com/questions/9439256/ for possibly updated fixes:
if(getRversion() >= "2.15.1")  utils::globalVariables(c("modifications", "pair", "ratio", "start"))
