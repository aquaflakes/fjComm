largeListToDf <- function(list)
{
  n <- length(list[[1]])
  structure(list, row.names = c(NA, -n), class = "data.frame")
}

letter2num <- function(x) { sapply(x, function(x) {utf8ToInt(x) - utf8ToInt("A") + 1L}) }
num2letter <- function(x) {  sapply(x, function(x) {intToUtf8 (x + utf8ToInt("A") - 1L)})  }

CN_sources <- function(){
  options("repos"=c(CRAN="https://mirrors.tuna.tsinghua.edu.cn/CRAN/"))
  options(BioC_mirror="https://mirrors.tuna.tsinghua.edu.cn/bioconductor")
}
