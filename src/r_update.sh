echo "update.packages(repos=\"https://cran.rstudio.com\", ask=FALSE)" | R --no-save
echo "BiocManager::install(ask=FALSE)" | R --no-save
echo "If update failed, try to turn off anaconda environment"
