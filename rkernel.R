install.packages('RCurl')
library(devtools)
git clone https://github.com/armstrtw/rzmq.git --recursive
install_local('./rzmq')
install_github('IRkernel/repr')
install_github('IRkernel/IRdisplay')
install_github('IRkernel/IRkernel')
IRkernel::installspec()
