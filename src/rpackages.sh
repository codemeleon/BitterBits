# https://blog.sellorm.com/2017/10/21/quick-script-to-install-an-r-package-from-the-command-line/
chmod +x rpkginstall

cp ./rpkginstall ~/bin/


rpkginstall dplyr
rpkginstall ggplot2
rpkginstall ggtree
rpkginstall devtools
rpkginstall tidyverse
rpkginstall esquisse # https://github.com/dreamRs/esquisse
rpkginstall ggThemeAssist
rpkginstall ggedit
rpkginstall reticulate
rpkginstall plumber
rpkginstall flexdashboard
rpkginstall purrr
rpkginstall shiny
rpkginstall miniUI
rpkginstall phylotools
rpkginstall ggsignif
rpkginstall gggenes

# For Metagenomics course
rpkginstall phyloseq
rpkginstall vegan
rpkginstall tidyverse
rpkginstall here
rpkginstall dada2
rpkginstall decontam
rpkginstall ggforce
rpkginstall ggtext
rpkginstall magrittr
rpkginstall ggrepel
rpkginstall RColorBrewer
rpkginstall scales
rpkginstall Rmarkdown
rpkginstall metagenomeSeq
rpkginstall fpc
rpkginstall ggdendro
rpkginstall dendextend
rpkginstall cowplot



# remotes::install_github("ThinkR-open/remedy")
# install.packages("colourpicker")
# install.packages("ggThemeAssist")
# devtools::install_github("dracodoc/mischelper")
# devtools::install_github("tylermorganwall/rayshader")
