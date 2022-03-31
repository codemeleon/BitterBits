# https://blog.sellorm.com/2017/10/21/quick-script-to-install-an-r-package-from-the-command-line/
chmod u+x rpkginstall
cp ./rpkginstall ~/bin/

chmod u+x rbiocinstall
cp ./rbiocinstall ~/bin/

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
rpkginstall ggnewscale
rpkginstall phylocanvas
rpkginstall tidymodels
rpkginstall modeltime
rpkginstall timetk
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
rpkginstall micropan
rpkginstall ggdendro
rpkginstall plotly
rpkginstall readxl

# Some ggplot2 packages
rpkginstall ggplot2
rpkginstall ggthemes
rpkginstall svglite
rpkginstall TDbook

# Geolocation
rpkginstall OpenStreetMap
rpkginstall DT
rpkginstall RColorBrewer
rpkginstall mapproj
rpkginstall sf
rpkginstall RgoogleMaps
rpkginstall scales
rpkginstall rworldmap
rpkginstall maps
rpkginstall tidyverse
rpkginstall rnaturalearth
rpkginstall rnaturalearthdata
rpkginstall rgeos
rpkginstall ggspatial
rpkginstall maptools
rpkginstall leaflet
rpkginstall sf
rpkginstall tmap
rpkginstall here
rpkginstall rgdal
rpkginstall scales
rpkginstall flextable

# install package from github
devtools::install_github("dkahle/ggmap", ref = "tidyup")
# install klippy for copy-to-clipboard button in code chunks
remotes::install_github("rlesur/klippy")


# Biocundoctor packages

rbiocinstall msa

# remotes::install_github("ThinkR-open/remedy")
# install.packages("colourpicker")
# install.packages("ggThemeAssist")
# devtools::install_github("dracodoc/mischelper")
# devtools::install_github("tylermorganwall/rayshader")

remotes::install_github("business-science/modeltime", dependencies = TRUE)
