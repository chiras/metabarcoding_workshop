
# VSEARCH
https://github.com/torognes/vsearch/releases/tag/v2.23.0
Download & unzip respective binary (bin/vsearch) and save in this projects ```./bin``` folder

# SeqFilter

```sh
git clone https://github.com/BioInf-Wuerzburg/SeqFilter
cd SeqFilter/
make
cp bin/SeqFilter ../metabarcoding_workshop/bin/
cp -R lib ../metabarcoding_workshop/
```
# R/Rstudio
https://cran.rstudio.com
https://posit.co/download/rstudio-desktop/

# R packages vegan, phyloseq, bipartite,ggplot2

```R
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("phyloseq")

install.packages("vegan","bipartite","ggplot2")
```
