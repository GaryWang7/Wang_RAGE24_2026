# RHOME=/home/rstudio
# docker run -it \
# --cpus=14 \
# --memory="600g" \
# --workdir $HOME \
# --name ragegex \
# -v /mnt/y:$HOME/y \
# -v /mnt/j:$HOME/j \
# -v /mnt/g:$HOME/g \
# -v /mnt/e:$HOME/e \
# -v /mnt/d/garyw:$HOME/garyw \
# -v /mnt/v:$HOME/v \
# -v $HOME:$HOME \
# -v /mnt/y:$RHOME/y \
# -v /mnt/j:$RHOME/j \
# -v /mnt/g/:$RHOME/g \
# -v /mnt/e:$RHOME/e \
# -v /mnt/d/garyw:$RHOME/garyw \
# -v /mnt/v:$RHOME/v \
# -v /var/run/docker.sock:/var/run/docker.sock \
# -e PASSWORD=garywang \
# -e DISABLE_AUTH=TRUE \
# -p 8787:8787 \
# garywang7/ragemultiome:1.0.4

# install.packages("ggpubr")
# install.packages("svglite")
# install.packages("msigdbr")
# install.packages("logger")
library(svglite)
library(tidyverse)
library(here)
library(future)
library(fs)
library(Seurat)
library(scCustomize)
library(clusterProfiler)
library(fgsea)


