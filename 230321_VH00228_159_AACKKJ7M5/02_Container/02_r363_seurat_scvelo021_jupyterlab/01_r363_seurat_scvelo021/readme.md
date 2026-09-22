
# r363_seurat_scvelo021

Based on rocker/r-ver:3.6.3 (no Rstudio) distribution: https://github.com/rocker-org/rocker-versioned
Add:
 - R packages for plotting and reporting 
 - Seurat
 - scVelo (python3)
 - Velocyto cli (creation of loom files)
 - Seurat tools for interactions with Python and scVelo (seurat-wrapper, seurat-disk, velocyto.R)

## Build

docker build . -t rfenouil/r363_seurat_scvelo021

## Run

```
docker run -u $(id -u ${USER}):$(id -g ${USER}) \
           -v /mnt:/mnt \
           rfenouil/r363_seurat_scvelo021
```
