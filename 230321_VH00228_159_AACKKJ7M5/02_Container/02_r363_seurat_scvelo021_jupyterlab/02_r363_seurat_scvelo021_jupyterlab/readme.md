
# Jupyter lab

This dockerfile imports the docker image of interest (FROM) and adds jupyter lab
 to it for ineractive tests.

This version is added to a R-based image so it also adds:
 - IRkenel for running R sessions in jupyter
 - rpy2 for running R code from python session

## Build

Here, jupyter layer is added to existing 'r3.6.3_seurat_scvelo' Dockerfile
docker build . -t rfenouil/r363_seurat_scvelo021_jupyterlab

## Save

docker save rfenouil/r363_seurat_scvelo021_jupyterlab | gzip > r363_seurat_scvelo021_jupyterlab.tar.gz

## Run

Must forward a port for browser connection, and a token used for initial login:

```
docker run -d -u $(id -u ${USER}):$(id -g ${USER}) \
           -p 8888:8888 \
           -e TOKEN=myPass \
           -v /mnt:/mnt \
           rfenouil/r363_seurat_scvelo021_jupyterlab
```

Then a browser must be used to connect locally: `https://127.0.0.1:8888`
