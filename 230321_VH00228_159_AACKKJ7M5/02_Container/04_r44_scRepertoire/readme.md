
# r44_screpertoire

Based on rocker/tidyverse distribution (includes Rstudio): 
Added:
 - Misc packages for data manipulation, figures, reports
 - Seurat
 - scRepertoire for clonotype analyses



## Build

docker build -t r44_screpertoire .



## Save

docker save r44_screpertoire > r44_screpertoire.tar.gz



## Run Rstudio

Give user details to internal script which sets user and permissions:

```
docker run -d --name r42_sccomp -p 9898:8787 -e PASSWORD=yourPass -e USER=$(whoami) -e USERID=$(id -u) -e GROUPID=$(id -g) -v /mnt:/mnt r42_sccomp_shazam
```

Then connect to the machine running docker (localhost) on mapped port (8787):
http://127.0.0.1:9898


