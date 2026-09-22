
# r42_sccomp_shazam

Based on rocker/tidyverse distribution (includes Rstudio): 
Added:
 - Misc packages for data manipulation, figures, reports
 - Seurat
 - sccomp for cell identity proportion comparison
 - shazam for hypermutation analysis



## Build

docker build -t r42_sccomp_shazam .



## Save

docker save r42_sccomp_shazam > r42_sccomp_shazam.tar.gz



## Run Rstudio

Give user details to internal script which sets user and permissions:

```
docker run -d --name r42_sccomp -p 9898:8787 -e PASSWORD=yourPass -e USER=$(whoami) -e USERID=$(id -u) -e GROUPID=$(id -g) -v /mnt:/mnt r42_sccomp_shazam
```

Then connect to the machine running docker (localhost) on mapped port (8787):
http://127.0.0.1:9898


