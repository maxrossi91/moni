# Utils command

# Build the docker image

```console
docker build --platform linux/amd64 --no-cache -t maxrossi91/moni . 
```

# Pseudo system test
```console
docker run --platform linux/amd64 -v `pwd`/data:/data  -it maxrossi91/moni bash

mkdir -p out
moni build -r data/SARS-CoV2/SARS-CoV2.1k.fa.gz -o out/sars-cov2 -f
moni mems -i out/sars-cov2 -p data/SARS-CoV2/reads.fastq.gz -o out/reads -s
```
