docker build -t samtools-bcftools:bash .
docker run --rm samtools-bcftools:bash bash -lc "bcftools --version; samtools --version"
