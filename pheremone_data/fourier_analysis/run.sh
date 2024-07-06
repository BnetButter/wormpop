#docker build -t fourier-analysis . > /dev/null
docker run --rm -v $(pwd):/app fourier-analysis