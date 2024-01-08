# Use an official base image (e.g., Debian)
FROM debian:10

# Set the working directory
WORKDIR /poisson

# Install dependencies (adjust based on your application)
RUN apt-get update && \
    apt-get install -y libblas-dev liblapacke-dev gcc gfortran make apt-utils fish bash

# Copy the local code to the container
COPY . .

# sudo docker build -t poisson .
# sudo docker run -it --rm -v /home/xepho/Bureau/poisson:/poisson poisson