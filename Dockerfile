# Use an official base image (e.g., Debian)
FROM debian:latest

# Set the working directory
WORKDIR /poisson

# Install dependencies (adjust based on your application)
RUN apt-get update && \
    apt-get install -y make libblas-dev liblapacke-dev liblapack-dev libatlas-base-dev gcc gfortran apt-utils

# Copy the local code to the container
COPY . .

# sudo docker build -t poisson .
# sudo docker run -it --rm -v /home/xepho/Bureau/poisson:/poisson poisson