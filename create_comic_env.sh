#!/bin/bash
# Clone the GitHub repository

git clone https://github.com/nis2034/comics.git comics-test

# Pull the docker image from KS GitHub repository

docker.exe pull ghcr.io/karstensuhre/tensordocker:2.0

# Run the container if not already
if ! docker.exe ps -a --format '{{.Names}}' | grep tensor2_test; then

	docker.exe run -v D:\\nis2034\\metabolomics\\comics-test\\comics-test:/home/rstudio/host -v D:\\nis2034\\python:/tf/nis2034 -it --detach --name tensor2_test.0 -p8881:8888 -p8784:8787 ghcr.io/karstensuhre/tensordocker:2.0
else
        echo "Container tensor2_test already exists. Skipping container creation."
fi

# Run rstudio inside the container

docker.exe exec -it tensor2_test.0 rstudio-server start

# Install shinydashboard package

docker.exe exec -it tensor2_test.0 R -e "install.packages('shinydashboard',repos='https://cran.rstudio.com/')"

