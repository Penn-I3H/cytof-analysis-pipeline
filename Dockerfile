FROM rocker/r-ver:4.3.2

WORKDIR /service

RUN apt clean && apt-get update && apt-get -y install alien

# R program dependencies
RUN apt-get install -y libudunits2-dev && apt-get install -y libgeos-dev && apt-get install -y libproj-dev && apt-get -y install libnlopt-dev && apt-get -y install pkg-config && apt-get -y install gdal-bin && apt-get install -y libgdal-dev
RUN apt-get -y install libcurl4-openssl-dev libfontconfig1-dev libxml2-dev libharfbuzz-dev libfribidi-dev libfreetype6-dev libpng-dev libtiff5-dev libjpeg-dev
RUN apt-get -y install glpk-utils libglpk-dev glpk-doc

RUN R --version

COPY dependencies ./dependencies

## Add additional program specific dependencies below ...
# tidyverse and its prerequisites
RUN Rscript -e "install.packages(c('BH'), Ncpus = 10, repos = 'https://cloud.r-project.org/', dependencies = TRUE)"
RUN Rscript -e "install.packages(c('tidyverse'), Ncpus = 10, dependencies=TRUE)"
RUN Rscript -e "library(tidyverse)" # sanity check

# fix issues with irlba
RUN install.packages("Matrix", type = "source")
RUN install.packages("irlba", type = "source")

# various dependencies
RUN Rscript -e "install.packages(c('KernSmooth', 'patchwork', 'uwot', 'ash', 'RColorBrewer', 'reshape2'), Ncpus = 10, dependencies=TRUE)"
RUN Rscript -e "install.packages(c('Rcpp', 'RcppArmadillo', 'RcppHNSW'), Ncpus = 10, dependencies=TRUE)"

# igraph
RUN Rscript -e "install.packages(c('igraph'), Ncpus = 10, dependencies=TRUE)"

# flowCore and its prerequisites
RUN Rscript -e "install.packages(c('BiocManager'), Ncpus=10)"
RUN Rscript -e "BiocManager::install('RProtoBufLib')"
RUN Rscript -e "BiocManager::install('cytolib', version='3.18')"
RUN Rscript -e "library(cytolib)" # sanity check
RUN Rscript -e "BiocManager::install('flowCore')"

# install CL2
RUN Rscript -e "install.packages('./dependencies/CL2', repos=NULL, type='source')"

# install old architecture dependencies

RUN apt-get -y install wget && apt-get -y install gnupg && apt-get -y install curl
RUN wget -O - https://packages.adoptium.net/artifactory/api/gpg/key/public | apt-key add -
RUN echo "deb https://packages.adoptium.net/artifactory/deb $(awk -F= '/^VERSION_CODENAME/{print$2}' /etc/os-release) main" | tee /etc/apt/sources.list.d/adoptium.list

RUN apt clean && apt-get update
# next flow dependencies
RUN apt-get -y install temurin-17-jdk
# install nextflow
RUN wget -qO- https://get.nextflow.io | bash && chmod +x nextflow && cp ./nextflow /usr/local
RUN apt-get -y install graphviz
ENV PATH="${PATH}:/usr/local/"
# cleanup
RUN rm -f /service/nextflow

# install Go
RUN wget https://go.dev/dl/go1.21.0.linux-amd64.tar.gz
RUN  rm -rf /usr/local/go && tar -C /usr/local -xzf go1.21.0.linux-amd64.tar.gz
ENV PATH="${PATH}:/usr/local/go/bin"
# cleanup
RUN rm -f go1.21.0.linux-amd64.tar.gz

# entrypoint
COPY . ./

RUN ls /service

RUN mkdir -p data

RUN go build -o /service/main main.go

ENTRYPOINT [ "/service/main" ]
