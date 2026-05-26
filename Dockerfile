FROM rocker/r-ver:4.3.0

RUN apt-get update && apt-get install -y \
    git \
    gfortran \
    libcurl4-openssl-dev \
    libssl-dev \
    libxml2-dev \
    libgit2-dev \
    libhdf5-dev \
    libpng-dev \
    libjpeg-dev \
    libtiff5-dev \
    libfontconfig1-dev \
    libfreetype6-dev \
    libharfbuzz-dev \
    libfribidi-dev \
    libmagick++-dev \
    pandoc \
    make \
    && rm -rf /var/lib/apt/lists/*

WORKDIR /project

COPY renv.lock renv.lock
COPY .Rprofile .Rprofile
COPY renv/activate.R renv/activate.R
COPY renv/settings.json renv/settings.json

RUN R -e "install.packages(c('renv', 'BiocManager')); BiocManager::install(version = '3.17', ask = FALSE); renv::restore(prompt = FALSE)"

COPY . .

CMD ["R", "-e", "targets::tar_make()"]