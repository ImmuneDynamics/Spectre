# Get base image. Tidyverse R studio
FROM rocker/geospatial:4.1.0
LABEL maintainer="Thomas Ashhurst <thomas.ashhurst@sydney.edu.au>"

# Set Locale to utf-8
RUN sed -i -e 's/# en_US.UTF-8 UTF-8/en_US.UTF-8 UTF-8/' /etc/locale.gen && \
    locale-gen
ENV LANG en_US.UTF-8
ENV LANGUAGE en_US:en
ENV LC_ALL en_US.UTF-8

# Install debian specific packages needed by Spectre
RUN apt-get update --fix-missing && \
    apt-get install -y libjpeg-dev && \
    apt-get clean && \
    rm -rf /var/lib/apt/lists/*






# Install all the packages Spectre depends on
RUN install2.r --error  devtools \
			data.table \
                        plyr \
                        dplyr \
                        tidyr \
                        rstudioapi \
                        Rtsne \
                        umap \
                        ggplot2 \
                        ggthemes \
                        scales \
                        colorRamps \
                        RColorBrewer \
                        gridExtra \
                        reticulate \
                        pheatmap \
                        ggpointdensity \
                        ggpubr \
                        viridis \
                        BiocManager \
			gtools \
                        tidygraph \
                        ggraph \
                        stringr \
                        caret \
                        FNN \
                        factoextra \
			raster \
                        tiff \
                        rgeos \
                        tidyr \
                        sp \
                        stars \
                        qs \
                        s2 \
                        exactextractr \
			sf

# Install Spectre dependencies from Bioconductor
RUN R -e "library(BiocManager); BiocManager::install(c('flowCore', 'Biobase', 'flowViz', 'FlowSOM', 'rhdf5', 'HDF5Array'), dependencies=TRUE)"
     
                 
# Install Spectre package. Master branch
RUN R -e "library(devtools); install_github('ImmuneDynamics/Spectre', 'pre-release')"


# Install cytonorm package.
RUN R -e "library(devtools); install_github('saeyslab/CytoNorm')"


# Install miniconda for ChronoClust
ENV LANG=C.UTF-8 LC_ALL=C.UTF-8
ENV PATH /opt/conda/bin:$PATH

RUN apt-get update --fix-missing && \
    apt-get install -y wget bzip2 ca-certificates curl git && \
    apt-get clean && \
    rm -rf /var/lib/apt/lists/*

RUN wget --quiet https://repo.anaconda.com/miniconda/Miniconda3-py37_4.8.2-Linux-x86_64.sh -O ~/miniconda.sh && \
    /bin/bash ~/miniconda.sh -b -p /opt/conda && \
    rm ~/miniconda.sh && \
    /opt/conda/bin/conda clean -tipsy && \
    ln -s /opt/conda/etc/profile.d/conda.sh /etc/profile.d/conda.sh && \
    echo ". /opt/conda/etc/profile.d/conda.sh" >> ~/.bashrc && \
    echo "conda activate base" >> ~/.bashrc

# Install ChronoClust dependencies
RUN conda install -y pandas && \
    conda install -y numpy && \
    conda install -y scipy && \
    conda install -y scikit-learn && \
    conda install -y tqdm && \
    conda install -y numba

# Install ChronoClust
RUN git clone https://github.com/ghar1821/Chronoclust.git chronoclust_repo && \
    cd chronoclust_repo && \
    python3 setup.py install
