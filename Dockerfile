# docker build -t reviewer:0.2.7 .
# docker run reviewer:0.2.7 REViewer --help
# work from latest LTS ubuntu release
FROM debian:buster

# run update and install necessary tools
RUN apt-get update -y && apt-get install -y \
    libcurl4-openssl-dev \
    g++ \
    libzmq3-dev \
    zlib1g-dev \
    libbz2-dev \
    liblzma-dev \
    cmake \
    git \
    wget \
    unzip \
    && rm -rf /var/lib/apt/lists/*

WORKDIR /usr/bin
# install REViewer
RUN wget https://github.com/Illumina/REViewer/archive/refs/tags/v0.2.7.zip; \
    unzip v0.2.7.zip; \
    cd REViewer-0.2.7/; \
    mkdir build; cd build; \
    cmake ..; \
    make; \
    ln -s /usr/bin/REViewer-0.2.7/build/install/bin/REViewer /usr/bin/REViewer
CMD ["bash"]
