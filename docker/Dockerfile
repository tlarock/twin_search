# Linux build environment for twin_search, primarily so that sanitizers
# (ASan/TSan/UBSan) and valgrind are available, which they largely are not
# on Apple silicon macOS.
FROM ubuntu:24.04

ENV DEBIAN_FRONTEND=noninteractive

RUN apt-get update && apt-get install -y --no-install-recommends \
        build-essential \
        g++-14 \
        gcc-14 \
        cmake \
        ninja-build \
        git \
        ca-certificates \
        libboost-all-dev \
        libtbb-dev \
        libgtest-dev \
        valgrind \
        gdb \
        python3 \
        python3-numpy \
        less \
        vim \
    && rm -rf /var/lib/apt/lists/*

RUN update-alternatives --install /usr/bin/gcc gcc /usr/bin/gcc-14 100 \
 && update-alternatives --install /usr/bin/g++ g++ /usr/bin/g++-14 100

WORKDIR /work/twin_search
