# Linux build environment for twin_search.
#
# Two jobs. First, sanitizers (ASan/TSan/UBSan) and valgrind, which are largely
# unavailable on Apple silicon macOS. Second - and this is what makes it worth
# keeping up to date - it is a stand-in for della's SAMPLING. std::uniform_int_
# distribution is implementation-defined, so a --seed run reproduces only
# within one standard library; macOS libc++ and della's libstdc++ draw
# different hypergraphs from the same mt19937 stream. g++-14 here gives the
# same libstdc++ as della's gcc-toolset-14, and this was verified against real
# della output: dry-run fingerprints for n=9 m=47 k=3 seed 20260930 accounted
# for all 453 rows of the della-produced CSV, with no unmatchable row and no
# surplus projection. Dry runs and sampling checks can therefore be done here
# rather than queued on the cluster.
#
# `time` is GNU /usr/bin/time, which the slurm-scripts require for peak RSS;
# the shell builtin is not a substitute (it has no -v).
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
        time \
        python3 \
        python3-numpy \
        less \
        vim \
    && rm -rf /var/lib/apt/lists/*

RUN update-alternatives --install /usr/bin/gcc gcc /usr/bin/gcc-14 100 \
 && update-alternatives --install /usr/bin/g++ g++ /usr/bin/g++-14 100

WORKDIR /work/twin_search
