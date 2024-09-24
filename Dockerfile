FROM ubuntu:latest as builder

WORKDIR /build


RUN apt-get update -qq && \
    apt-get install -y zlib1g-dev \
                    git \
                    cmake \
                    build-essential \
                    python3 \
                    gcc-9 \
                    g++-9 \
                    && \
    update-alternatives --install /usr/bin/gcc gcc /usr/bin/gcc-9 9 && \
    update-alternatives --install /usr/bin/g++ g++ /usr/bin/g++-9 9

COPY . /workspace
RUN cd /workspace; rm -rf build; mkdir build; cd build; cmake ..; make -j 8; make install;

# # Cleanup cmake and git
# RUN apt remove -y cmake git && apt autoremove -y

FROM ubuntu:latest

LABEL org.opencontainers.image.authors="maxrossi91@gmail.com"
RUN apt-get update -qq && \
    apt-get install -y zlib1g-dev \
    python3

COPY --from=builder /usr/local/bin /bin
CMD ["/bin/moni"]
