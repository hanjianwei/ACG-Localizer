FROM python:0.0.2
MAINTAINER Jianwei Han <hanjianwei@gmail.com>

RUN apt-get update
RUN apt-get upgrade -y

RUN apt-get install -y cmake libgmm++-dev liblapack-dev libf2c2-dev unzip

ADD http://www.cs.umd.edu/~mount/ANN/Files/1.1.2/ann_1.1.2.tar.gz /

RUN tar -zxf ann_1.1.2.tar.gz

RUN sed -i '158s/double/float/' ann_1.1.2/include/ANN/ANN.h

RUN cd ann_1.1.2/src && make linux-g++ && mv ../include/* /usr/include && mv ../lib/* /usr/lib

ADD http://people.cs.ubc.ca/~mariusm/uploads/FLANN/flann-1.6.11-src.zip /

RUN unzip flann-1.6.11-src.zip

COPY flann_modification/kmeans_index.h flann-1.6.11-src/src/cpp/flann/algorithms/kmeans_index.h

RUN cd flann-1.6.11-src && cmake -D BUILD_MATLAB_BINDINGS=OFF -D BUILD_PYTHON_BINDINGS=OFF . && make && make install

RUN rm -rf ann* flann*

RUN apt-get clean && rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*

WORKDIR /src/github.com/hanjianwei/ACG-Localizer

ENV PATH /src/github.com/hanjianwei/ACG-Localizer/build/bin:$PATH

CMD ["bash"]
