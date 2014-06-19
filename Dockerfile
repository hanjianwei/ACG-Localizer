FROM hanjianwei/python:2.7

RUN apt-get update

RUN apt-get install -y cmake libgmm++-dev liblapack-dev libf2c2-dev

RUN apt-get install -y wget unzip

RUN wget http://www.cs.umd.edu/~mount/ANN/Files/1.1.2/ann_1.1.2.tar.gz

RUN tar -zxf ann_1.1.2.tar.gz

RUN sed -i '158s/double/float/' /ann_1.1.2/include/ANN/ANN.h

RUN cd ann_1.1.2/src && make linux-g++ && mv ../include/* /usr/include && mv ../lib/* /usr/lib

RUN wget http://people.cs.ubc.ca/~mariusm/uploads/FLANN/flann-1.6.11-src.zip

RUN unzip flann-1.6.11-src.zip

RUN cd flann-1.6.11-src && cmake -D BUILD_MATLAB_BINDINGS=OFF -D BUILD_PYTHON_BINDINGS=OFF . && make && make install

COPY flann_modification/kmeans_index.h /usr/local/include/flann/algorithms/kmeans_index.h

RUN rm -rf ann* flann*

CMD ["bash"]
