FROM ubuntu:jammy

RUN apt-get update && \
	apt-get install -y build-essential git cmake autoconf libtool pkg-config libboost-all-dev libnetcdf-c++4-dev subversion nvidia-cuda-dev nvidia-cuda-toolkit
	
RUN svn checkout https://mmoorkamp@svn.code.sf.net/p/jif3d/jif3dsvn/trunk/jif3D .
RUN rm boost-build.jam
RUN b2 release jointinv
CMD ["Joint/bin/gcc-11/release/parallelization-openmp/threading-multi/jointinv"]

