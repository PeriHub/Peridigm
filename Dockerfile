FROM registry.gitlab-test.dlr.de/hess_ja/images/trilinos  
MAINTAINER Jan-Timo Hesse <jan-timo.hesse@dlr.de>

ENV HOME /root

WORKDIR /

RUN apt-get -yq update
RUN apt-get -yq install openmpi-bin
RUN apt-get -yq install openssh-server
RUN apt-get -yq install libboost1.55

#Build Peridigm
RUN mkdir Peridigm
WORKDIR /Peridigm
ADD src src
ADD test test
ADD scripts scripts
ADD examples examples 
ADD CMakeLists.txt .
RUN mkdir /Peridigm/build

WORKDIR /Peridigm/build/
RUN cmake \
    -D CMAKE_BUILD_TYPE:STRING=Release \
    -D CMAKE_INSTALL_PREFIX:PATH=/usr/local/peridigm \
    -D CMAKE_CXX_FLAGS:STRING="-O2 -Wall -std=c++14 -pedantic -Wno-long-long -ftrapv -Wno-deprecated" \
    -D TRILINOS_DIR:PATH=/usr/local/trilinos \
    -D USER_LIBRARY_DIRS:PATH=/Peridigm/src/materials/umats \
    -D CMAKE_CXX_COMPILER:STRING="mpicxx" \
    -D PERFORMANCE_TEST_MACHINE:STRING="WSL" \
    -D ALBANY_SEACAS:BOOL=ON \
    ..; \
    make && make install; \
    cd ..
    # rm -rf peridigm

VOLUME /output
WORKDIR /output
ENV LD_LIBRARY_PATH /usr/local/netcdf/lib
ENV PATH /usr/local/peridigm/bin:$PATH
ENV PATH /usr/local/trilinos/bin:$PATH

# Allow mpirun as root, should only be used in container
ENV OMPI_ALLOW_RUN_AS_ROOT 1
ENV OMPI_ALLOW_RUN_AS_ROOT_CONFIRM 1

RUN mkdir /var/run/sshd

RUN  echo 'root:root' | chpasswd
RUN sed -i'' -e's/^#PermitRootLogin prohibit-password$/PermitRootLogin yes/' /etc/ssh/sshd_config \
        && sed -i'' -e's/^#PasswordAuthentication yes$/PasswordAuthentication yes/' /etc/ssh/sshd_config \
        && sed -i'' -e's/^#PermitEmptyPasswords no$/PermitEmptyPasswords yes/' /etc/ssh/sshd_config \
        && sed -i'' -e's/^UsePAM yes/UsePAM no/' /etc/ssh/sshd_config
RUN service ssh start

WORKDIR /app/

EXPOSE 22
CMD    ["/usr/sbin/sshd", "-D"]
