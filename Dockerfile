FROM intel/oneapi-hpckit

RUN apt update -y; apt upgrade -y; exit 0
RUN apt install gdb -y

RUN curl -LO https://github.com/ClementTsang/bottom/releases/download/0.9.6/bottom_0.9.6_amd64.deb
RUN dpkg -i bottom_0.9.6_amd64.deb

WORKDIR /
RUN curl -LO https://ftp.gnu.org/gnu/make/make-4.4.tar.gz
RUN tar -xvf make-4.4.tar.gz
WORKDIR /make-4.4
RUN ./configure
RUN make 
RUN make install

WORKDIR /app/
ENTRYPOINT ["/bin/bash"]
