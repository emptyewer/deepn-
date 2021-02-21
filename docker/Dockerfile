FROM ubuntu:20.04

ENV TZ=America/Los_Angeles
RUN ln -snf /usr/share/zoneinfo/$TZ /etc/localtime && echo $TZ > /etc/timezone


RUN apt-get update
RUN apt-get install -y software-properties-common lsb-release
RUN apt-key adv --keyserver keyserver.ubuntu.com --recv-keys 86B72ED9
RUN add-apt-repository "deb [arch=amd64] https://pkg.mxe.cc/repos/apt `lsb_release -sc` main"
RUN apt-get update

RUN apt-get update
RUN apt-get install -y build-essential gcc qtcreator qt5-default
RUN apt-get install -y autoconf automake autopoint bash bison bzip2 flex g++ g++-multilib
RUN apt-get install -y lzip make openssl p7zip-full patch perl python ruby sed unzip wget xz-utils
RUN apt-get install -y gettext git gperf intltool libc6-dev-i386 libgdk-pixbuf2.0-dev libltdl-dev libssl-dev libtool-bin libxml-parser-perl

RUN apt-get install -y mxe-x86-64-w64-mingw32.static-cc mxe-i686-w64-mingw32.static-cc
RUN apt-get install -y mxe-x86-64-w64-mingw32.static-qt5 mxe-x86-64-w64-mingw32.static-armadillo
RUN apt-get install -y mxe-i686-w64-mingw32.static-qt5 mxe-i686-w64-mingw32.static-armadillo
RUN apt-get install -y python3-pip
RUN pip3 install gitlab-auto-release

RUN echo "PATH=/usr/lib/mxe/usr/bin:\$PATH" >> /etc/environment
COPY entrypoint.sh /entrypoint.sh

ENTRYPOINT [ "/bin/bash", "entrypoint.sh" ]

# RUN mkdir /install && cd /install
# RUN git clone https://github.com/mxe/mxe.git
# RUN cd mxe && make MXE_TARGETS='x86_64-w64-mingw32.static' qt5 qtcharts qtimageformats qtmultimedia qtnetworkauth qtwinextras qttools tiff armadillo
