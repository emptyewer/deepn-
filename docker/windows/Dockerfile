FROM ubuntu:20.04

ENV TZ=America/Los_Angeles
RUN ln -snf /usr/share/zoneinfo/$TZ /etc/localtime && echo $TZ > /etc/timezone

RUN apt-get update
RUN apt-get install -y gnupg2 lsb-release software-properties-common
RUN apt-key adv --keyserver keyserver.ubuntu.com --recv-keys 86B72ED9
RUN add-apt-repository "deb [arch=amd64] https://pkg.mxe.cc/repos/apt `lsb_release -sc` main"

# MXE Installs
RUN apt-get update
RUN apt-get install -y mxe-x86-64-w64-mingw32.static-cc mxe-x86-64-w64-mingw32.static-qt5 mxe-x86-64-w64-mingw32.static-armadillo p7zip-full python3-pip vim

RUN pip3 install gitlab-auto-release
# Remove all apt pacakges
RUN apt-get clean autoclean && apt-get autoremove --yes && rm -rf /var/lib/{apt,dpkg,cache,log}/

# Install conan
RUN pip3 install conan
RUN mkdir -p /root/.conan/profiles
COPY docker/windows/default /root/.conan/profiles/

RUN echo "PATH=/usr/lib/mxe/usr/bin:\$PATH" >> /etc/environment
COPY docker/windows/entrypoint.sh /entrypoint.sh 
ENTRYPOINT [ "/bin/bash", "/entrypoint.sh" ]
