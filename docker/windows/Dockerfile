FROM openscad/mxe-x86_64-gcc:latest


ENV TZ=America/Los_Angeles
RUN ln -snf /usr/share/zoneinfo/$TZ /etc/localtime && echo $TZ > /etc/timezone

RUN apt-get update
RUN DEBIAN_FRONTEND=noninteractive apt-get install -y -qq keyboard-configuration
RUN apt-get install -y -qq gnupg2 lsb-release software-properties-common > /dev/null
RUN apt-key adv --keyserver keyserver.ubuntu.com --recv-keys 86B72ED9
RUN add-apt-repository "deb [arch=amd64] https://pkg.mxe.cc/repos/apt `lsb_release -sc` main" 

RUN apt-get update
RUN apt-get install -y p7zip-full python3-pip vim
RUN apt-get clean autoclean && apt-get autoremove --yes && rm -rf /var/lib/{apt,dpkg,cache,log}/

# MXE installs
RUN make -j8 qt5 armadillo 

# Install conan
RUN pip3 install conan gitlab-auto-release
RUN mkdir -p /root/.conan/profiles
COPY docker/windows/default /root/.conan/profiles/

RUN echo "PATH=/mxe/usr/bin:\$PATH" >> /etc/environment
COPY docker/windows/entrypoint.sh /entrypoint.sh
ENTRYPOINT [ "/bin/bash", "/entrypoint.sh" ]
