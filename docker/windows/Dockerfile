FROM venkykrishna/mxe:latest

# MXE installs
RUN make -j8 qt5 armadillo cmake

# Install conan
RUN pip3 install conan gitlab-auto-release
RUN mkdir -p /root/.conan/profiles
COPY docker/windows/default /root/.conan/profiles/

RUN echo "PATH=/mxe/usr/bin:\$PATH" >> /etc/environment
COPY docker/windows/entrypoint.sh /entrypoint.sh
ENTRYPOINT [ "/bin/bash", "/entrypoint.sh" ]
