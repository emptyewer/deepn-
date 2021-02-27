FROM emptyewer/deepn-plus/base:latest

# QT Compile and Install
RUN make -j6 install
WORKDIR /
RUN rm -rf /qt

# Install conan
RUN pip3 install conan gitlab-auto-release
RUN mkdir -p /root/.conan/profiles
COPY docker/linux/default /root/.conan/profiles/
COPY docker/linux/entrypoint.sh /entrypoint.sh
ENTRYPOINT [ "/bin/bash", "/entrypoint.sh" ]
