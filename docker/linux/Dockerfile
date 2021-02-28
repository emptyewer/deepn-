FROM registry.gitlab.com/emptyewer/deepn-plus/linux64-stage3:latest
WORKDIR /
RUN rm -rf /qt

# Install conan
RUN pip3 install conan gitlab-auto-release
RUN mkdir -p /root/.conan/profiles
COPY docker/linux/default /root/.conan/profiles/
COPY docker/linux/entrypoint.sh /entrypoint.sh
ENTRYPOINT [ "/bin/bash", "/entrypoint.sh" ]
