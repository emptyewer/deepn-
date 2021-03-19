FROM registry.gitlab.com/emptyewer/deepn-plus/base:latest

# QT Compile and Install
RUN make -j8 -s install

WORKDIR /
RUN rm -rf /qt

RUN useradd -m gitlab-runner

# Install conan 
RUN pip3 install conan gitlab-auto-release
RUN mkdir -p /home/gitlab-runner/.conan/profiles
COPY docker/windows/default /home/gitlab-runner/.conan/profiles/

RUN chown -R gitlab-runner:gitlab-runner /home/gitlab-runner

COPY docker/linux/entrypoint.sh /entrypoint.sh
ENTRYPOINT [ "/bin/bash", "/entrypoint.sh" ]