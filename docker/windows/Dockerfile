FROM venkykrishna/mxe:latest

# MXE installs 
RUN make -j8 qt5 armadillo cmake

RUN useradd -m gitlab-runner

# Install conan 
RUN pip3 install conan gitlab-auto-release
RUN mkdir -p /home/gitlab-runner/.conan/profiles
COPY docker/windows/default /home/gitlab-runner/.conan/profiles/

RUN chown -R gitlab-runner:gitlab-runner /home/gitlab-runner

RUN echo "PATH=/mxe/usr/bin:\$PATH" >> /etc/environment
COPY docker/windows/entrypoint.sh /entrypoint.sh
ENTRYPOINT [ "/bin/bash", "/entrypoint.sh" ]
