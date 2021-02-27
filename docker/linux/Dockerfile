FROM $CI_REGISTRY/emptyewer/deepn-plus/base:latest

# QT Compile and Install
WORKDIR /qt/qt-everywhere-src-5.15.2
RUN ./configure -prefix /opt/qt5.15.2 -sysconfdir /etc/xdg -confirm-license -opensource  -static -release -platform linux-g++ -openssl-linked -dbus-linked -no-opengl -system-harfbuzz -system-sqlite -nomake examples -no-rpath -skip qtwebengine -skip qt3d -skip qtlocation
RUN make -j6 install
WORKDIR /
RUN rm -rf /qt

# Install conan
RUN pip3 install conan gitlab-auto-release
RUN mkdir -p /root/.conan/profiles
COPY docker/linux/default /root/.conan/profiles/
COPY docker/linux/entrypoint.sh /entrypoint.sh
ENTRYPOINT [ "/bin/bash", "/entrypoint.sh" ]
