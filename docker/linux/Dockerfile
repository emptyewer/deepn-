FROM ubuntu:20.04

ENV TZ=America/Los_Angeles
RUN ln -snf /usr/share/zoneinfo/$TZ /etc/localtime && echo $TZ > /etc/timezone

# Ubuntu Installs chan
RUN apt-get update
RUN apt-get install -y build-essential gcc axel lzip make openssl p7zip-full patch perl python ruby sed unzip wget xz-utils vim
RUN apt-get install -y autoconf automake autopoint bash bison bzip2 flex g++ g++-multilib gettext git gperf intltool python3-pip
RUN apt-get install -y libc6-dev-i386 libgdk-pixbuf2.0-dev libltdl-dev libssl-dev libtool-bin libxml-parser-perl libclang-dev

# Remove all apt pacakges
RUN apt-get install -y libfontconfig1-dev libfreetype6-dev libx11-dev libx11-xcb-dev libxext-dev libxfixes-dev libxi-dev libxrender-dev libxcb1-dev libxcb-glx0-dev libxcb-keysyms1-dev libxcb-image0-dev libxcb-shm0-dev libxcb-icccm4-dev libxcb-sync-dev libxcb-xfixes0-dev libxcb-shape0-dev libxcb-randr0-dev libxcb-render-util0-dev libxcb-util-dev libxcb-xinerama0-dev libxcb-xkb-dev libxkbcommon-dev libxkbcommon-x11-dev libglu1-mesa-dev
RUN apt-get install -y libharfbuzz-dev libicu-dev libjpeg-dev libtiff-dev libmng-dev libwebp-dev libsqlite-dev libinput-dev ibus postgresql libproxy-dev cups
RUN apt-get install -y llvm libdbus-1-dev libdbus-glib-1-dev libdbus-glib1.0-cil-dev libdbus-glib2.0-cil-dev libdbus1.0-cil-dev libdbus2.0-cil-dev libdbusada5-dev libdbuskit-dev libdbusmenu-glib-dev libdbusmenu-gtk-dev libdbusmenu-gtk3-dev libsqlite3-dev libcppdb-sqlite3-0 libvsqlitepp-dev libwxsqlite3-3.0-dev libzstd-dev libgnatcoll-syslog2-dev libdouble-conversion-dev libglobus-gssapi-gsi-dev  libglobus-gss-assist-dev libglobus-gssapi-error-dev libdrm-dev libgbm-dev libmd4c-dev libxcb-xinput-dev libgulkan-dev libvulkan-dev vulkan-validationlayers-dev libvkd3d-dev libts-dev libhts-dev libcups2-dev libcddb2-dev libmmdb2-dev libmysql++-dev libmysqlclient-dev libmysqlcppconn-dev libsoci-dev unixodbc-dev libpqxx-dev libodb-pgsql-dev freetds-dev libsdl2-dev libsdl1.2-dev libbz2-dev
RUN apt-get clean autoclean && apt-get autoremove --yes && rm -rf /var/lib/{apt,dpkg,cache,log}/

# QT Static Install
RUN mkdir -p qt
WORKDIR /qt
RUN axel -n 10 http://mirrors.ocf.berkeley.edu/qt/archive/qt/5.15/5.15.2/single/qt-everywhere-src-5.15.2.tar.xz
RUN tar xvf qt-everywhere-src-5.15.2.tar.xz
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
