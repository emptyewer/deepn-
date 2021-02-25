# Build Windows Compilers
cp -r conan ./windows
cd windows
docker build --rm --no-cache -t venkykrishna/windows64:latest .
rm -rf conan
cd ..
docker push venkykrishna/windows64:latest


# Build Linux Compilers
cp -r conan ./linux
cd linux
docker build --rm --no-cache -t venkykrishna/linux64:latest .
rm -rf conan
cd ..
docker push venkykrishna/linux64:latest
