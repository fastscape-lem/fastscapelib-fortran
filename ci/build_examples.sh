mkdir build
cd build
cmake -DBUILD_FASTSCAPELIB_SHARED=ON -DBUILD_EXAMPLES=ON -DCMAKE_BUILD_TYPE=Debug ..
make