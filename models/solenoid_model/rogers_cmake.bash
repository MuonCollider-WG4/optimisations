set -e

rm -r build
mkdir build
cmake CMakeLists.txt \
      -Bbuild/ -H./ \
      -DVERBOSE=1 \
      -DPython_LIBRARY=/usr/lib/x86_64-linux-gnu/libpython3.10.so \
      -DPython_INCLUDE_DIR=/usr/include/python3.10/ \
      -DCMAKE_INSTALL_PREFIX=${OPAL_THIRD_PARTY} \
      -DPYTHON_INSTALL_DIR=${OPAL_THIRD_PARTY}/lib/python/site-packages/
cd build
make -j20

# comments:-
# I had to apply manually the patch as in https://github.com/boostorg/phoenix/issues/111
# there was a version in script 110-build-boost but it didn't appear to work
#
# I had to build manually the libgtest.so and libgtest_main.so as shared libraries - gtest
# cmake only built static ".a" libraries
