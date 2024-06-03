set -e

python_version=3.10

rm -r build
mkdir build
cmake CMakeLists.txt \
      -Bbuild/ -H./ \
      -DVERBOSE=1 \
      -DPython_LIBRARY=/usr/lib/x86_64-linux-gnu/libpython${python_version}.so \
      -DPython_INCLUDE_DIR=/usr/include/python${python_version}/ \
      -DCMAKE_INSTALL_PREFIX=${OPAL_THIRD_PARTY} \
      -DPYTHON_INSTALL_DIR=${OPAL_THIRD_PARTY}/lib/python/site-packages/
cd build
make -j20

