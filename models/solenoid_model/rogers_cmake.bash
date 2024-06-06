set -e

python_version=3.10
build_dir=build

if [ -d ${build_dir} ]; then
    rm -r ${build_dir}
fi
mkdir ${build_dir}
cmake CMakeLists.txt \
      -B${build_dir}/ -H./ \
      -DVERBOSE=1 \
      -DPython_LIBRARY=/usr/lib/x86_64-linux-gnu/libpython${python_version}.so \
      -DPython_INCLUDE_DIR=/usr/include/python${python_version}/ \
      -DCMAKE_INSTALL_PREFIX=${OPAL_THIRD_PARTY} \
      -DPYTHON_INSTALL_DIR=${OPAL_THIRD_PARTY}/lib/python/site-packages/
cd build
make -j20

