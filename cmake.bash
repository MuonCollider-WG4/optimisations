set -e

python_version=3.10
build_dir=build

if [ -d ${build_dir} ]; then
    rm -r ${build_dir}
fi
mkdir ${build_dir}
cmake CMakeLists.txt \
      -B${build_dir}/ -H./ \
      -DPython_INCLUDE_DIR=/usr/include/python${python_version}/ \
      -DPython_LIBRARY=/usr/lib/x86_64-linux-gnu/libpython${python_version}.so \
      -DGSL_INCLUDE_DIR=${HOME}/Software/install/include/ \
      -DGSL_LIBRARY=${HOME}/Software/install/lib/libgsl.a \
      -DCMAKE_INSTALL_PREFIX=${OPAL_THIRD_PARTY}
#-DPYTHON_INSTALL_DIR=${OPAL_THIRD_PARTY}/lib/python/site-packages/
#-DVERBOSE=1 \

cd build
make -j20

