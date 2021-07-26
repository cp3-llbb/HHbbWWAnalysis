# How to compile the HME library 
To be run inside `HMEStudy`
```
mkdir build
cd build/
cmake -DCMAKE_BUILD_TYPE=Release -Dbamboo_DIR=$(python -c 'from pkg_resources import resource_filename; print(resource_filename("bamboo", "cmake"))')  ..
make -j 4
```
... and profit !
Note : for debugging purposes you may add `-DCMAKE_CXX_FLAGS='-DHME_DEBUG=1'` but be aware that the log dump is huge and therefore should not be applied when running multiple jobs
