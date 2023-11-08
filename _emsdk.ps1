git submodule update --recursive
cd emsdk
$execution_policy = Get-ExecutionPolicy -Scope Process
Set-ExecutionPolicy RemoteSigned -Scope Process
./emsdk install latest
./emsdk activate latest
Set-ExecutionPolicy $execution_policy -Scope Process
cd ..
emcc fftw-3.3.10/libs/libfftw3.a main.cpp -o main.js -s USE_BOOST_HEADERS=1 -std=c++20 -l embind -pthread -s PTHREAD_POOL_SIZE=1 -s ALLOW_MEMORY_GROWTH -O3
Pause