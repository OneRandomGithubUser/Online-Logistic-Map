powershell -NoExit Add-Content -path (Get-PSReadlineOption).HistorySavePath 'cls; emcc fftw-3.3.10/libs/libfftw3.a main.cpp -o main.js -s USE_BOOST_HEADERS=1 -std=c++20 -lembind -g -sNO_DISABLE_EXCEPTION_CATCHING -sALLOW_MEMORY_GROWTH -sPTHREAD_POOL_SIZE=2 -pthread'
:: powershell -NoExit opens powershell without exiting after the command is complete
:: Add-Content -path (Get-PSReadlineOption).HistorySavePath adds the command to be run to the console history so that it can be run more easily
:: cls; to clear the output of the previous compile
:: emcc fftw-3.3.10/libs/libfftw3.a main.cpp -o main.js -s USE_BOOST_HEADERS=1 -std=c++20 -lembind -g -sNO_DISABLE_EXCEPTION_CATCHING -sALLOW_MEMORY_GROWTH -sPTHREAD_POOL_SIZE=2 -pthread compiles the file with debug commands, without optimizations
:: NOTE: emsdk must be activated permanently by running .\emsdk activate latest --permanent