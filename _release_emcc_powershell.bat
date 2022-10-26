powershell -NoExit Add-Content -path (Get-PSReadlineOption).HistorySavePath 'cls; emcc main.cpp -o main.js -s USE_BOOST_HEADERS=1 -std=c++20 -lembind -sALLOW_MEMORY_GROWTH -O3 -sPTHREAD_POOL_SIZE=2 -pthread'
:: powershell -NoExit opens powershell without exiting after the command is complete
:: Add-Content -path (Get-PSReadlineOption).HistorySavePath adds the command to be run to the console history so that it can be run more easily
:: cls; to clear the output of the previous compile
:: emcc main.cpp -o main.js -s USE_BOOST_HEADERS=1 -std=c++20 -lembind -sALLOW_MEMORY_GROWTH -O3 -sPTHREAD_POOL_SIZE=2 -pthread compiles the file without debug commands, with optimization
:: NOTE: emsdk must be activated permanently by running .\emsdk activate latest --permanent