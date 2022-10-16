wt -d %~dp0 powershell -NoExit Add-Content -path (Get-PSReadlineOption).HistorySavePath 'cls\; emcc main.cpp -o main.js -s USE_BOOST_HEADERS=1 -std=c++20 -lembind -g -sNO_DISABLE_EXCEPTION_CATCHING -sALLOW_MEMORY_GROWTH'
:: wt -d %~dp0 opens Windows Terminal to the current directory
:: powershell -NoExit opens powershell without exiting after the command is complete
:: NOTE: emsdk must be activated permanently by running .\emsdk activate latest --permanent
:: Add-Content -path (Get-PSReadlineOption).HistorySavePath adds the command to be run to the console history so that it can be run more easily
:: cls\; to clear the output of the previous compile, with the semicolon escaped with a backslash so that windows terminal does not interpret it as another command
:: 'emcc path\to\cpp\file -o path\of\output\js\file -s USE_BOOST_HEADERS=1 -std=c++20 -lembind -g -sNO_DISABLE_EXCEPTION_CATCHING -fmodule-map-file=''kynedraw.modulemap'' -fmodules' compiles the file with debug commands