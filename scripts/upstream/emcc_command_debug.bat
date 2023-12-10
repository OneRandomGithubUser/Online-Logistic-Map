SET emcc_main=emcc fftw-3.3.10/libs/libfftw3.a src/main.cpp -o build/main.js
SET emcc_lib=-s USE_BOOST_HEADERS=1 -std=c++20 -l embind -pthread
SET emcc_options=-s PTHREAD_POOL_SIZE=1 -s ALLOW_MEMORY_GROWTH -s NO_DISABLE_EXCEPTION_CATCHING
SET emcc_optimization=-g
SET emcc_simdDebug=-Rpass=loop-vectorize -Rpass-missed=loop-vectorize -Rpass-analysis=loop-vectorize -gline-tables-only -gcolumn-info
SET emcc_command=%emcc_main% %emcc_lib% %emcc_options% %emcc_optimization%