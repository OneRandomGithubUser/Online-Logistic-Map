cd ..
call scripts/upstream/emcc_command_debug.bat
call scripts/upstream/emsdk.bat
echo %cd%
echo %~dp0
%emcc_command_powershell%