:: NOTE: requires %emcc_command% to be defined
:: NOTE: must be called through %emcc_powershell% or %emcc_windows_terminal% in a batch file

SET emcc_color_text_path=color_output_powershell.ps1

:: Take the content of %emcc_color_text_path% and turn it into one multiline string (Out-String) and then run it (Invoke-Expression)
:: ^^^| resolves to ^| when %emcc_color_text% is set, then resolves to | once %emcc_color_text% is called
SET emcc_color_text=Get-Content -Path %emcc_color_text_path% ^^^| Out-String ^^^| Invoke-Expression

SET emcc_update=^
Write-ColorOutput \"Updating submodules...\" Black Green;^
git submodule update --recursive

SET emcc_activate=^
Write-ColorOutput \"Activating EMCC environment...\" Black Green;^
cd emsdk;^
$execution_policy = Get-ExecutionPolicy -Scope Process;^
Set-ExecutionPolicy RemoteSigned -Scope Process;^
./emsdk install latest;^
./emsdk activate latest;^
Set-ExecutionPolicy $execution_policy -Scope Process;^
cd ..

SET emcc_initialize_message=Write-ColorOutput \"Press up to reenter the compilation command\" Black Green

SET emcc_initialize=%emcc_color_text%;%emcc_update%;%emcc_activate%

SET emcc_clear=cls

:: this requires %emcc_command% to be set beforehand
SET emcc_compile=Write-ColorOutput \"Compiling started...\" Black Green;%emcc_command%
SET emcc_recompile=Write-ColorOutput \"Compiling restarted...\" Black Green;%emcc_command%

SET emcc_add_recompile_command=Add-Content -path (Get-PSReadlineOption).HistorySavePath '%emcc_clear%;%emcc_recompile%'

SET emcc_powershell=powershell -NoExit -Command "%emcc_initialize%;%emcc_compile%;%emcc_initialize_message%;%emcc_add_recompile_command%"
SET emcc_windows_terminal=wt -d %~dp0 %emcc_powershell:;=\;%
SET emcc_windows_terminal=%emcc_windows_terminal:\"=\\\"%
:: %emcc_powershell:;=\;% replaces ; with \; in the command line so that Windows Terminal does not interpret them as new tabs.
:: %emcc_windows_terminal:\"=\\\"% then replaces \" with \\\" so that Windows Terminal interprets it as \"
:: See https://ss64.com/nt/syntax-replace.html for how this is used