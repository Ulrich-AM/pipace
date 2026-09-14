@echo off
setlocal
cd /d "%~dp0"
call build.bat
if errorlevel 1 exit /b 1
if exist build\pipace.exe start "" /high build\pipace.exe
endlocal
