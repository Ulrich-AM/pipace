@echo off
setlocal
cd /d "%~dp0"
call build.bat
if errorlevel 1 exit /b 1
if exist build\water_test.exe start "" /high build\water_test.exe
endlocal
