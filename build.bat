@echo off
setlocal
cd /d "%~dp0"
if not exist build mkdir build
g++ -std=c++17 -O2 -Wall -Wextra -pedantic -static-libgcc -static-libstdc++ -municode -mwindows -I. main.cpp fluid\FluidEngine.cpp fluid\WorkerPool.cpp rigid\RigidBodyEngine.cpp ui\UiShell.cpp ui\UiLanguage.cpp ui\UiAssets.cpp gas\GasEngine.cpp thermal\ThermalEngine.cpp -o build\water_test.exe -lgdi32 -luser32 -lgdiplus
if errorlevel 1 (
  echo.
  echo Build failed. Close PIPACE if it is still running, then try again.
  echo Or install the MSYS2 UCRT64 toolchain / build with CMake.
  pause
  exit /b 1
)
if not exist build\assets mkdir build\assets
xcopy /y /q assets\*.png build\assets\ >nul
echo Built: %CD%\build\water_test.exe
endlocal
