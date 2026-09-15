@echo off
setlocal
cd /d "%~dp0"
if not exist build mkdir build
g++ -std=c++17 -O2 -Wall -Wextra -pedantic -static-libgcc -static-libstdc++ -municode -mwindows -I. main.cpp fluid\FluidEngine.cpp fluid\WorkerPool.cpp rigid\RigidBodyEngine.cpp ui\UiShell.cpp ui\UiLanguage.cpp ui\UiAssets.cpp ui\ToolWindow.cpp gas\GasEngine.cpp thermal\ThermalEngine.cpp render\WorldRenderer.cpp substance\SubstanceRegistry.cpp world\WorldQuery.cpp -o build\pipace.exe -lgdi32 -luser32 -lgdiplus
if errorlevel 1 (
  echo.
  echo Build failed. Close PIPACE if it is still running, then try again.
  echo Or install the MSYS2 UCRT64 toolchain / build with CMake.
  pause
  exit /b 1
)
if not exist build\assets mkdir build\assets
xcopy /y /q assets\*.png build\assets\ >nul
echo Built: %CD%\build\pipace.exe
endlocal
