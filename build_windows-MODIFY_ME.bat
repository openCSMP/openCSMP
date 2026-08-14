@echo off
:: 1. Initialize Intel oneAPI Environment
call "C:\Program Files (x86)\Intel\oneAPI\setvars.bat" intel64

:: 2. Create and enter build directory
if not exist build mkdir build
cd build

:: 3. Run CMake using the Intel Toolset for Visual Studio
:: -T identifies the Intel Compiler to the Visual Studio Generator
cmake .. -G "Visual Studio 18 2026" -A x64 -T "Intel C++ Compiler 2026" -Dsamg=ON

:: 4. Build the project
cmake --build . --config Release
pause