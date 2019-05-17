@echo off

echo --- Start process
echo --- Setting initial variables

set OUT_DIR=out
set MPI_DIR=C:\work\MPI
set PROCESS_NUMBER=4
set COMPILE_FILE_NAME=3d_heat_spread
set OUT_FILE_NAME=input1.txt
set OUT_FILE_DIR=data

where cl.exe >nul 2>nul
if %ERRORLEVEL% NEQ 0 (
    echo --- Source compiler
    call "%VS2017INSTALLDIR%\VC\Auxiliary\Build\vcvars64.bat"
) else (
    echo --- Compiler already sourced
)

echo --- Start building
echo.

rem /openmp
cl %COMPILE_FILE_NAME%.c /EHsc                 ^
   /Fo%OUT_DIR%\ /Fe%OUT_DIR%\                 ^
   /I%MPI_DIR%\Include /I%MPI_DIR%\Include\x64 ^
   /link                                       ^
   /LIBPATH:%MPI_DIR%\Lib\x64 msmpi.lib

echo.
echo --- Run
echo.

%MPI_DIR%\bin\mpiexec.exe -n %PROCESS_NUMBER% %OUT_DIR%\%COMPILE_FILE_NAME%.exe %OUT_FILE_DIR%\%OUT_FILE_NAME%
