@echo off

cd Source
mkdir CompileFolder
copy *.f90 CompileFolder 

@call "C:\Program Files (x86)\Intel\oneAPI\setvars.bat" intel64
@call "C:\Program Files (x86)\Intel\oneAPI\compiler\latest\env\vars.bat" intel64

    
echo Note: For REDIST, add these to 'PATH' in Environment Variables: 
echo			%\Intel\oneAPI\compiler\latest\windows\redist\intel64_win\compiler
echo			%\Intel\oneAPI\mkl\latest\redist\intel64
echo       OR: Run the .exe from the Intel oneAPI command prompt for Intel 64 for VS20xx in the W10 Start Menu
     
cd CompileFolder

ifort mod_CH.f90 mod_jit_dgemm.f90 *.f90  /extend-source:132 /fpp /real-size:64 /D_AVX512 /Qipo /libs:dll /threads /Qmkl:parallel /fast /c
ifort *.obj /exe:_CH

copy *.exe ..\..



del *.obj
del *.mod
del *.exe
del *.f90

exit

ifort mod_CH.f90 mod_jit_dgemm.f90 *.f90  /extend-source:132 /fpp /real-size:64 /Qipo /fast /libs:dll /threads /Qmkl:parallel /QaxCORE-AVX2 /QxCORE-AVX2 /c

ifort mod_CH.f90 mod_jit_dgemm.f90 *.f90  /extend-source:132 /fpp /real-size:64 /D_AVX512 /Qipo /fast /libs:dll /threads /Qmkl:parallel /QxCOMMON-AVX512 /QaxCOMMON-AVX512 /c

echo try also the flags:/QxCORE-AVX512 /QaxCORE-AVX512 /Qaxrocketlake /Qxrocketlake  /QxCOMMON-AVX512 /QaxCOMMON-AVX512

