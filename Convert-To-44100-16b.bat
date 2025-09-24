
@echo off
title FLAG_SRC 44100-16

REM go to the directory of the bat file
cd /d "%~dp0"

:again
if "%~1" == "" goto done

flag_resampler.exe -i "%~1" -o "%~1".44100.16b.wav -r 44100 -f 16

shift
goto again

:done
pause
exit


