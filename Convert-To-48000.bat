
@echo off
title FLAG_SRC 48000

REM go to the directory of the bat file
cd /d "%~dp0"

:again
if "%~1" == "" goto done

flag_resampler.exe -i "%~1" -o "%~1".48000.wav -r 48000

shift
goto again

:done
pause
exit


