# flag_resampler
High quality sample rate converter (src) a.k.a resampler. Loads wav and flac and produces wav.

* Overlap-Save Convolution
* Self-Convolved Kernel
* Large Filter-Kernel (depending on source and target rate)
* Support very long files (much larger than your RAM could hold)


### Algorithm / Approach ###

* Integer upsampling is done with zero-padding and convolution with a filter kernel that creates a nice interpolator.  
* Integer downsampling is done with convolution of with a filter kernel that band-limits the signal followed by a trivial decimator (just skip N samples for every sample output).  
* Rational ratios up to a numerator of 31 are handled with a combination of first upsample and bandlimit with a single filter kernel and then a trivial decimation step (the divisor)
* More uneven ratios (such as 147/320, needed for 96k -> 44k1) there is first an upsample and bandlimit-step to get to 32x and then an interpolating sampler is used based on "Polynomial Interpolators for
High-Quality Resampling of
Oversampled Audio" by Olli Niemitalo (2001)


### Quality ###

Quality measurements by "Hydrogen Audio" here:
https://src.hydrogenaudio.org/compareresults?id1=7ad149e6-7495-493b-a9c0-f9b17d8d8a34


### How do I get set up? ###

(For Windows, install Msys2, also the executable becomes "flag_resampler.exe")

```bash
git clone https://github.com/j-norberg/flag_resampler.git
cd flag_resampler
make
```

### Easiest way to use ###

In windows, drag and drop wav or flac files onto the bat-files named "Convert-To-44100.bat" and similar


### Command-line Options ###

use like this:

flag_resampler.exe -i infile.wav -o outfile.wav -r 44100 -f [16, 24, f32, f64]

or:

flag_resampler.exe --input infile.wav --output outfile.wav --sample_rate 44100 --format [16, 24, f32, f64]

-i / --input is followed by a path to wave/flac file

-o / --output is followed by a path to wave file that will write to

-r / --sample_rate followed by desired sample rate (for instance 44100, 48000, 96000, 192000)

-q / --quality followed by a percentage (default is 100)

-f / --format are one of these:

16 = 16 bit dithered (useful for CD)

24 = 24 bit dithered (useful for DVD/Blu-ray)

f32 = 32 bit IEEE float (this is the default)

f64 = 64 bit IEEE float (this is the best quality)


