# flag_resampler
High quality sample rate converter (src) a.k.a resampler. Loads wav and flac and produces wav.

* Overlap-Save Convolution
* Self-Convolved Kernel
* Large Filter-Kernel, up to 32000 taps (depending on source and target rate)
* Support very long files (much larger than your RAM could hold)
