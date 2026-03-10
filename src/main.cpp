#pragma warning( disable : 4514 ) // unref. inline
#pragma warning( disable : 4710 ) // not inlined
#pragma warning( disable : 4711 ) // inline expansion
#pragma warning( disable : 4820 ) // padding
#pragma warning( disable : 5045 ) // spectre

#include <chrono>
#include <cstdint>
#include <cstdio>
#include <cstring>
#include <iostream>
#include <string>
#include <vector>

#define _USE_MATH_DEFINES 
#include <math.h>

#include "streamer.h"

#include "memory_reader.h"

#include "file_reader.h"
#include "writer.h"

#include "options.inl"
#include "timer.inl"


// some constants
const char* VERSION_STRING = "FLAG Sample Rate Converter (flag_resampler) version 0.88";

const char* USAGE_STRING =
R"(---------------------------
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
)";

void put_progress(int pct)
{
	for (int i = 0; i < pct; ++i)
		putchar('>');
	for (int i = pct; i < 100; ++i)
		putchar('.');
	putchar('\r');
	fflush(stdout);
}

bool perform_conversion(const std::string& out_file, int out_sr, int quality, Writer::OutFormat format, ISampleProducer* reader, int64_t in_frame_count, bool is_flac)
{
	// time the whole thing, setup and all
	Timer t1;
	t1.reset();

	int in_channel_count = reader->get_channel_count();
	int in_sample_rate = reader->get_sample_rate();
	double in_seconds = (double)in_frame_count / (double)in_sample_rate;

	printf("INPUT: channel-count=%d, rate=%d, frames=%lld (%.3fs) (flac=%d)\n", in_channel_count, in_sample_rate, in_frame_count, in_seconds, is_flac ? 1 : 0);

	// calculate out-file samples, round to closest
	int64_t out_frame_count = (out_sr * in_frame_count + (in_sample_rate / 2)) / in_sample_rate;

	printf("OUTPUT: rate=%d, frames=%lld\n", out_sr, out_frame_count);

	ISampleProducer* streamer = nullptr;

	// check special case
	if (in_sample_rate == out_sr)
	{
		puts("note: sample-rates same, will still output file");
		streamer = reader;
	}
	else
	{
		//
		streamer = streamer_factory(reader, out_sr, quality);
	}

	if (streamer == nullptr)
	{
		// issue?
		puts("error: can not create converter?");
		return 0;
	}

	Writer writer(out_file.c_str(), out_frame_count, streamer, format);

	while (writer.update())
		put_progress(writer.get_progress_percent());

	float elapsed1 = (float)t1.elapsed_ms();

	put_progress(100);

	if (elapsed1 > 1000.0f)
		printf("\nConversion Done in %f s \n", elapsed1 / 1000.0f);
	else
		printf("\nConversion Done in %f ms \n", elapsed1);

	delete streamer; // deletes their input

	return true;
}

std::vector<double> generate_sine_sweep(
	double start_freq = 20.0,
	double end_freq = 48000.0,
	double duration = 10.0,
	double sample_rate = 96000.0)
{
	const size_t numSamples = static_cast<size_t>(duration * sample_rate);
	std::vector<double> sweep(numSamples);

	// Rate of frequency increase per second
	const double chirp_rate = (end_freq - start_freq) / duration;

	for (size_t i = 0; i < numSamples; ++i)
	{
		const double t = static_cast<double>(i) / sample_rate;
		const double phase = 2.0 * M_PI * (start_freq * t + 0.5 * chirp_rate * t * t);
		sweep[i] = std::cos(phase);
	}

	return sweep;
}

void run_tests()
{
	puts("run tests"); fflush(stdout);

	// generate sweep at 96k, 44k1
	std::vector<double> v0 = generate_sine_sweep(20,48000,60, 96000);

	// check performance doing 96k->44k1
	MemoryReader* reader = new MemoryReader(v0, 96000, 1);
	bool is_flac = false;
	int64_t in_frame_count = reader->get_frame_count();
	bool ok = perform_conversion("test.wav", 44100, 100, Writer::eFmtFloat, reader, in_frame_count, is_flac);

	// compare

//	std::vector<double> v1 = generate_sine_sweep(20, 48000, 60, 44100); // will alias, fade out as post-processing

	puts("done"); fflush(stdout);
}

// fixme maybe other tools, like generate sweeps?
int main(int argc, const char** argv)
{
	Options s = get_options(argc, argv);

	if (s.test)
	{
		run_tests();
		return 0;
	}

	if (s.show_version || s.show_usage)
	{
		puts(VERSION_STRING);
	}

	if (s.show_usage)
	{
		puts(USAGE_STRING);
	}

#if 0
	if (s.test)
	{
		// just run tests
		run_tests();
		return 0;
	}
#endif

	if (s.out_sr < 1)
	{
		// error about bad sample-rate
		puts("bad sample rate");
		return -1;
	}

	if (s.quality < 1)
	{
		// error about bad quality
		puts("bad quality");
		return -1;
	}


	//
	if (s.out_file.empty())
	{
		// error about missing output
		return -1;
	}

	// check input file
	FileReader* reader = new FileReader(s.in_file.c_str());
	if (reader->get_error_code() != FileReader::eOk)
	{
		printf("can not open file %s\n", s.in_file.c_str());
		// can not read file
		return -1;
	}

	perform_conversion(s.out_file, s.out_sr, s.quality, s.format, reader, reader->get_frame_count(), reader->is_flac());

	return 0;
}
