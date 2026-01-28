#include <chrono>
#include <cstdint>
#include <cstdio>
#include <cstring>
#include <iostream>
#include <string>
#include <vector>

#include "streamer.h"

#include "memory_reader.h"
#include "file_reader.h"
#include "writer.h"

#define _USE_MATH_DEFINES 
#include <math.h>

// some constants
const char* VERSION_STRING = "FLAG Sample Rate Converter (flag_resampler) version 0.88";

const char* USAGE_STRING =
R"(---------------------------
use like this:

flag_resampler.exe -i infile.wav -o outfile.wav -r 44100 -f [16, 24, f32]

or:

flag_resampler.exe --input infile.wav --output outfile.wav --sample_rate 44100 --format [16, 24, f32]

-i / --input is followed by a path to wave file
-o / --output is followed by a path to wave file that will write to
-r / --sample_rate followed by desired sample rate (for instance 44100, 48000, 96000)

-f / --format are one of these:
16 = 16 bit dithered (useful for CD)
24 = 24 bit dithered (useful for DVD/Blu-ray)
f32 = 32 bit IEEE float (this is the default)
f64 = 64 bit IEEE float (this is the best quality)
)";

class Timer
{
public:
	Timer() :
		_beg(clock_::now())
	{
	}
	void reset()
	{
		_beg = clock_::now();
	}

	double elapsed_ms() const
	{
		std::chrono::duration<double, std::milli> ms = clock_::now() - _beg;
		return ms.count();
	}

private:
	typedef std::chrono::high_resolution_clock clock_;
	typedef std::chrono::duration<double, std::ratio<1> > second_;
	std::chrono::time_point<clock_> _beg;
};


















struct Options
{
	// bool test = false;
	bool show_version = false;
	bool show_usage = false;

	std::string in_file;
	std::string out_file;
	Writer::OutFormat _format = Writer::eFmtFloat;

	int out_sr = 44100;
};

const static char* next_arg(int argc, const char** argv, int& index, const char* error_msg)
{
	if (index < (argc - 1))
	{
		++index;
		return argv[index];
	}

	puts(error_msg);
	return nullptr;
}

bool handle_flag(Options& o, int argc, const char** argv, int& index)
{
	std::string arg(argv[index]);

#if 0
	if (arg == "-t" || arg == "--test")
	{
		o.test = true;
		return true;
	}
#endif

	if (arg == "-v" || arg == "--version")
	{
		o.show_version = true;
		return true;
	}

	if (arg == "-i" || arg == "--input")
	{
		o.in_file = next_arg(argc, argv, index, "error: -i should be followed by filename\n");
		return o.in_file.empty() == false;
	}

	if (arg == "-o" || arg == "--output")
	{
		o.out_file = next_arg(argc, argv, index, "error: -o should be followed by filename\n");
		return o.out_file.empty() == false;
	}

	if (arg == "-r" || arg == "--sample_rate")
	{
		const char* sr_str = next_arg(argc, argv, index, "error: -r should be followed by sample rate\n");
		o.out_sr = sr_str ? atoi(sr_str) : 0; // set to 0 to indicate issue
		return o.out_sr > 0;
	}

	if (arg == "-f" || arg == "--format")
	{
		const char* format_error_message = "error: -f should be followed by format (16, 24, f32, f64)\n";
		const char* fmt_str = next_arg(argc, argv, index, format_error_message);
		if (fmt_str == nullptr)
			return false;

		if (0 == strcmp(fmt_str, "16"))
		{
			o._format = Writer::eFmtInt16Dithered;
			return true;
		}

		if (0 == strcmp(fmt_str, "24"))
		{
			o._format = Writer::eFmtInt24Dithered;
			return true;
		}

		if (0 == strcmp(fmt_str, "f32"))
		{
			o._format = Writer::eFmtFloat;
			return true;
		}

		if (0 == strcmp(fmt_str, "f64"))
		{
			o._format = Writer::eFmtDouble;
			return true;
		}

		puts(format_error_message);
		return false;
	}

	printf("unhandled flag %s\n\n", arg.c_str());
	return false;
}


Options get_options(int argc, const char** argv)
{
	// fixme more reporting about invalid flags
	Options o;

	for (int i = 1; i < argc; ++i)
		if (!handle_flag(o, argc, argv, i))
		{
			o.out_file.clear(); // make invalid
			break;
		}

	// when do we need to show usage?
//	if (!o.test)
	{
		if (o.in_file.empty() || o.out_file.empty())
			o.show_usage = true;
	}

	return o;
};





















void put_progress(int pct)
{
	for (int i = 0; i < pct; ++i)
		putchar('>');
	for (int i = pct; i < 100; ++i)
		putchar('.');
	putchar('\r');
	fflush(stdout);
}


// fixme maybe other tools, like generate sweeps?
int main(int argc, const char** argv)
{
	Options s = get_options(argc, argv);

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

	int64_t in_frame_count = reader->get_frame_count();
	int in_channel_count = reader->get_channel_count();
	int in_sample_rate = reader->get_sample_rate();
	double in_seconds = (double)in_frame_count / (double)in_sample_rate;

	printf("INPUT: rate=%d, frames=%lld, seconds=%.3f, channel-count=%d\n", in_sample_rate, in_frame_count, in_seconds, in_channel_count);

	// calculate out-file samples, round to closest
	int64_t out_frame_count = (s.out_sr * in_frame_count + (in_sample_rate/2)) / in_sample_rate;

	printf("OUTPUT: rate=%d, frames=%lld\n", s.out_sr, out_frame_count);

	ISampleProducer* streamer = nullptr;

	// check special case
	if (in_sample_rate == s.out_sr)
	{
		puts("sample-rates same, will still output file");
		streamer = reader;
	}
	else
	{
		//
		streamer = streamer_factory(reader, s.out_sr);
	}

	if (streamer == nullptr)
	{
		// issue?
		puts("can not create converter?");
		return 0;
	}

	int64_t hacked_frames = out_frame_count + 50000; // fixme
	Writer writer(s.out_file.c_str(), hacked_frames, streamer, s._format);

	Timer t1;
	t1.reset();

	while (writer.update())
		put_progress(writer.get_progress_percent());

	put_progress(100);

	float elapsed1 = (float)t1.elapsed_ms();
	if (elapsed1 > 1000.0f)
		printf("\nConversion Done in %f s \n", elapsed1/1000.0f);
	else
		printf("\nConversion Done in %f ms \n", elapsed1);

	delete streamer; // deletes their input

	return 0;
}

// test downsample (back to 44k1)
// -r 44100 -i sweep_1_44100-48000.wav -o sweep_1_44100-48000_downsample.wav

// test upsample (to 48k)
// -r 48000 -i sweep_0.wav -o sweep_0_upsample.wav
