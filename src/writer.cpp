#include "sample_producer.h"

#include "writer.h"

#include <cstdlib>
#include <cassert>

// Write to a large-ish buffer before calling fwrite
// also this is how much dither noise is created and looped
enum { k_writer_buf_frames = 4 * 1024 * 1024 };





// Simple but good sounding white noise
// pre-calculate a bunch of this...
struct WhiteNoise
{
	static const int WN_RM = 0x7fff;
	long int next = 1;

	int next_int()
	{
		next = next * 1103515245 + 12345;
		unsigned int v = next >> 16;
		v &= WN_RM;
		return (int)v;
	}

	inline int next_triangle()
	{
		int a = next_int();
		int b = next_int();
		int sum = b - a;
		return sum;
	}
};



inline int ftoi(float v)
{
	v += (v >= 0) ? 0.5f : -0.5f;
	return (int)v;
}

// fixme sse
inline short clamp_short(float v)
{
	const float lim_f = 32767.0f;
	const int lim_i = 32767;

	int vi = ftoi(v * lim_f);
	if (vi < lim_i)
	{
		if (vi > -lim_i)
		{
			return (short)vi;
		}
		return -lim_i;
	}

	return lim_i;
}

inline int32_t clamp_24(float v)
{
	const float lim_f = 8388607.0f;
	const int lim_i = 8388607;

	int32_t vi = ftoi(v * lim_f);
	if (vi < lim_i)
	{
		if (vi > -lim_i)
		{
			return vi;
		}
		return -lim_i;
	}

	return lim_i;
}

inline void write_24(uint8_t* dest, int32_t q)
{
	uint32_t qi = (uint32_t)q;
	dest[0] = (qi >>  0) & 255;
	dest[1] = (qi >>  8) & 255;
	dest[2] = (qi >> 16) & 255;
}

int FormatToBitsPerSample(Writer::OutFormat fmt)
{
	switch (fmt)
	{
	case Writer::eFmtInt16Dithered:
		return 16;

	case Writer::eFmtInt24Dithered:
		return 24;

	case Writer::eFmtFloat:
		return 32;

	default:
		break;
	}

	// error
	assert(false);
	return -1;
}

int FormatToWavFormat(Writer::OutFormat fmt)
{
	switch (fmt)
	{
	case Writer::eFmtInt16Dithered:
	case Writer::eFmtInt24Dithered:
		return DR_WAVE_FORMAT_PCM;

	case Writer::eFmtFloat:
		return DR_WAVE_FORMAT_IEEE_FLOAT;

	default:
		break;
	}

	// error
	assert(false);
	return -1;
}

float* allocate_dither_noise(Writer::OutFormat fmt, int64_t frame_count)
{
	float* dest = nullptr;
	WhiteNoise wn;

	float recip = 1;

	switch (fmt)
	{
	case Writer::eFmtFloat:
		return nullptr;

	case Writer::eFmtInt24Dithered:
		recip = 0.8f / (WhiteNoise::WN_RM * 8388607.0f);
	break;

	case Writer::eFmtInt16Dithered:
		recip = 0.8f / (WhiteNoise::WN_RM * 32767.0f);
		break;

	default:
		break;
	}

	// fixme, little optional notch around 4k? and a shelf at 10k?

	dest = new float[frame_count];
	for (int64_t i = 0; i < frame_count; ++i)
	{
		dest[i] = (float)wn.next_triangle() * recip;
	}

	return dest;
}

Writer::Writer(const char* file_name, int64_t total_frame_count, ISampleProducer* input, OutFormat fmt)
{
	_fmt = fmt;

	_frame_index = 0;
	_total_frame_count = total_frame_count;

	_sample_producer = input;

	// open the file here
	drwav_data_format format;
	format.container = drwav_container_riff;

	format.bitsPerSample = FormatToBitsPerSample(fmt);
	format.format = FormatToWavFormat(fmt);

	int channel_count = input->get_channel_count();

	format.channels = channel_count;
	format.sampleRate = input->get_sample_rate();

	_wav = drwav_open_file_write(file_name, &format);
	// fixme deal with write error

	_buf_interleaved = new float[k_writer_buf_frames * channel_count];

	_buf_interleaved_quantized = nullptr;

	// allocate quantizer-buffer
	switch (fmt)
	{
	case Writer::eFmtFloat:
		break;

	case Writer::eFmtInt16Dithered:
		_buf_interleaved_quantized = new uint8_t[k_writer_buf_frames * channel_count * 2];
		break;

	case Writer::eFmtInt24Dithered:
		_buf_interleaved_quantized = new uint8_t[k_writer_buf_frames * channel_count * 3];
		break;

	default:
		break;
	}

	_buf_dither_noise = allocate_dither_noise(fmt, k_writer_buf_frames);
}

Writer::~Writer()
{
	// close
	drwav_close(_wav);

	delete[] _buf_interleaved;

	if (_buf_interleaved_quantized != nullptr)
		delete[] _buf_interleaved_quantized;

	if (_buf_dither_noise != nullptr)
		delete[] _buf_dither_noise;
}

bool Writer::update()
{
	// are we done yet?
	if (_frame_index >= _total_frame_count)
		return false;

	bool keep_going = true;
	size_t prev_index = _frame_index;
	_frame_index += k_writer_buf_frames;
	if (_frame_index >= _total_frame_count)
	{
		keep_going = false;
		_frame_index = _total_frame_count;
	}

	size_t frame_count = _frame_index - prev_index;

	// read/resample
	_sample_producer->get_next(_buf_interleaved, frame_count);

	//
	int channel_count = _sample_producer->get_channel_count();

	uint8_t* buf_to_write = nullptr;

	// straight write, alt. quantize and write
	switch (_fmt)
	{
	case Writer::eFmtFloat:
		{
			buf_to_write = (uint8_t*)_buf_interleaved;
		}
		break;

	case Writer::eFmtInt24Dithered:
		{
			// fixme SSE
			size_t read_index = 0;
			uint8_t* dest_char = _buf_interleaved_quantized;
			for (size_t fi = 0; fi < frame_count; ++fi)
			{
				float dn = _buf_dither_noise[fi];
				for (int c = 0; c < channel_count; ++c)
				{
					int32_t q = clamp_24(_buf_interleaved[read_index] + dn);
					write_24(dest_char, q);
					dest_char += 3;
					++read_index;
				}
			}

			buf_to_write = _buf_interleaved_quantized;
		}
		break;

	case Writer::eFmtInt16Dithered:
		{
			// fixme SSE
			short* dest = (short*)_buf_interleaved_quantized;
			size_t index = 0;
			for (size_t fi = 0; fi < frame_count; ++fi)
			{
				float dn = _buf_dither_noise[fi];
				for (int c = 0; c < channel_count; ++c)
				{
					dest[index] = clamp_short(_buf_interleaved[index] + dn);
					++index;
				}
			}

			buf_to_write = _buf_interleaved_quantized;
		}
		break;

	default:
		break;
	}

	// write buffer
	drwav_uint64 dummy_for_now = drwav_write(_wav, frame_count * channel_count, buf_to_write);
	(void)dummy_for_now; // fixme check for issues

	// return
	return keep_going;
}

void Writer::update_all()
{
	for (;;)
	{
		if (!update())
			break;
	}
}


bool simple_wav_write_mono(const char* file_name, float* buf, int64_t total_frame_count)
{
	drwav_data_format format;
	format.container = drwav_container_riff;
	format.format = DR_WAVE_FORMAT_IEEE_FLOAT;
	format.channels = 1;
	format.sampleRate = 44100;
	format.bitsPerSample = 32;

	drwav* wav = drwav_open_file_write(file_name, &format);
	if (wav == nullptr)
		return false;

	drwav_uint64 dummy_for_now = drwav_write(wav, total_frame_count, buf);
	drwav_close(wav);

	return true;
}
