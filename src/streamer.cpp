#include "dep/pffft.h"

#include <vector>
#include <cassert>
#include <cstring> // for memset
#include <string> // for debugging

#define _USE_MATH_DEFINES
#include <math.h> // for M_PI

#include <limits> // for epsilon

#include "ok_interp.inl"
#include "streamer.h"
#include "sample_producer.h"


// fixme todo jnorberg don't need this in reality
extern bool simple_wav_write_mono(const char* file_name, float* buf, int64_t total_frame_count);

// hamming is bad

// would like to test:
// Blackman–Harris window
// Dolph–Chebyshev window
// Generalized adaptive polynomial (GAP) window

#if 0
// blackman
inline double win_funcion(double f)
{
	double p = 0.5 * f + 0.5; // translate from whole window to half-window
	double v =
		+ 0.42
		- 0.50 * cos(p * M_PI * 2.0) +
		+ 0.08 * cos(p * M_PI * 4.0);
	return v;
}
#endif


#if 1
//_blackman_nuttall
inline double win_funcion(double f)
{
	double p = 0.5 * f + 0.5; // translate from whole window to half-window
	double v =
		+0.3635819
		- 0.4891775 * cos(p * M_PI * 2.0)
		+ 0.1365995 * cos(p * M_PI * 4.0)
		- 0.0106411 * cos(p * M_PI * 6.0);

	return v;
}
#endif






// 
void create_windowed_sinc(float* dest, int count, double window_time_scale, double window_amp_scale)
{
	constexpr double epsilon = std::numeric_limits<double>::epsilon();
	double recip = 1.0 / count;
	double mul = M_PI * count * 0.5 * window_time_scale;
	for (int i = 0; i < count; ++i)
	{
		double p01 = (i + 0.5) * recip; // 0..1
		double x = (p01 * 2.0) - 1.0;   // -1 .. 1
		double win = win_funcion(x) * window_amp_scale;
		double pN = x * mul; // normalized -> sin
		if (abs(pN) > epsilon)
		{
			double sinc = sin(pN) / pN;
			double v = sinc * win;
			dest[i] = (float)v;
		}
		else
		{
			dest[i] = (float)win;
		}
	}
}







template<int K_BUFS>
struct BufsT
{
	// this is the only one the is not fully temp
	float _in_buf_time[K_BUFS];	// input goes here and is sliding in to this buffer from the back

	// to ensure no aliasing, all buffers have to be explicit
	float _buf_filter[K_BUFS]; // should only allocate this one
	float _buf_work[K_BUFS]; // this one points into buf above

	float _in_buf_freq[K_BUFS];	// input is transformed into frequences

	float _out_buf_freq[K_BUFS];	// result of convolve goes here
	float _out_buf_time[K_BUFS];	// IFFT here
};

// 
template<int TRANSFORM_SIZE>
struct OverlapSaveStructT
{
	typedef BufsT<TRANSFORM_SIZE> Bufs;

	int _filter_size = 0;
	int _inblock_size = 0;
	float _transform_recip = 1.0f;

	PFFFT_Setup* _fft; // could be shared btw all instances

	Bufs *_bufs;

	OverlapSaveStructT()
	{
		_fft = nullptr;
		_bufs = nullptr;
	}

	~OverlapSaveStructT()
	{
		pffft_aligned_free(_bufs);
		pffft_destroy_setup(_fft);
	}

	void init(float* temp_windowed_sinc, int filter_size)
	{
		assert(_bufs == nullptr);
		assert(_fft == nullptr);

		_filter_size = filter_size;
		_inblock_size = TRANSFORM_SIZE - filter_size;
		_transform_recip = 1.f / TRANSFORM_SIZE;

		_fft = pffft_new_setup(TRANSFORM_SIZE, PFFFT_REAL);

		_bufs = (Bufs*)pffft_aligned_malloc(sizeof(Bufs));
		memset(_bufs, 0, sizeof(Bufs));

		pffft_transform(_fft, temp_windowed_sinc, _bufs->_buf_filter, _bufs->_buf_work, PFFFT_FORWARD);
	}

	// 1. fill buffer (call this K_INBUF_SIZE times)
	inline void write(int index, float v)
	{
		_bufs->_in_buf_time[_filter_size + index] = v;
	}

	// 2. call this
	void convolve()
	{
		// FFT
		pffft_transform( _fft, _bufs->_in_buf_time, _bufs->_in_buf_freq, _bufs->_buf_work, PFFFT_FORWARD);

		// scoot the input-buffer
		memmove(_bufs->_in_buf_time, _bufs->_in_buf_time + _inblock_size, _filter_size * sizeof(float));

		// multiply with "kernel"
		pffft_zconvolve_no_accumulate(_fft, _bufs->_in_buf_freq, _bufs->_buf_filter, _bufs->_out_buf_freq, _transform_recip);

		// IFFT
		pffft_transform(_fft, _bufs->_out_buf_freq, _bufs->_out_buf_time, _bufs->_buf_work, PFFFT_BACKWARD);
	}

	// 3. read from here (call this K_INBUF_SIZE times)
	inline float read(int index)
	{
		return _bufs->_out_buf_time[_filter_size + index];
	}
};








//
template<int TRANSFORM_SIZE>
struct StreamerUpT : ISampleProducer
{
	typedef OverlapSaveStructT<TRANSFORM_SIZE> OverlapSaveStruct;

	// Odd filter-length makes it easier to align upsampling
	int _filter_size;
	int _inblock_size;
	int _up;

	int _pad_ttl; // will zero-pad the input

	int _out_index;

	ISampleProducer* _input;
	int _channel_count;
	int _sr;

	std::vector<float> _channels_buf;
	OverlapSaveStruct* _oss;
	
	StreamerUpT(float* filter_kernel, int filter_size, int up, ISampleProducer* input)
	{
		_filter_size = filter_size;
		_inblock_size = TRANSFORM_SIZE - filter_size;
		_up = up;

		_pad_ttl = 0;
		_out_index = _inblock_size; // force initial feed to buffer

		_input = input;
		_channel_count = input->get_channel_count();
		_sr = _input->get_sample_rate() * up;

		// shared thing
		_oss = new OverlapSaveStruct[_channel_count];
		for (int i = 0; i < _channel_count; ++i)
			_oss[i].init(filter_kernel, filter_size);

		_channels_buf.resize(_channel_count);

		// pre-feed to counteract internal latency
		const int kInterpolatedSamplerHalf = 3;
		int pre_feed = _filter_size / 2 - kInterpolatedSamplerHalf;
		skip_next(pre_feed);
	}

	~StreamerUpT()
	{
		delete [] _oss;
		delete _input;
	}

	int get_sample_rate() override { return _sr; };
	int get_channel_count() override { return _channel_count; }

	inline void fill_oss_buffer()
	{
		// 1. de-interleave
		// 2. zero-padding

		for (int i = 0; i < _inblock_size; ++i)
		{
			if (_pad_ttl < 1)
			{
				// get one single (interleaved) frame for all channels
				_input->get_next(_channels_buf.data(), 1);
				
				// populate oss-structs
				for (int c = 0; c < _channel_count; ++c)
					_oss[c].write(i, _channels_buf[c]);

				// reset ttl
				_pad_ttl = _up;
			}
			else
			{
				// zero-pad
				for (int c = 0; c < _channel_count; ++c)
					_oss[c].write(i, 0);
			}

			--_pad_ttl;
		}

		for (int c = 0; c < _channel_count; ++c)
			_oss[c].convolve();
	}

	void get_next(float* buf_interleaved, int64_t frame_count) override
	{
		int write_index = 0;
		for (int i = 0; i < frame_count; ++i)
		{
			if (_out_index >= _inblock_size)
			{
				fill_oss_buffer();
				_out_index = 0;
			}

			for (int c = 0; c < _channel_count; ++c)
			{
				buf_interleaved[write_index] = _oss[c].read(_out_index);
				++write_index;
			}

			++_out_index;
		}
	}

	void skip_next(int64_t frame_count) override
	{
		for (int i = 0; i < frame_count; ++i)
		{
			if (_out_index >= _inblock_size)
			{
				fill_oss_buffer();
				_out_index = 0;
			}

			++_out_index;
		}
	}

};











// very simple
struct TrivialDecimator : ISampleProducer
{
	ISampleProducer* _input;
	int _channel_count;
	int _out_sr;

	int _read_ttl;
	int _skipFrames;

	TrivialDecimator(int multiplier, ISampleProducer* input)
	{
		_input = input;
		_channel_count = input->get_channel_count();
		_out_sr = input->get_sample_rate() / multiplier;

		_read_ttl = 0;
		_skipFrames = multiplier - 1;

		// match latency to interpolated sampler
		_input->skip_next(3);
	}

	~TrivialDecimator()
	{
		delete _input;
	}

	int get_sample_rate() override { return _out_sr; };
	int get_channel_count() override { return _channel_count; }

	void get_next(float* buf_interleaved, int64_t frame_count) override
	{
		for (int i = 0; i < frame_count; ++i)
		{
			// read one frame
			_input->get_next(buf_interleaved, 1);
			buf_interleaved += _channel_count;

			_input->skip_next(_skipFrames);
		}
	}

	void skip_next(int64_t frame_count) override
	{
		int64_t skip = frame_count * (_skipFrames + 1);
		_input->skip_next(skip);
	}

};



struct InterpolatedSampler : ISampleProducer
{
	ISampleProducer* _input;
	int _channel_count;
	int _out_sr;

	enum {
		k_bufs = 64,
		k_mask = k_bufs - 1,
		k_half = k_bufs >> 1
	};

	struct Buf {
		float _b[k_bufs];
	};

	std::vector<float> _channel_buffer;

	Buf* _bufs;

	int _read_index;
	int _read_index_frac;
	int _read_index_frac_add;
	int _read_index_frac_limit;
	float _read_index_frac_limit_recip;

	InterpolatedSampler(int sr, ISampleProducer* input)
	{
		_input = input;
		_channel_count = input->get_channel_count();
		_out_sr = sr;

		_read_index = 0;

		_read_index_frac = 0;
		_read_index_frac_add = input->get_sample_rate();
		_read_index_frac_limit = sr;
		_read_index_frac_limit_recip = 1.0f / sr;

		_channel_buffer.resize(_channel_count);
		_bufs = new Buf[_channel_count];

		// fill whole buffer
		for (int i = 0 ; i < k_bufs ; ++i)
		{
			input->get_next(_channel_buffer.data(), 1);
			for (int c = 0; c < _channel_count; ++c)
				_bufs[c]._b[i] = _channel_buffer[c];
		}
	}

	~InterpolatedSampler()
	{
		delete[] _bufs;
		delete _input;
	}

	int get_sample_rate() override { return _out_sr; };
	int get_channel_count() override { return _channel_count; }

	void internal_fill_buffers(int ofs)
	{
		if (_channel_count < 2)
		{
			// no need for de-interleave
			_input->get_next(_bufs[0]._b + ofs, k_half);
			return;
		}

		// de-interleave
		for (int i = 0; i < k_half; ++i)
		{
			_input->get_next(_channel_buffer.data(), 1);
			for (int c = 0; c < _channel_count; ++c)
				_bufs[c]._b[ofs+i] = _channel_buffer[c];
		}
	}

	void get_next(float* buf_interleaved, int64_t frame_count) override
	{
		int write_index = 0;
		for (int i = 0; i < frame_count; ++i)
		{
			// calculate frac (common for all channels)
			float frac = (float)_read_index_frac * _read_index_frac_limit_recip;

			// write interleaved
			for (int c = 0; c < _channel_count; ++c)
			{
				buf_interleaved[write_index] = sample_16x_6p_5z(_bufs[c]._b, k_mask, _read_index, frac);
				++write_index;
			}

			// move index
			_read_index_frac += _read_index_frac_add;

			// handle wrapping (maybe trust modulo here)
			while (_read_index_frac >= _read_index_frac_limit)
			{
				++_read_index;
				_read_index_frac -= _read_index_frac_limit;

				switch (_read_index)
				{
				case k_half:
					internal_fill_buffers(0);
					break;
				case k_bufs:
					internal_fill_buffers(k_half);
					_read_index = 0;
					break;
				}
			}

		}
	}

	void skip_next(int64_t frame_count) override
	{
		for (int i = 0; i < frame_count; ++i)
		{
			// move index
			_read_index_frac += _read_index_frac_add;

			// handle wrapping (maybe trust modulo here)
			while (_read_index_frac >= _read_index_frac_limit)
			{
				++_read_index;
				_read_index_frac -= _read_index_frac_limit;

				switch (_read_index)
				{
				case k_half:
					internal_fill_buffers(0);
					break;
				case k_bufs:
					internal_fill_buffers(k_half);
					_read_index = 0;
					break;
				}
			}

		}
	}

};


void create_filter(float* out_kernel_buffer, int half_len, int kernel_len, float scale_len, float scale_amp)
{
	// create initial filter
	int len1 = 1 + half_len * 2;
	int len2 = 1 + len1 * 2;

	// pot (power of two)
	int transform_size = 1024;
	while (transform_size < len2)
		transform_size <<= 1;

	float transform_recip = 1.0f / transform_size;

	assert(len2 == kernel_len);

	float* filter = (float*)pffft_aligned_malloc(transform_size * sizeof(float));

	create_windowed_sinc(filter, len1, scale_len, scale_amp);
	int zero_count = (transform_size - len1);
	memset(filter + len1, 0, zero_count * sizeof(float));

	// debug
//	simple_wav_write_mono("flt_1.wav", filter, transform_size);

	// self convolve
	float* buf_freq = (float*)pffft_aligned_malloc(transform_size * sizeof(float));
	float* buf_work= (float*)pffft_aligned_malloc(transform_size * sizeof(float));

	PFFFT_Setup* fft = pffft_new_setup(transform_size, PFFFT_REAL);

	pffft_transform(fft, filter, buf_freq, buf_work, PFFFT_FORWARD);
	pffft_zconvolve_no_accumulate(fft, buf_freq, buf_freq, buf_freq, transform_recip); // self-convolve
	pffft_transform(fft, buf_freq, out_kernel_buffer, buf_work, PFFFT_BACKWARD);

	pffft_destroy_setup(fft);
	pffft_aligned_free(buf_work);
	pffft_aligned_free(buf_freq);
	pffft_aligned_free(filter);
}









enum
{
	k_bits = 18, // 88 ms  
	k_transform_len = 1 << k_bits
};

ISampleProducer* make_integer_upsampler(int up, float bw, ISampleProducer* input)
{
	int filter_1_half_len = 640 + up * 460; // does this make sense?
	int filter_1_len = filter_1_half_len * 2 + 1;
	int filter_2_len = filter_1_len * 2 + 1;

	// fixme this should go away
	printf("\nOSS transform=%d, filter=%d, inbuf=%d\n", k_transform_len, filter_2_len, (k_transform_len - filter_2_len));

	assert(k_transform_len > (filter_2_len + 1024));

	float* filter_kernel = (float*)pffft_aligned_malloc(k_transform_len * sizeof(float)); // filter_kernel is only used in init
	memset(filter_kernel, 0, k_transform_len * sizeof(float));

	float up_sqrt = sqrtf((float)up);
	create_filter(filter_kernel, filter_1_half_len, filter_2_len, bw / up, bw / up_sqrt );

	// debug
//	simple_wav_write_mono("flt_2.wav", filter_kernel, k_transform_len);


	ISampleProducer* upsampler = new StreamerUpT<k_transform_len>(filter_kernel, filter_2_len, up, input);
	pffft_aligned_free(filter_kernel); // filter_kernel is only used in init

	return upsampler;
}








// 11 -> 2048
// 12 -> 4096
// 13 -> 8192
// 14 -> 16k
// 15 -> 32k
// 16 -> 64k
// 17 -> 128k
// 18 -> 256k

ISampleProducer* streamer_factory(ISampleProducer* input, int sr_out)
{
	int sr_in = input->get_sample_rate();

	// fixme bw can be input
	float bw = 0.99f;
	if (sr_out < sr_in)
		bw *= (float)sr_out / (float)sr_in;

	// try and find integer ratio (too high becomes very expensive)
	for (int up = 1; up < 16; ++up)
	{
		for (int decimate = 1; decimate < 256; ++decimate)
		{
			if (sr_out * decimate == sr_in * up)
			{
				printf("Rational: %d / %d\n", up, decimate);
				// found rational, very cool, avoids interpolation
				return new TrivialDecimator(decimate, make_integer_upsampler(up, bw, input));
			}
		}

	}

	return new InterpolatedSampler(sr_out, make_integer_upsampler(16, bw, input) );
}


