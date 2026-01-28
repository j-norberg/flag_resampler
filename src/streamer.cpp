#include "dep/pffft-double.h"

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


// hamming is bad
// Blackman–Harris window
// Dolph–Chebyshev window
// Generalized adaptive polynomial (GAP) window

#if 0
// blackman
inline double win_funcion(double p)
{
	double v =
		+ 0.42
		- 0.50 * cos(p * M_PI * 2.0) +
		+ 0.08 * cos(p * M_PI * 4.0);
	return v;
}
#endif


#if 1
//_blackman_nuttall
inline double win_funcion_bn(double p)
{
	double v =
		+ 0.3635819
		- 0.4891775 * cos(p * M_PI * 2.0)
		+ 0.1365995 * cos(p * M_PI * 4.0)
		- 0.0106411 * cos(p * M_PI * 6.0);

	return v;
}
#endif





// 
void create_windowed_sinc(double* dest, int count, double window_time_scale, double window_amp_scale)
{
	constexpr double epsilon = std::numeric_limits<double>::epsilon();
	double recip = 1.0 / count;
	double mul = M_PI * count * 0.5 * window_time_scale;
	for (int i = 0; i < count; ++i)
	{
		double p01 = (i + 0.5) * recip; // 0..1
		double x = (p01 * 2.0) - 1.0;   // -1 .. 1
		double win = win_funcion_bn(p01) * window_amp_scale;
		double pN = x * mul; // normalized -> sin
		if (abs(pN) > epsilon)
		{
			double sinc = sin(pN) / pN;
			double v = sinc * win;
			dest[i] = v;
		}
		else
		{
			dest[i] = win;
		}
	}
}





template<int K_BUFS>
struct BufsT
{
	// this is the only one the is not fully temp
	double _in_buf_time[K_BUFS];	// input goes here and is sliding in to this buffer from the back

	// to ensure no aliasing, all buffers have to be explicit
	double _buf_filter[K_BUFS];
	double _buf_work[K_BUFS];

	double _in_buf_freq[K_BUFS];	// input is transformed into frequences

	double _out_buf_freq[K_BUFS];	// result of convolve goes here
	double _out_buf_time[K_BUFS];	// IFFT here
};

// 
template<int TRANSFORM_SIZE>
struct OverlapSaveStructT
{
	typedef BufsT<TRANSFORM_SIZE> Bufs;

	int _filter_size = 0;
	int _inblock_size = 0;
	double _transform_recip = 1.0f;

	PFFFTD_Setup* _fft; // could be shared btw all instances

	Bufs *_bufs;

	OverlapSaveStructT()
	{
		_fft = nullptr;
		_bufs = nullptr;
	}

	~OverlapSaveStructT()
	{
		pffftd_aligned_free(_bufs);
		pffftd_destroy_setup(_fft);
	}

	void init(double* temp_windowed_sinc, int filter_size)
	{
		assert(_bufs == nullptr);
		assert(_fft == nullptr);

		_filter_size = filter_size;
		_inblock_size = TRANSFORM_SIZE - filter_size;
		_transform_recip = 1.f / TRANSFORM_SIZE;

		_fft = pffftd_new_setup(TRANSFORM_SIZE, PFFFTD_REAL);

		_bufs = (Bufs*)pffftd_aligned_malloc(sizeof(Bufs));
		memset(_bufs, 0, sizeof(Bufs));

		pffftd_transform(_fft, temp_windowed_sinc, _bufs->_buf_filter, _bufs->_buf_work, PFFFTD_FORWARD);
	}

	// 1. fill buffer (call this K_INBUF_SIZE times)
	inline void write(int index, double v)
	{
		_bufs->_in_buf_time[_filter_size + index] = v;
	}

	// 2. call this
	void convolve()
	{
		// FFT
		pffftd_transform( _fft, _bufs->_in_buf_time, _bufs->_in_buf_freq, _bufs->_buf_work, PFFFTD_FORWARD);

		// scoot the input-buffer
		memmove(_bufs->_in_buf_time, _bufs->_in_buf_time + _inblock_size, _filter_size * sizeof(double));

		// multiply with "kernel"
		pffftd_zconvolve_no_accumulate(_fft, _bufs->_in_buf_freq, _bufs->_buf_filter, _bufs->_out_buf_freq, _transform_recip);

		// IFFT
		pffftd_transform(_fft, _bufs->_out_buf_freq, _bufs->_out_buf_time, _bufs->_buf_work, PFFFTD_BACKWARD);
	}

	// 3. read from here (call this K_INBUF_SIZE times)
	inline double read(int index)
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

	std::vector<double> _channels_buf;
	OverlapSaveStruct* _oss;
	
	StreamerUpT(double* filter_kernel, int filter_size, int up, ISampleProducer* input)
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

	virtual int get_padding_frame_count() override
	{
		int external_upsampled_padding = _input->get_padding_frame_count() * _up;
		int internal_padding = (_filter_size - 1) / 2;
		internal_padding -= 1; // where does this come from?
		return external_upsampled_padding + internal_padding;
	}

	void get_next(double* buf_interleaved, int64_t frame_count) override
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

	int _skipFrames;

	TrivialDecimator(int multiplier, ISampleProducer* input)
	{
		_input = input;
		_channel_count = input->get_channel_count();
		_out_sr = input->get_sample_rate() / multiplier;

		_skipFrames = multiplier - 1;

		// throw away all padding
		input->skip_next(input->get_padding_frame_count());
	}

	~TrivialDecimator()
	{
		delete _input;
	}

	int get_sample_rate() override { return _out_sr; };
	int get_channel_count() override { return _channel_count; }

	virtual int get_padding_frame_count() override
	{
		// padding skipped in contructor
		return 0;
	}

	void get_next(double* buf_interleaved, int64_t frame_count) override
	{
		if (_skipFrames == 0)
		{
			// pass-through (common for integer-upscaling)
			_input->get_next(buf_interleaved, frame_count);
			return;
		}

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
		int64_t skip_frames = frame_count * (_skipFrames + 1);
		_input->skip_next(skip_frames);
	}

};



struct InterpolatedSampler : ISampleProducer
{
	ISampleProducer* _input;
	int _channel_count;
	int _out_sr;

	enum {
		k_bufs = 256, // interpolation only uses 6 samples 
		k_mask = k_bufs - 1,
		k_half = k_bufs >> 1
	};

	struct Buf {
		double _b[k_bufs] = { -0.99 }; // fixme should never be visited
	};

	std::vector<double> _channel_buffer;

	Buf* _bufs;

	int _read_index;
	int _read_index_frac;
	int _read_index_frac_add;
	int _read_index_frac_limit;
	double _read_index_frac_limit_recip;

	InterpolatedSampler(int sr, ISampleProducer* input)
	{
		_input = input;
		_channel_count = input->get_channel_count();
		_out_sr = sr;

		_read_index = 0;

		_read_index_frac = 0;
		_read_index_frac_add = input->get_sample_rate();
		_read_index_frac_limit = sr;
		_read_index_frac_limit_recip = 1.0f / (double)sr;

		_channel_buffer.resize(_channel_count);

		_bufs = new Buf[_channel_count];

		// note: induced 2 samples of latency
		// but skip all up to that point
		int padding = input->get_padding_frame_count();
		padding -= 2;
		input->skip_next(padding);

		// fill first half of buffer with actual data
		internal_fill_buffers(0);
		internal_fill_buffers(k_half);
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

	int get_padding_frame_count() override
	{
		// the padding is eaten in constructor and return 0 here
		return 0;
	}

	void get_next(double* buf_interleaved, int64_t frame_count) override
	{
		int write_index = 0;
		for (int frame_index = 0; frame_index < frame_count; ++frame_index)
		{
			// calculate frac (common for all channels)
			double frac = (double)_read_index_frac * _read_index_frac_limit_recip;

			// write interleaved
			for (int c = 0; c < _channel_count; ++c)
			{
				buf_interleaved[write_index] = sample_32x_6p_5z(_bufs[c]._b, k_mask, _read_index, frac);
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
		for (int frame_index = 0; frame_index < frame_count; ++frame_index)
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


void create_filter(double* out_kernel_buffer, int len1, int kernel_len, int transform_size, double bw, double up)
{
	// create initial filter
	double transform_recip = 1.0 / transform_size;

	double* filter = (double*)pffftd_aligned_malloc(transform_size * sizeof(double));

	double scale_len = bw / up;
	double scale_amp = bw / sqrt(up);
	create_windowed_sinc(filter, len1, scale_len, scale_amp);
	int zero_count = (transform_size - len1);

	// clear the rest of the transform
	memset(filter + len1, 0, zero_count * sizeof(double));

	// self convolve
	double* buf_work= (double*)pffftd_aligned_malloc(transform_size * sizeof(double));
	double* buf_freq = (double*)pffftd_aligned_malloc(transform_size * sizeof(double));
	double* buf_freq2 = (double*)pffftd_aligned_malloc(transform_size * sizeof(double));

	PFFFTD_Setup* fft = pffftd_new_setup(transform_size, PFFFTD_REAL);
	pffftd_transform(fft, filter, buf_freq2, buf_work, PFFFTD_FORWARD);

	// create a second sligtly lower (maybe switch back to pure self convolve)
	scale_len = (bw * 0.999) / up;
	scale_amp = (bw * 0.999) / sqrt(up);
	create_windowed_sinc(filter, len1, scale_len, scale_amp);
	pffftd_transform(fft, filter, buf_freq, buf_work, PFFFTD_FORWARD);
	pffftd_zconvolve_no_accumulate(fft, buf_freq2, buf_freq, buf_freq, transform_recip); // convolve
	pffftd_transform(fft, buf_freq, out_kernel_buffer, buf_work, PFFFTD_BACKWARD);

	pffftd_destroy_setup(fft);
	pffftd_aligned_free(buf_freq2);
	pffftd_aligned_free(buf_freq);
	pffftd_aligned_free(buf_work);
	pffftd_aligned_free(filter);
}

// maybe normalize each phase?
// or just the whole kernel
void normalize_filter(double* filter_kernel, int transform_len, int up)
{
	for (int phase = 0; phase < up; ++phase)
	{
		double sum = 0.0;
		for (int i = phase; i < transform_len; i+=up)
			sum += filter_kernel[i];

		double mul = 1.0 / sum;
		
//		uint64_t* mul_p = (uint64_t*)&mul;
//		printf("phase=%i/%i, mul=%g %llx\n", phase, up, mul, *mul_p);

		for (int i = phase; i < transform_len; i += up)
			filter_kernel[i] *= mul;
	}

	// finally full normalize
	{
		double sum = 0.0;
		for (int i = 0; i < transform_len; ++i)
			sum += filter_kernel[i];

		double mul = (double)up / sum;

//		uint64_t* mul_p = (uint64_t*)&mul;
//		printf("FINAL, mul=%g %llx\n", mul, *mul_p);

		for (int i = 0; i < transform_len; ++i)
			filter_kernel[i] *= mul;
	}

}







// 11 -> 2048
// 12 -> 4096
// 13 -> 8192
// 14 -> 16k
// 15 -> 32k
// 16 -> 64k
// 17 -> 128k
// 18 -> 256k
// 19 -> 512k
// 20 -> 1M
// 21 -> 2M
// 22 -> 4M
// 23 -> 8M

enum
{
	k_bits = 20,
	k_transform_len = 1 << k_bits
};

ISampleProducer* make_integer_upsampler(int up, double bw, ISampleProducer* input)
{
//	int filter_1_half_len = 200 + up * 200;
//	int filter_1_half_len = 640 + up * 460;
	int filter_1_half_len = 2000 + up * 2000;

	// clamp
	const int limit = 100000;
	if (filter_1_half_len > limit)
		filter_1_half_len = limit;

	int filter_1_len = filter_1_half_len * 2 + 1;
	int filter_2_len = filter_1_len * 2 + 1;

	// fixme this should go away
	int inbuf_len = k_transform_len - filter_2_len;
	printf("\nOSS transform-len=%dK, filter-len=%dK, inbuf-len=%dK\n", k_transform_len/1024, filter_2_len / 1024, inbuf_len / 1024);

	// don't even run if the inbuf is below N
//	assert(inbuf_len > 100000);
	assert(inbuf_len > 500);

	double* filter_kernel = (double*)pffftd_aligned_malloc(k_transform_len * sizeof(double)); // filter_kernel is only used in init
	memset(filter_kernel, 0, k_transform_len * sizeof(double));

	create_filter(filter_kernel, filter_1_len, filter_2_len, k_transform_len, bw, up );
	normalize_filter(filter_kernel, k_transform_len, up);

	ISampleProducer* upsampler = new StreamerUpT<k_transform_len>(filter_kernel, filter_2_len, up, input);
	pffftd_aligned_free(filter_kernel); // filter_kernel is only used in init

	return upsampler;
}









ISampleProducer* streamer_factory(ISampleProducer* input, int sr_out)
{
	int sr_in = input->get_sample_rate();

	// fixme bw can be input
	double bw = 0.999f;
	if (sr_out < sr_in)
		bw *= (double)sr_out / (double)sr_in;

	// try and find integer ratio (too high becomes very expensive)
	// should do 147 / 320 to suppoer 96 -> 44.1k
	for (int up = 1; up < 32; ++up)
	{
		for (int decimate = 1; decimate < 513; ++decimate)
		{
			if ( (sr_out * decimate) == (sr_in * up))
			{
				// found rational, very cool, avoids interpolation

				if (1 == decimate)
				{
					// if decimate with 1, just pass-through instead...
					printf("Integer Upsampling: up = %d\n", up);
					// found rational, very cool, avoids interpolation
					return make_integer_upsampler(up, bw, input);
				}
				else
				{
					printf("Rational: %d / %d\n", up, decimate);
					return new TrivialDecimator(decimate, make_integer_upsampler(up, bw, input));
				}

			}
		}
	}

	return new InterpolatedSampler(sr_out, make_integer_upsampler(32, bw, input));
}


