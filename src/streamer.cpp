#include "dep/pffft-double.h"

#include <vector>
#include <cassert>
#include <cstring> // for memset
#include <string> // for debugging

#define _USE_MATH_DEFINES
#include <math.h> // for M_PI

#include <limits> // for epsilon

#pragma warning( disable : 4514 ) // unused inl. removed
#pragma warning( disable : 4820 ) // padding
#pragma warning( disable : 5045 ) // spectre

#include "ok_interp.inl"
#include "streamer.h"
#include "sample_producer.h"

#include "writer.h" // for simple write mono


// hamming is bad
// Blackman–Harris window
// Dolph–Chebyshev window
// Generalized adaptive polynomial (GAP) window

#if 0
// blackman
inline double win_function(double p)
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
inline double win_function(double p)
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
		double win = win_function(p01) * window_amp_scale;
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

	void init_oss(double* temp_windowed_sinc, int filter_size)
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

	// 1. fill buffer (call this _inblock_size times)
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

		// multiply with the freq-version of the kernel
		pffftd_zconvolve_no_accumulate(_fft, _bufs->_in_buf_freq, _bufs->_buf_filter, _bufs->_out_buf_freq, _transform_recip);

		// IFFT
		pffftd_transform(_fft, _bufs->_out_buf_freq, _bufs->_out_buf_time, _bufs->_buf_work, PFFFTD_BACKWARD);
	}

	// 3. read from here (call this _inblock_size times)
	inline double read(int index)
	{
		return _bufs->_out_buf_time[_filter_size + index];
	}
};








//
template<int BITS>
struct StreamerUpT : ISampleProducer
{
	enum { TRANSFORM_SIZE = 1 << BITS };
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
		// fixme, share the freq-space-filter-kernel
		_oss = new OverlapSaveStruct[(size_t)_channel_count];
		for (int i = 0; i < _channel_count; ++i)
			_oss[i].init_oss(filter_kernel, filter_size);

		_channels_buf.resize((size_t)_channel_count);
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
				for (size_t c = 0; c < (size_t)_channel_count; ++c)
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

		_channel_buffer.resize((size_t)_channel_count);

		_bufs = new Buf[(size_t)_channel_count];

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

	void internal_fill_buffers(size_t ofs)
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
			for (size_t c = 0; c < (size_t)_channel_count; ++c)
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
		for (int64_t frame_index = 0; frame_index < frame_count; ++frame_index)
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
		for (int64_t frame_index = 0; frame_index < frame_count; ++frame_index)
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

// maybe normalize each phase?
// or just the whole kernel
void normalize_filter(double* filter_kernel, int filter_len, int up)
{
#if 0
	// ensure symmetry
	{
		int a = 0;
		int b = filter_len - 1;
		for (; a < b; ++a, --b)
		{
			double diff = (filter_kernel[a] - filter_kernel[b]);
			if (fabs(diff) > 0.000001)
				puts("no symmetry, internal bug");

//			double avg = (filter_kernel[a] + filter_kernel[b]) * 0.5;
	//		filter_kernel[a] = avg;
		//	filter_kernel[b] = avg;
		}
	}
#endif

#if 1
	for (int phase = 0; phase < up; ++phase)
	{
		double sum = 0.0;
		for (int i = phase; i < filter_len; i += up)
			sum += filter_kernel[i];

		double mul = 1.0 / sum;

		//uint64_t* mul_p = (uint64_t*)&mul;
		//printf("phase=%i/%i, mul=%g %llx\n", phase, up, mul, *mul_p);

		for (int i = phase; i < filter_len; i += up)
			filter_kernel[i] *= mul;
	}
#endif

#if 0
	// finally full normalize
	{
		double sum = 0.0;
		for (int i = 0; i < filter_len; ++i)
			sum += filter_kernel[i];

		double mul = (double)up / sum;

		//		uint64_t* mul_p = (uint64_t*)&mul;
		//		printf("FINAL, mul=%g %llx\n", mul, *mul_p);

		for (int i = 0; i < filter_len; ++i)
			filter_kernel[i] *= mul;
	}
#endif

}

void create_filter_self_convolved(double* out_kernel_buffer, int len1, int transform_size, double bw, double up)
{
	PFFFTD_Setup* fft = pffftd_new_setup(transform_size, PFFFTD_REAL);
	double* buf_work = (double*)pffftd_aligned_malloc(transform_size * sizeof(double));

	// create initial filters
	double* filter1 = (double*)pffftd_aligned_malloc(transform_size * sizeof(double));
	memset(filter1, 0, transform_size * sizeof(double));
	double scale_len = bw / up;
	double scale_amp = bw / sqrt(up);
	create_windowed_sinc(filter1, len1, scale_len, scale_amp);
	double* buf_freq1 = (double*)pffftd_aligned_malloc(transform_size * sizeof(double));
	pffftd_transform(fft, filter1, buf_freq1, buf_work, PFFFTD_FORWARD);

	// create a second sligtly lower (maybe switch back to pure self convolve)
	double* filter2 = (double*)pffftd_aligned_malloc(transform_size * sizeof(double));
	memset(filter2, 0, transform_size * sizeof(double));
	scale_len = (bw * 0.999) / up;
	scale_amp = (bw * 0.999) / sqrt(up);
	create_windowed_sinc(filter2, len1, scale_len, scale_amp);
	double* buf_freq2 = (double*)pffftd_aligned_malloc(transform_size * sizeof(double));
	pffftd_transform(fft, filter2, buf_freq2, buf_work, PFFFTD_FORWARD);

	// self convolve
	double* buf_freq3 = (double*)pffftd_aligned_malloc(transform_size * sizeof(double));
	double transform_recip = 1.0 / transform_size;
	pffftd_zconvolve_no_accumulate(fft, buf_freq1, buf_freq2, buf_freq3, transform_recip); // convolve
	pffftd_transform(fft, buf_freq3, out_kernel_buffer, buf_work, PFFFTD_BACKWARD);

	pffftd_aligned_free(buf_freq3);
	pffftd_aligned_free(buf_freq2);
	pffftd_aligned_free(buf_freq1);
	pffftd_aligned_free(filter2);
	pffftd_aligned_free(filter1);
	pffftd_aligned_free(buf_work);
	pffftd_destroy_setup(fft);
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

ISampleProducer* make_integer_upsampler(int up, double bw, ISampleProducer* input, int quality_percentage, int limit = 100000)
{
	if (quality_percentage < 1)
		puts("ERROR");

	if (quality_percentage > 100)
		puts("ERROR");

	// formula to deal with BW-limiting and upsampling
//	int filter_1_half_len = 2 + up * 1;
//	int filter_1_half_len = 160 + up * 115;

	int filter_1_half_len = 320 + up * 230;

//	int filter_1_half_len = 640 + up * 460;

	// filter_1_half_len = 1280 + up * 920; // no amplification at nyqvist

	if (quality_percentage > 95)
		filter_1_half_len = 3200 + up * 2300;
	

	//	int filter_1_half_len = 2560 + up * 1840;


	// clamp
	if (filter_1_half_len > limit)
	{
//		printf("limit half-len to %d\n", limit);
		filter_1_half_len = limit;
	}

	int filter_1_len = filter_1_half_len * 2 + 1;

	// IMPORTANT: the length of self-convolving?
	int filter_2_len = filter_1_len * 2 - 1;

	// create kernel in a temp-location
	int bits = 6;
	int desired_transform_len = filter_2_len * 4; // fixme how much larger do we need?
	while ((1 << bits) < desired_transform_len)
		++bits;
	int transform_len = 1 << bits;

	int inbuf_len = transform_len - filter_2_len;
	printf("FFT bits=%d, transform-len=%dK, filter-len=%dK, inbuf-len=%dK, BW=%f\n", bits, transform_len / 1024, filter_2_len / 1024, inbuf_len / 1024, bw);
	fflush(stdout);

	double* filter_kernel = (double*)pffftd_aligned_malloc(transform_len * sizeof(double)); // filter_kernel is only used in init
	memset(filter_kernel, 0, transform_len * sizeof(double));
	create_filter_self_convolved(filter_kernel, filter_1_len, transform_len, bw, up);
	
//	simple_wav_write_mono_f64("kernel2.wav", filter_kernel, transform_len);

	normalize_filter(filter_kernel, filter_2_len, up);
	
	// how large transform do we really need?

	ISampleProducer* upsampler = nullptr;
	switch (bits)
	{
	case  6: upsampler = new StreamerUpT< 6>(filter_kernel, filter_2_len, up, input); break;
	case  7: upsampler = new StreamerUpT< 7>(filter_kernel, filter_2_len, up, input); break;
	case  8: upsampler = new StreamerUpT< 8>(filter_kernel, filter_2_len, up, input); break;
	case  9: upsampler = new StreamerUpT< 9>(filter_kernel, filter_2_len, up, input); break;
	case 10: upsampler = new StreamerUpT<10>(filter_kernel, filter_2_len, up, input); break;
	case 11: upsampler = new StreamerUpT<11>(filter_kernel, filter_2_len, up, input); break;
	case 12: upsampler = new StreamerUpT<12>(filter_kernel, filter_2_len, up, input); break;
	case 13: upsampler = new StreamerUpT<13>(filter_kernel, filter_2_len, up, input); break;
	case 14: upsampler = new StreamerUpT<14>(filter_kernel, filter_2_len, up, input); break;
	case 15: upsampler = new StreamerUpT<15>(filter_kernel, filter_2_len, up, input); break;
	case 16: upsampler = new StreamerUpT<16>(filter_kernel, filter_2_len, up, input); break;
	case 17: upsampler = new StreamerUpT<17>(filter_kernel, filter_2_len, up, input); break;
	case 18: upsampler = new StreamerUpT<18>(filter_kernel, filter_2_len, up, input); break;
	case 19: upsampler = new StreamerUpT<19>(filter_kernel, filter_2_len, up, input); break;
	case 20: upsampler = new StreamerUpT<20>(filter_kernel, filter_2_len, up, input); break;
	case 21: upsampler = new StreamerUpT<21>(filter_kernel, filter_2_len, up, input); break;
	case 22: upsampler = new StreamerUpT<22>(filter_kernel, filter_2_len, up, input); break;
	case 23: upsampler = new StreamerUpT<23>(filter_kernel, filter_2_len, up, input); break;

	default:
		printf("requested %d bits of FFT, too much\n", bits);
		fflush(stdout);
		assert(false);
		break;
	}

	pffftd_aligned_free(filter_kernel); // filter_kernel is only used in init

	return upsampler;
}








ISampleProducer* make_upsampler_pair(int up1, int up2, double bw, ISampleProducer* input, int quality_percentage)
{
	printf("Up-pair of %dx(HQ) and %dx\n", up1, up2);
	int limit = 600;
	if (quality_percentage < 90)
		limit = 100;

	return make_integer_upsampler(up2, 0.99, make_integer_upsampler(up1, bw, input, quality_percentage), quality_percentage, limit);
}

ISampleProducer* make_upsampler_chain(int up, double bw, ISampleProducer* input, int quality_percentage)
{
	// special case 1
	if (up < 2)
		return make_integer_upsampler(1, bw, input, quality_percentage);

	// try split in 2
	int factors[] = { 2,3,5,7 };
	int fcount = sizeof(factors) / sizeof(int);
	for (int i = 0 ; i < fcount; ++i)
	{
		// find the lowest factor and use HQ for that
		int f = factors[i];
		if (f >= up)
			break;

		if ((up % f) == 0)
			return make_upsampler_pair(f, up / f, bw, input, quality_percentage);
	}

	// need to use HQ for whole thing
	return make_integer_upsampler(up, bw, input, quality_percentage);
}

ISampleProducer* streamer_factory(ISampleProducer* input, int sr_out, int quality_percentage)
{
	int sr_in = input->get_sample_rate();

	// fixme bw can be input
	double bw = 0.999f;
	
	if (quality_percentage < 90)
		bw = 0.995;

	if (sr_out < sr_in)
		bw *= (double)sr_out / (double)sr_in;

	// try and find integer ratio (too high becomes very expensive)
	for (int up = 1; up <= 32; ++up)
	{
		for (int decimate = 1; decimate < 513; ++decimate)
		{
			if ( (sr_out * decimate) == (sr_in * up))
			{
				// found rational, very cool, avoids interpolation
				ISampleProducer* upsampler = make_upsampler_chain(up, bw, input, quality_percentage);
				if (1 == decimate)
				{
					// if decimate with 1, just pass-through instead...
					printf("Integer Upsampling: up = %d\n", up);
					return upsampler;
				}

				printf("Rational: %d / %d\n", up, decimate);
				return new TrivialDecimator(decimate, upsampler);
			}
		}
	}

	// interpolated
	printf("Not rational enough: %d / %d\n", sr_out, sr_in);
	return new InterpolatedSampler(sr_out, make_upsampler_chain(32, bw, input, quality_percentage));
}


