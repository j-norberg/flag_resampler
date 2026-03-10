#pragma warning( disable : 4514 ) // unused inl. removed
#pragma warning( disable : 4710 ) // not inlined
#pragma warning( disable : 4711 ) // inline expansion
#pragma warning( disable : 4820 ) // padding
#pragma warning( disable : 5045 ) // spectre

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

	// to ensure no aliasing (of memory), all buffers have to be explicit
	double _buf_filter[K_BUFS];
	double _buf_work[K_BUFS];

	double _in_buf_freq[K_BUFS];	// input is transformed into frequences

	double _out_buf_freq[K_BUFS];	// result of convolve goes here
	double _out_buf_time[K_BUFS];	// IFFT here
};

// 
template<int TRANSFORM_SIZE>
struct OverlapSaveStructT : Job
{
	typedef BufsT<TRANSFORM_SIZE> Bufs;

	int _filter_size = 0;
	int _inblock_size = 0;
	double _transform_recip = 1.0;

	PFFFTD_Setup* _fft; // could be shared btw all instances

	Bufs* _bufs;

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

	void init_oss(const double* temp_windowed_sinc, int filter_size)
	{
		assert(_bufs == nullptr);
		assert(_fft == nullptr);
		assert((filter_size & 1) == 1); // require odd filter length

		_filter_size = filter_size;
		_inblock_size = TRANSFORM_SIZE - filter_size;
		_transform_recip = 1.0 / TRANSFORM_SIZE; // need to be double here

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

	void init_from_previous(Bufs* prev)
	{
		// scoot the input-buffer (or copy from prev)
		memmove(_bufs->_in_buf_time, prev->_in_buf_time + _inblock_size, _filter_size * sizeof(double));
	}

	virtual void run()
	{
		// FFT
		pffftd_transform(_fft, _bufs->_in_buf_time, _bufs->_in_buf_freq, _bufs->_buf_work, PFFFTD_FORWARD);

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

	// a block is including all channels
	struct Block
	{
		Block(int channel_count, const double* filter_kernel, int filter_size )
		{
			_channel_count = channel_count;
			_channel = new OverlapSaveStruct[(size_t)channel_count];

			for (int i = 0; i < _channel_count; ++i)
				_channel[i].init_oss(filter_kernel, filter_size);
		}

		~Block()
		{
			delete[] _channel;
		}

		void init_from_previous(Block* prev)
		{
			for (int i = 0; i < _channel_count; ++i)
				_channel[i].init_from_previous(prev->_channel[i]._bufs);
		}

		// are all jobs in this block done?
		bool all_channels_done()
		{
			for (int i = 0; i < _channel_count; ++i)
			{
				if (!_channel[i].done)
					return false;
			}
			return true;
		}

		void give_to_workers_or_run(WorkQueue* work_queue)
		{
			if (work_queue != nullptr)
			{
				// give to workers
				for (int i = 0; i < _channel_count; ++i)
					work_queue->add(&_channel[i]);
			}
			else
			{
				// work instantly
				for (int i = 0; i < _channel_count; ++i)
				{
					_channel[i].run();
					_channel[i].done = true;
				}
			}
		}

		OverlapSaveStruct* _channel; // one per channel
		int _channel_count;
	};

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

	int _rotations = 0; // just for debugging

	// these are queued-up and worked on by threads
	// the front one is the current one, the last one the last one queued
	std::queue<Block*> _queue;

	WorkQueue* _work_queue;

	StreamerUpT(WorkQueue* work_queue, const double* filter_kernel, int filter_size, int up, ISampleProducer* input)
	{
		_work_queue = work_queue;

		_filter_size = filter_size;
		_inblock_size = TRANSFORM_SIZE - filter_size;
		_up = up;

		_pad_ttl = 0;
		_out_index = _inblock_size; // force initial feed to buffer

		_input = input;
		_channel_count = input->get_channel_count();
		_channels_buf.resize((size_t)_channel_count);

		_sr = _input->get_sample_rate() * up;

		// build up the queue, how long do we need?
		const int in_flight = 24;
		for (int i = 0; i < in_flight; ++i)
		{
			Block* block = new Block(_channel_count, filter_kernel, filter_size);
			if (!_queue.empty())
				block->init_from_previous(_queue.back()); // copy shifted data from previous

			fill_block_from_input(block);
			block->give_to_workers_or_run(_work_queue);
			_queue.push(block);
		}
	}

	~StreamerUpT()
	{
		// delete all blocks
		while (!_queue.empty())
		{
			delete _queue.front();
			_queue.pop();
		}

		delete _input;
	}

	int get_sample_rate() override { return _sr; };
	int get_channel_count() override { return _channel_count; }

	inline void fill_block_from_input(Block* b)
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
					b->_channel[c].write(i, _channels_buf[c]);

				// reset ttl
				_pad_ttl = _up;
			}
			else
			{
				// zero-pad
				for (int c = 0; c < _channel_count; ++c)
					b->_channel[c].write(i, 0);
			}

			--_pad_ttl;
		}
	}

	virtual int get_padding_frame_count() override
	{
		int external_upsampled_padding = _input->get_padding_frame_count() * _up;
		int internal_padding = (_filter_size - 1) / 2;
		return external_upsampled_padding + internal_padding;
	}

	// call this when the front of the queue is all read
	// and ready to be reused
	// return current
	Block* rotate_queue()
	{
		++_rotations;
//		printf("rotate %d\n", _rotations); fflush(stdout);
		Block* current_block = _queue.front();
		current_block->init_from_previous(_queue.back()); // copy shifted data from previous
		fill_block_from_input(current_block);

		current_block->give_to_workers_or_run(_work_queue);

		_queue.push(current_block);
		_queue.pop(); // remove current from queue

		// get the next one
		current_block = _queue.front();
		while (!current_block->all_channels_done())
		{
			// if the main thread would wait on the workers, do work in main thread too
			Job* j = _work_queue->try_take();
			if (j != nullptr)
			{
				// puts("main-working"); fflush(stdout);
				j->run();
				j->done = true;
			}
			else {
				puts("main-waiting"); fflush(stdout);
				std::this_thread::yield();
			}
		}
		return current_block;
	}

	void get_next(double* buf_interleaved, int64_t frame_count) override
	{
		Block* current_block = _queue.front();

		int write_index = 0;
		for (int i = 0; i < frame_count; ++i)
		{
			if (_out_index >= _inblock_size)
			{
				current_block = rotate_queue();
				_out_index = 0;
			}

			for (int c = 0; c < _channel_count; ++c)
			{
				buf_interleaved[write_index] = current_block->_channel[c].read(_out_index); // current-block
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
				Block* current_block = rotate_queue();
				_out_index = 0;
			}

			++_out_index;
		}
	}

};

#include "trivial_decimator.inl"
#include "interpolated_sampler.inl"



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

#if 1
	// finally full normalize
	{
		double sum = 0.0;
		for (int i = 0; i < filter_len; ++i)
			sum += filter_kernel[i];

		double mul = (double)up / sum;

		// uint64_t* mul_p = (uint64_t*)&mul;
		// printf("FINAL, mul=%g %llx\n", mul, *mul_p);

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

	double* buf_freq3 = (double*)pffftd_aligned_malloc(transform_size * sizeof(double));
	double transform_recip = 1.0 / transform_size;

	// self convolve
	pffftd_zconvolve_no_accumulate(fft, buf_freq1, buf_freq1, buf_freq3, transform_recip); // convolve
	pffftd_transform(fft, buf_freq3, out_kernel_buffer, buf_work, PFFFTD_BACKWARD);

	pffftd_aligned_free(buf_freq3);
	pffftd_aligned_free(buf_freq1);
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

ISampleProducer* make_integer_upsampler(WorkQueue* workQueue, int up, double bw, ISampleProducer* input, int quality_percentage, int limit = 100000)
{
	if (quality_percentage < 1)
		puts("ERROR");

	if (quality_percentage > 100)
		puts("ERROR");

	// formula to deal with BW-limiting and upsampling
	int filter_1_half_len = 500 + up * 500;
	// filter_1_half_len = 1000 + up * 1000; // no amplification at nyqvist
	// filter_1_half_len = 2000 + up * 2000;

	if (quality_percentage > 95)
		filter_1_half_len = 3000 + up * 3000;

	// filter_1_half_len = 4000 + up * 4000;


	// clamp
	if (filter_1_half_len > limit)
	{
		printf("limit half-len to %d\n", limit);
		filter_1_half_len = limit;
	}

	int filter_1_len = filter_1_half_len * 2 + 1;

	// IMPORTANT: the length of self-convolving?
	int filter_2_len = filter_1_len * 2 - 1;

	// create kernel in a temp-location
	int bits = 6;
	int desired_transform_len = filter_2_len * 3; // fixme how much larger do we need?
	while ((1 << bits) < desired_transform_len)
		++bits;
	int transform_len = 1 << bits;

	int inbuf_len = transform_len - filter_2_len;
	printf("FFT bits=%d, transform-len=%dK, filter-len=%dK, inbuf-len=%dK, BW=%f\n", bits, transform_len / 1024, filter_2_len / 1024, inbuf_len / 1024, bw);
	fflush(stdout);

	double* filter_kernel = (double*)pffftd_aligned_malloc(transform_len * sizeof(double)); // filter_kernel is only used in init
	memset(filter_kernel, 0, transform_len * sizeof(double));
	create_filter_self_convolved(filter_kernel, filter_1_len, transform_len, bw, up);
	
	//
//	if (up < 3)
//	{
		// only first
//		simple_wav_write_mono_f64("kernel_self_convolved_t.wav", filter_kernel, transform_len);
//		simple_wav_write_mono_f64("kernel_self_convolved.wav", filter_kernel, filter_2_len);
//	}
	
	normalize_filter(filter_kernel, filter_2_len, up);
	
	// if up is large, might be better off doing
	// naive convolution

	// how large transform do we really need?

	ISampleProducer* upsampler = nullptr;
	switch (bits)
	{
	case  6: upsampler = new StreamerUpT< 6>(workQueue, filter_kernel, filter_2_len, up, input); break;
	case  7: upsampler = new StreamerUpT< 7>(workQueue, filter_kernel, filter_2_len, up, input); break;
	case  8: upsampler = new StreamerUpT< 8>(workQueue, filter_kernel, filter_2_len, up, input); break;
	case  9: upsampler = new StreamerUpT< 9>(workQueue, filter_kernel, filter_2_len, up, input); break;
	case 10: upsampler = new StreamerUpT<10>(workQueue, filter_kernel, filter_2_len, up, input); break;
	case 11: upsampler = new StreamerUpT<11>(workQueue, filter_kernel, filter_2_len, up, input); break;
	case 12: upsampler = new StreamerUpT<12>(workQueue, filter_kernel, filter_2_len, up, input); break;
	case 13: upsampler = new StreamerUpT<13>(workQueue, filter_kernel, filter_2_len, up, input); break;
	case 14: upsampler = new StreamerUpT<14>(workQueue, filter_kernel, filter_2_len, up, input); break;
	case 15: upsampler = new StreamerUpT<15>(workQueue, filter_kernel, filter_2_len, up, input); break;
	case 16: upsampler = new StreamerUpT<16>(workQueue, filter_kernel, filter_2_len, up, input); break;
	case 17: upsampler = new StreamerUpT<17>(workQueue, filter_kernel, filter_2_len, up, input); break;
	case 18: upsampler = new StreamerUpT<18>(workQueue, filter_kernel, filter_2_len, up, input); break;
	case 19: upsampler = new StreamerUpT<19>(workQueue, filter_kernel, filter_2_len, up, input); break;
	case 20: upsampler = new StreamerUpT<20>(workQueue, filter_kernel, filter_2_len, up, input); break;
	case 21: upsampler = new StreamerUpT<21>(workQueue, filter_kernel, filter_2_len, up, input); break;
	case 22: upsampler = new StreamerUpT<22>(workQueue, filter_kernel, filter_2_len, up, input); break;
	case 23: upsampler = new StreamerUpT<23>(workQueue, filter_kernel, filter_2_len, up, input); break;

	default:
		printf("requested %d bits of FFT, too much\n", bits);
		fflush(stdout);
		assert(false);
		break;
	}

	pffftd_aligned_free(filter_kernel); // filter_kernel is only used in init

	return upsampler;
}








ISampleProducer* make_upsampler_pair(WorkQueue* workQueue, int up1, int up2, double bw, ISampleProducer* input, int quality_percentage)
{
	printf("Up-pair of %dx(HQ) and %dx\n", up1, up2);
	int limit = 600;
	if (quality_percentage < 90)
		limit = 80;

	ISampleProducer* step1 = make_integer_upsampler(workQueue, up1, bw, input, quality_percentage);
	ISampleProducer* step2 = make_integer_upsampler(workQueue, up2, 0.75, step1, quality_percentage, limit);
	return step2;
}

ISampleProducer* make_upsampler_chain(WorkQueue* workQueue, int up, double bw, ISampleProducer* input, int quality_percentage)
{
	// try find pairs for "large" ish factors
	if (up > 3)
	{
		// try split in factors
		int factors[] = { 2,3,5,7 };
		int fcount = sizeof(factors) / sizeof(int);
		for (int i = 0; i < fcount; ++i)
		{
			// find the lowest factor and use HQ for that
			int f = factors[i];
			if (f >= up)
				break;

			if ((up % f) == 0)
				return make_upsampler_pair(workQueue, f, up / f, bw, input, quality_percentage);
		}
	}

	// need to use HQ for whole thing
	return make_integer_upsampler(workQueue, up, bw, input, quality_percentage);
}

ISampleProducer* streamer_factory(WorkQueue* workQueue, ISampleProducer* input, int sr_out, int quality_percentage)
{
	int sr_in = input->get_sample_rate();

	// fixme bw can be input
	// but currently is kind of connected to the filter-size
	double bw = 0.9992;
	
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
				ISampleProducer* upsampler = make_upsampler_chain(workQueue, up, bw, input, quality_percentage);
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
	printf("Not rational enough: %d / %d (will go through 32x)\n", sr_out, sr_in);
	ISampleProducer* upsampler32 = make_upsampler_chain(workQueue, 32, bw, input, quality_percentage);
	return new InterpolatedSampler(sr_out, upsampler32);
}


