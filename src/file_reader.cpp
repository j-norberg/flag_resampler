// implementation is in main

#pragma warning( disable : 4514 )
#pragma warning( disable : 4820 ) // padding
#pragma warning( disable : 5045 ) // whole file

#pragma warning( push )

#pragma warning( disable : 4365 )
#pragma warning( disable : 4820 )
#pragma warning( disable : 5039 )
#pragma warning( disable : 5219 )
#pragma warning( disable : 6262 )

#define DR_WAV_IMPLEMENTATION
#include "dep/dr_wav.h"

#define DR_FLAC_IMPLEMENTATION
#include "dep/dr_flac.h"

#pragma warning( pop )

#include <math.h> // for fabs

#include "file_reader.h"

// this struct holds the dr_wav and dr_flac structs
struct FileReaderImpl
{
	drwav* _wav = nullptr;
	drflac* _flac = nullptr;

	~FileReaderImpl()
	{
		if (_wav != nullptr)
			drwav_close(_wav);

		if (_flac != nullptr)
			drflac_close(_flac);
	}
};

class BufReader
{
	const double* _buf = nullptr;
	int _stride = 1;
public:
	BufReader(const double* buf, int stride) : _buf(buf), _stride(stride) { }
	double read(int i) { return _buf[i * _stride]; }
};

class BufWriter
{
	double* _buf = nullptr;
	int _stride = 1;
public:
	BufWriter(double* buf, int stride) : _buf(buf) , _stride(stride) { }
	void write(int i, double v) { _buf[i * _stride] = v; }
};

struct ACValue
{
	double slope0 = 0;
	double slope1 = 0;
	double slope2 = 0;
	double ofs = 0;

	ACValue(BufReader& buf, int index)
	{
		double a = buf.read(index + 0);
		double b = buf.read(index + 1);
		double c = buf.read(index + 2);
		double d = buf.read(index + 3);

		ofs = a;
		slope0 = b - a;
		slope1 = c - b;
		slope2 = d - c;
	}

	// Score, lower is closer
	double Score(ACValue& v)
	{
		return
			fabs(slope0 - v.slope0) +
			fabs(slope1 - v.slope1) +
			fabs(slope2 - v.slope2) +
			fabs(ofs - v.ofs) * 0.1; // weigh the difference in sum less
	}
};



// fill_count = how many samples need to be filled (starting at 0), counting into future
// history_count = how many samples of history we have (starting at 0) counting to the past
void predict(BufWriter& predicted, BufReader& history, int fill_count, int history_count)
{
#if 0
	// force extending last sample
	{
		double v = history.read(0);
		for (int i = 0; i < fill_count; ++i)
			predicted.write(i, v);	// write

		return;
	}
#endif

	// limit based on a small value on a cd
	double limit = 9.0 / 65000.0;

	// if 3 last very quiet, do linear
	if (
		fabs(history.read(0)) < limit &&
		fabs(history.read(1)) < limit &&
		fabs(history.read(2)) < limit
		)
	{
//		puts("Linear prediction");
		double slope = history.read(0) - history.read(1);

		double v = history.read(0);
		for (int i = 0; i < fill_count; ++i)
		{
			slope *= 0.95;					// limit slope severely
			v += slope;						// also aim slow toward 0 (allow oscillation)
			slope -= v * 0.3;
			v *= 0.9;						// clamp v (pull to 0)
			predicted.write(i,v);	// write
		}
		return;
	}

//	puts("Autocorrelation");

	// try autocorrelation?
	// look at 4 samples care about the slopes a lot and the avg value a little
	// look in entire history and find best match
	// copy-paste the best match over and over
	ACValue edge(history, 0);
	int best_ind = 0;
	double best_score = 1000; //
	double best_ofs = 0; //
	for (int i = 4; i < history_count - 3; ++i)
	{
		ACValue v(history, i);
		double score = edge.Score(v);
		// printf("i=%d, score=%f\n", i, score);

		if (score < best_score)
		{
			best_ind = i;
			best_score = score;
			best_ofs = v.ofs;
		}
	}

//	printf("fill_count=%d, best_ind=%d best_ofs=%f edge_ofs=%f\n", fill_count, best_ind, best_ofs, edge.ofs);

	double lf_ofs = edge.ofs - best_ofs;
//	printf("lf_ofs=%f\n", lf_ofs);

	// copy all
	int read_i = best_ind - 1;
	for (int i = 0; i < fill_count; ++i)
	{
		double x = (i + 0.5) / (double)fill_count;
		double scale = 1.0 - x*x*(3.0-2.0*x);
		double v = (history.read(read_i) + lf_ofs) * scale;
		predicted.write(i, v);

		read_i -= 1;
		if (read_i < 0)
			read_i = best_ind - 1;
	}
}

void check_all(double* /* buf */, int /* count */)
{
#if false
	for (int i = 0; i < count; ++i)
	{
		double v = buf[i];

		// simple check
		if (v > 4.0)
			puts("bad hi");
		if (v < -4.0)
			puts("bad lo");
	}
#endif
}

FileReader::FileReader(const char* file_name)
{
	// assume no file
	_channel_count = 1;
	_total_frame_size = 1;
	_sample_rate = 44100;
	_error_code = eNoFile;

	_impl = new FileReaderImpl();

	// 1. try wav
	_impl->_wav = drwav_open_file(file_name);
	if (_impl->_wav != nullptr)
	{
		// get stats from wav
		_channel_count = _impl->_wav->channels;
		_total_frame_size = (int64_t)_impl->_wav->totalSampleCount / _channel_count;
		_sample_rate = (int)_impl->_wav->sampleRate;
		_error_code = eOk;
	}
	else
	{
		// 2. try flac
		_impl->_flac = drflac_open_file(file_name, nullptr);
		if (_impl->_flac != nullptr)
		{
			// get stats from wav
			_channel_count = _impl->_flac->channels;
			_total_frame_size = (int64_t)_impl->_flac->totalPCMFrameCount / _channel_count;
			_sample_rate = (int)_impl->_flac->sampleRate;
			_error_code = eOk;
		}
	}

	// allocate buffers
	_buf_interleaved_f64 = new double[(size_t)(k_buf_frames * _channel_count)];

	// set to zero to ensure no bad values start out in here (needed for very short files)
	memset(_buf_interleaved_f64, 0, k_buf_frames * _channel_count * sizeof(double));

	// check whole buffer
	check_all(_buf_interleaved_f64, k_buf_frames * _channel_count);

	// fill half of buffer and 
	read_data_from_file(k_half, k_half);

	// check whole buffer
	check_all(_buf_interleaved_f64, k_buf_frames * _channel_count);

	// predict each channel
	// buf is half-way through the whole buffer
	double* buf = _buf_interleaved_f64 + k_half * _channel_count;
	for (int i = 0; i < _channel_count; ++i)
	{
		BufWriter w(buf - _channel_count + i, -_channel_count);
		BufReader r(buf + i, _channel_count);
		predict(w, r, k_half, k_half);
	}

	// check whole buffer
	check_all(_buf_interleaved_f64, k_buf_frames * _channel_count);

}


// from float to double
void UpConvert(double* dest, const float* source, size_t count)
{
	for (size_t i = 0; i < count; ++i)
		dest[i] = (double)source[i];
}


void FileReader::read_data_from_file(size_t frame_offset, size_t frame_count)
{
	size_t frame_size = sizeof(double) * _channel_count;

	if (_buf_predicted_tail != nullptr)
	{
		// handle tail, copy tail or zero
		double* dst = _buf_interleaved_f64 + frame_offset * _channel_count;
		for (int i = 0 ; i < k_buf_frames; ++i)
		{
			if (_tail_frame_index < k_buf_frames)
				memcpy(dst, _buf_predicted_tail + _tail_frame_index * _channel_count, frame_size);
			else
				memset(dst, 0, frame_size); // zero
			dst += _channel_count;
			++_tail_frame_index;
		}

		// check both buffers
		check_all(_buf_predicted_tail, k_buf_frames * _channel_count);
		check_all(_buf_interleaved_f64, k_buf_frames * _channel_count);
		return;
	}

	size_t samples_offset = frame_offset * _channel_count;
	size_t samples_to_read = frame_count * _channel_count;
	size_t samples_read = 0;

	// write 32bit data to second half and explode later 
	double* buf_f64 = _buf_interleaved_f64 + samples_offset;
	float* buf_f32 = ((float*)buf_f64) + samples_to_read;

	// read N samples from the file
	if (_impl->_wav)
	{
		drwav* pWav = _impl->_wav;
		if (pWav->translatedFormatTag == DR_WAVE_FORMAT_IEEE_FLOAT && pWav->bytesPerSample == 8)
		{
			// 64bit proper
			samples_read = drwav_read(pWav, samples_to_read, buf_f64);
		}
		else
		{
			samples_read = drwav_read_f32(_impl->_wav, samples_to_read, buf_f32);
			UpConvert(buf_f64, buf_f32, samples_read);
		}
	}
	else if (_impl->_flac)
	{
		samples_read = drflac_read_pcm_frames_f32(_impl->_flac, samples_to_read, buf_f32);
		UpConvert(buf_f64, buf_f32, samples_read);
	}
	else
	{
		puts("did not expect to get here");
	}

	// check whole buffer
	check_all(_buf_interleaved_f64, k_buf_frames * _channel_count);

	// filled buffer, good, nothing more to do
	if (samples_read == samples_to_read)
		return;

	// create prediction here
	_buf_predicted_tail = new double[(size_t)(k_buf_frames * _channel_count)];

	// copy last samples (half buffer) to tail (first half)
	// read half the buffer to end up at the last read sample
	int samples_to_copy = k_half * _channel_count;
	for (int i = 0; i < samples_to_copy; ++i)
		_buf_predicted_tail[i] = _buf_interleaved_f64[(samples_read + i - samples_to_copy) & k_mask];

	// predict each channel
	// buf is half-way through the whole buffer
	double* buf = _buf_predicted_tail + samples_to_copy;
	for (int i = 0; i < _channel_count; ++i)
	{
		BufReader r(buf - _channel_count + i, -_channel_count);
		BufWriter w(buf + i, _channel_count);
		predict(w, r, k_half, k_half);
	}
	_tail_frame_index = k_half;

	// create prediction in other half
	// copy rest of buffer from predicted tail
	{
		size_t frames_read = samples_read / _channel_count;
		size_t frames_to_read = samples_to_read / _channel_count;
		double* dst = _buf_interleaved_f64 + frames_read * _channel_count; // write after
		for (; frames_read < frames_to_read; ++frames_read)
		{
			if (_tail_frame_index < k_buf_frames)
				memcpy(dst, _buf_predicted_tail + _tail_frame_index * _channel_count, frame_size);
			else
				memset(dst, 0, sizeof(double) * _channel_count); // zero
			++_tail_frame_index;
			dst += _channel_count;
		}
	}

	// double check tail
	check_all(_buf_predicted_tail, k_buf_frames * _channel_count);

}


FileReader::~FileReader()
{
	delete _impl;
	delete [] _buf_interleaved_f64;

	if (_buf_predicted_tail != nullptr)
		delete[] _buf_predicted_tail;
}
