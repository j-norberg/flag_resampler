// implementation is in main

#define DR_WAV_IMPLEMENTATION
#include "dep/dr_wav.h"

#define DR_FLAC_IMPLEMENTATION
#include "dep/dr_flac.h"


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
		_total_frame_size = _impl->_wav->totalSampleCount / _channel_count;
		_sample_rate = _impl->_wav->sampleRate;
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
			_total_frame_size = _impl->_flac->totalPCMFrameCount / _channel_count;
			_sample_rate = _impl->_flac->sampleRate;
			_error_code = eOk;
		}
	}

	// allocate buffers
	_buf_interleaved_f64 = new double[k_reader_buf_frames * _channel_count];
	_last_frame_f64 = new double[_channel_count];
	
	// clear this (pretend we read 0)
	for (int i = 0; i < _channel_count; ++i)
		_last_frame_f64[i] = 0;

	read_more_data_from_file();
}


void UpConvert(double* dest, const float* source, size_t count)
{
	for (size_t i = 0; i < count; ++i)
	{
		dest[i] = (double)source[i];
	}
}

void FileReader::read_more_data_from_file()
{
	if (_did_final_clear)
		return;

	size_t samples_to_read = k_reader_buf_frames * _channel_count;

	size_t samples_read = 0;

	// write 32bit data to second half and explode later 
	float* buf_f32 = ((float*)_buf_interleaved_f64) + samples_to_read;

	// read N samples from the file
	if (_impl->_wav)
	{
		drwav* pWav = _impl->_wav;
		if (pWav->translatedFormatTag == DR_WAVE_FORMAT_IEEE_FLOAT && pWav->bytesPerSample == 8)
		{
			samples_read = drwav_read(pWav, samples_to_read, _buf_interleaved_f64);
		}
		else
		{
			samples_read = drwav_read_f32(_impl->_wav, samples_to_read, buf_f32);
			UpConvert(_buf_interleaved_f64, buf_f32, samples_read);
		}
	}
	else if (_impl->_flac)
	{
		samples_read = drflac_read_pcm_frames_f32(_impl->_flac, samples_to_read, buf_f32);
		UpConvert(_buf_interleaved_f64, buf_f32, samples_read);
	}
	else
	{
		_did_final_clear = true; // unexpected case
	}

	// filled buffer, good
	if (samples_read == samples_to_read)
		return;

	if (samples_read == 0)
		_did_final_clear = true;
	else
	{
		if (samples_read < (size_t)_channel_count)
		{
			// some issue, don't update last frame
		}
		else
		{
			size_t last_frames = samples_read - _channel_count;
			
			// store off last samples read
			for (int i = 0; i < _channel_count; ++i)
				_last_frame_f64[i] = _buf_interleaved_f64[last_frames + i];
		}
	}

	// clear rest of buffer with last valid sample
	for (; samples_read < samples_to_read; ++samples_read)
		_buf_interleaved_f64[samples_read] = _last_frame_f64[samples_read % _channel_count];

}


FileReader::~FileReader()
{
	delete _impl;
	delete [] _buf_interleaved_f64;
	delete [] _last_frame_f64;
}
