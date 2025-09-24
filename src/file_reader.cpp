// implementation is in main

#define DR_WAV_IMPLEMENTATION
#include "dep/dr_wav.h"

#define DR_FLAC_IMPLEMENTATION
#include "dep/dr_flac.h"


#include "file_reader.h"

enum { k_reader_buf_frames = 4 * 1024 * 1024 };

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
	_buf_interleaved = new float[k_reader_buf_frames * _channel_count];
	_buf_frame_size = k_reader_buf_frames;

	read_more_data_from_file();
}


void FileReader::read_more_data_from_file()
{
	if (_did_final_clear)
		return;

	size_t samples_to_read = k_reader_buf_frames * _channel_count;

	size_t samples_read = 0;
	// read N samples from the file
	if (_impl->_wav)
		samples_read = drwav_read_f32(_impl->_wav, samples_to_read, _buf_interleaved);
	else if (_impl->_flac)
		samples_read = drflac_read_pcm_frames_f32(_impl->_flac, samples_to_read, _buf_interleaved);
	else
		_did_final_clear = true; // unexpected case

	// filled buffer, good
	if (samples_read == samples_to_read)
		return;

	if (samples_read == 0)
		_did_final_clear = true;

	// clear rest of buffer
	for (; samples_read < samples_to_read; ++samples_read)
		_buf_interleaved[samples_read] = 0;

}


FileReader::~FileReader()
{
	delete _impl;
	drwav_free((void*)_buf_interleaved);
}
