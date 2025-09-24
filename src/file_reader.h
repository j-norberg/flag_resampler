#ifndef FILE_READER_H
#define FILE_READER_H

#include <cstdint> // types

#include "sample_producer.h"

struct FileReaderImpl;

// Represents an audio - stream
class FileReader : public ISampleProducer
{
public:

	enum ErrorCode
	{
		eOk,
		eNoFile,
	};

private:
	FileReaderImpl* _impl;

	// a buffer
	float* _buf_interleaved = nullptr;
	bool _did_final_clear = false; // if buffer is fully cleared because we read past the end

	int _sample_rate = 44100;
	int _channel_count = 1;

	int64_t _buf_frame_size = 0;
	int64_t _total_frame_size = 0;
	int64_t _frame_index = 0;
	ErrorCode _error_code = eOk;

	FileReader(const FileReader&); // non construction-copyable
	FileReader& operator=(const FileReader&); // non copyable

	inline void get_next_internal(float* buf_interleaved)
	{
		if (_frame_index >= _buf_frame_size)
		{
			// if buffered file, read more and restart from start of buffer
			read_more_data_from_file();
			_frame_index = 0;
		}

		int64_t read_index = _frame_index * _channel_count;
		for (int i = 0; i < _channel_count; ++i)
		{
			buf_interleaved[i] = _buf_interleaved[read_index + i];
		}

		++_frame_index;
	}

	inline void skip_next_internal()
	{
		if (_frame_index >= _buf_frame_size)
		{
			read_more_data_from_file();
			_frame_index = 0;
		}

		++_frame_index;
	}

public:

	FileReader(const char* file_name); // loads a whole file into memory

	ErrorCode get_error_code() { return _error_code; };

	int get_sample_rate() { return _sample_rate; }
	int get_channel_count() { return _channel_count; }

	int64_t get_frame_count() { return _total_frame_size; }

	~FileReader();

	void read_more_data_from_file();

	inline void get_next(float* buf_interleaved, int64_t frame_count)
	{
		// fixme can in theory be faster
		// basically just a memcpy until the end of the buffer
		for (int64_t i = 0; i < frame_count; ++i)
		{
			get_next_internal(buf_interleaved);
			buf_interleaved += _channel_count;
		}
	}

	inline void skip_next(int64_t frame_count)
	{
		// fixme can in theory be faster
		// basically just a memcpy until the end of the buffer
		for (int64_t i = 0; i < frame_count; ++i)
		{
			skip_next_internal();
		}
	}

};

#endif
