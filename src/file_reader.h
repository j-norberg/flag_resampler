#ifndef FILE_READER_H
#define FILE_READER_H

#include <cstdint> // types

#include "sample_producer.h"

struct FileReaderImpl;

// Represents an audio - stream
class FileReader : public ISampleProducer
{
	enum {
		k_buf_frames = 4 * 1024, // the longer this is the longer the "prediction" is
		k_mask = k_buf_frames-1,
		k_half = k_buf_frames / 2
	};

public:

	enum ErrorCode
	{
		eOk,
		eNoFile,
	};

private:
	FileReaderImpl* _impl;

	double* _buf_interleaved_f64 = nullptr;	// a buffer

	// when the file is fully read we use the last samples and predict some more
	double* _buf_predicted_tail = nullptr;
	int64_t _tail_frame_index = 0; // points into the tail-buffer, after it's exhausted, fill w. zero

	int _sample_rate = 44100;
	int _channel_count = 1;

	int64_t _total_frame_size = 0;
	int64_t _frame_index = 0;
	ErrorCode _error_code = eOk;

	FileReader(const FileReader&); // non construction-copyable
	FileReader& operator=(const FileReader&); // non copyable

	inline void get_next_internal(double* buf_interleaved)
	{
		if (_frame_index >= k_buf_frames)
		{
			// if buffered file, read more and restart from start of buffer
			read_data_from_file(0, k_buf_frames);
			_frame_index = 0;
		}

		int64_t read_index = _frame_index * _channel_count;
		for (int i = 0; i < _channel_count; ++i)
		{
			double v = _buf_interleaved_f64[read_index + i];
			buf_interleaved[i] = v;
		}

		++_frame_index;
	}

	inline void skip_next_internal()
	{
		if (_frame_index >= k_buf_frames)
		{
			read_data_from_file(0, k_buf_frames);
			_frame_index = 0;
		}

		++_frame_index;
	}

public:

	FileReader(const char* file_name); // loads a whole file into memory

	ErrorCode get_error_code() { return _error_code; };
	int64_t get_frame_count() { return _total_frame_size; }

	~FileReader();

	void read_data_from_file(size_t frame_offset, size_t frame_count);

	int get_sample_rate() override { return _sample_rate; }
	int get_channel_count() override { return _channel_count; }

	virtual int get_padding_frame_count() override
	{
		return k_half; // always create half a buffer of prediction
	}

	inline void get_next(double* buf_interleaved, int64_t frame_count)
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
