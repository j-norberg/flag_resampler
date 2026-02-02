#ifndef MEMORY_READER_H
#define MEMORY_READER_H

#include <cstdint> // types

#include "sample_producer.h"

// Represents an audio - stream
class MemoryReader : public ISampleProducer
{
private:
	const double* _buf_interleaved = nullptr;

	int _sample_rate = 44100;
	int _channel_count = 1;

	int64_t _frame_index = 0;
	int64_t _frame_count = 0;

	MemoryReader(const MemoryReader&); // non construction-copyable
	MemoryReader& operator=(const MemoryReader&); // non copyable

public:

	MemoryReader(const double* buf, int64_t frame_count, int sample_rate, int channel_count)
	{
		_buf_interleaved = buf;
		_channel_count = channel_count;
		_sample_rate = sample_rate;
		_frame_count = frame_count;
	}

	int get_sample_rate() override { return _sample_rate; }
	int get_channel_count() override { return _channel_count; }
	
	inline void get_next(double* buf_interleaved, int64_t frame_count) override 
	{
		// fixme can be simpler/faster
		int64_t fi = 0;
		int64_t sample_index = _frame_index * _channel_count;
		int64_t write_index = 0;
		for (; fi < frame_count; ++fi)
		{
			if (_frame_index >= _frame_count)
				break;

			for (int i = 0; i < _channel_count; ++i)
			{
				buf_interleaved[write_index] = _buf_interleaved[sample_index];
				++sample_index;
				++write_index;
			}

			++_frame_index;
		}

		for (; fi < frame_count; ++fi)
		{
			for (int i = 0; i < _channel_count; ++i)
			{
				buf_interleaved[write_index] = 0; // write 0 beyound end of buffer
				++write_index;
			}
		}
	}

	inline void skip_next(int64_t frame_count) override 
	{
		_frame_index += frame_count;
	}

};

// Represents an audio - stream
class MemoryReader32 : public ISampleProducer
{
private:
	const float* _buf_interleaved = nullptr;

	int _sample_rate = 44100;
	int _channel_count = 1;

	int64_t _frame_index = 0;
	int64_t _frame_count = 0;

	MemoryReader32(const MemoryReader&); // non construction-copyable
	MemoryReader32& operator=(const MemoryReader&); // non copyable

public:

	MemoryReader32(const float* buf, int64_t frame_count, int sample_rate, int channel_count)
	{
		_buf_interleaved = buf;
		_channel_count = channel_count;
		_sample_rate = sample_rate;
		_frame_count = frame_count;
	}

	int get_sample_rate() override { return _sample_rate; }
	int get_channel_count() override { return _channel_count; }
	
	inline void get_next(double* buf_interleaved, int64_t frame_count)
	{
		// fixme can be simpler/faster
		int64_t fi = 0;
		int64_t sample_index = _frame_index * _channel_count;
		int64_t write_index = 0;
		for (; fi < frame_count; ++fi)
		{
			if (_frame_index >= _frame_count)
				break;

			for (int i = 0; i < _channel_count; ++i)
			{
				buf_interleaved[write_index] = (double)_buf_interleaved[sample_index]; // cast here
				++sample_index;
				++write_index;
			}

			++_frame_index;
		}

		for (; fi < frame_count; ++fi)
		{
			for (int i = 0; i < _channel_count; ++i)
			{
				buf_interleaved[write_index] = 0; // write 0 beyound end of buffer
				++write_index;
			}
		}
	}

	inline void skip_next(int64_t frame_count)
	{
		_frame_index += frame_count;
	}

};

#endif
