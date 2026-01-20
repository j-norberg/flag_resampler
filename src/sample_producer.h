#ifndef SAMPLE_PRODUCER_H
#define SAMPLE_PRODUCER_H

#include <cstdint> // types

class ISampleProducer
{
public:
	virtual ~ISampleProducer() {};

	virtual int get_sample_rate() = 0;
	virtual int get_channel_count() = 0;

	virtual void peek(double* buf_interleaved) = 0;
	virtual void get_next(double* buf_interleaved, int64_t frame_count) = 0;
	virtual void skip_next(int64_t frame_count) = 0;
};

#endif
