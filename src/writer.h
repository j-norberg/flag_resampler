#ifndef WRITER_H
#define WRITER_H

#include <cstdint> // types

// implementation is in main
#include "dep/dr_wav.h"

class ISampleProducer;

// This is not a sample-producer
class Writer
{
public:

	enum OutFormat
	{
		eFmtFloat,
		eFmtInt24Dithered,
		eFmtInt16Dithered
	};

	enum State
	{
		eInit,
		eNoFile,
		eActive,
		eDone
	};

private:
	double* _buf_interleaved_f64 = nullptr;
	
	// for quantization (used for 16 and 24 bit)
	double* _buf_dither_noise = nullptr;

	// need to quantize from double to float too, but apply no dither noise
	uint8_t* _buf_interleaved_quantized = nullptr;

	int64_t _total_frame_count = 0;
	int64_t _frame_index = 0;

	State _state = eInit;
	drwav* _wav = nullptr;
	ISampleProducer* _sample_producer;

	OutFormat _fmt;

	Writer(const Writer&); // non construction-copyable
	Writer& operator=(const Writer&); // non copyable

public:

	Writer(const char* file_name, int64_t total_frame_count, ISampleProducer* input, OutFormat fmt);

	State get_error_code()
	{
		return _state;
	};

	~Writer();

	int get_progress_percent()
	{
		return (int)((_frame_index * 100) / _total_frame_count);
	};

	// converts one block of data
	bool update();

	// does all blocks of data
	void update_all();
};

extern bool simple_wav_write_mono(const char* file_name, float* buf, int64_t total_frame_count);


#endif
