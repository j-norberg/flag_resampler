
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

