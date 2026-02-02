struct InterpolatedSampler : ISampleProducer
{
	ISampleProducer* _input;
	int _channel_count;
	int _out_sr;

	enum {
		k_bufs = 64, // interpolation only uses 6 samples 
		k_mask = k_bufs - 1,
		k_half = k_bufs >> 1
	};

	struct Buf {
		double _b[k_bufs] = { -0.99 }; // fixme should never be visited
	};

	std::vector<double> _channel_buffer;

	Buf* _bufs;

	int _read_index;
	int _read_index_frac;
	int _read_index_frac_add;
	int _read_index_frac_limit;
	double _read_index_frac_limit_recip;

	InterpolatedSampler(int sr, ISampleProducer* input)
	{
		_input = input;
		_channel_count = input->get_channel_count();
		_out_sr = sr;

		_read_index = 0;

		_read_index_frac = 0;
		_read_index_frac_add = input->get_sample_rate();
		_read_index_frac_limit = sr;
		_read_index_frac_limit_recip = 1.0 / (double)sr;

		_channel_buffer.resize((size_t)_channel_count);

		_bufs = new Buf[(size_t)_channel_count];

		// note: induced 2 samples of latency
		// but skip all up to that point
		int padding = input->get_padding_frame_count();
		padding -= 2;
		input->skip_next(padding);

		// fill first half of buffer with actual data
		internal_fill_buffers(0);
		internal_fill_buffers(k_half);
	}

	~InterpolatedSampler()
	{
		delete[] _bufs;
		delete _input;
	}

	int get_sample_rate() override { return _out_sr; };
	int get_channel_count() override { return _channel_count; }

	void internal_fill_buffers(size_t ofs)
	{
		if (_channel_count < 2)
		{
			// no need for de-interleave
			_input->get_next(_bufs[0]._b + ofs, k_half);
			return;
		}

		// de-interleave
		for (int i = 0; i < k_half; ++i)
		{
			_input->get_next(_channel_buffer.data(), 1);
			for (size_t c = 0; c < (size_t)_channel_count; ++c)
				_bufs[c]._b[ofs + i] = _channel_buffer[c];
		}
	}

	int get_padding_frame_count() override
	{
		// the padding is eaten in constructor and return 0 here
		return 0;
	}

	void get_next(double* buf_interleaved, int64_t frame_count) override
	{
		int write_index = 0;
		for (int64_t frame_index = 0; frame_index < frame_count; ++frame_index)
		{
			// calculate frac (common for all channels)
			double frac = (double)_read_index_frac * _read_index_frac_limit_recip;

			// write interleaved
			for (int c = 0; c < _channel_count; ++c)
			{
				buf_interleaved[write_index] = sample_32x_6p_5z(_bufs[c]._b, k_mask, _read_index, frac);
				++write_index;
			}

			// move index
			_read_index_frac += _read_index_frac_add;

			// handle wrapping (maybe trust modulo here)
			while (_read_index_frac >= _read_index_frac_limit)
			{
				++_read_index;
				_read_index_frac -= _read_index_frac_limit;

				switch (_read_index)
				{
				case k_half:
					internal_fill_buffers(0);
					break;
				case k_bufs:
					internal_fill_buffers(k_half);
					_read_index = 0;
					break;
				}
			}

		}
	}

	void skip_next(int64_t frame_count) override
	{
		for (int64_t frame_index = 0; frame_index < frame_count; ++frame_index)
		{
			// move index
			_read_index_frac += _read_index_frac_add;

			// handle wrapping (maybe trust modulo here)
			while (_read_index_frac >= _read_index_frac_limit)
			{
				++_read_index;
				_read_index_frac -= _read_index_frac_limit;

				switch (_read_index)
				{
				case k_half:
					internal_fill_buffers(0);
					break;
				case k_bufs:
					internal_fill_buffers(k_half);
					_read_index = 0;
					break;
				}
			}

		}
	}
};
