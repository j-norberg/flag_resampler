
class Timer
{
public:
	Timer() :
		_beg(clock_::now())
	{
	}
	void reset()
	{
		_beg = clock_::now();
	}

	double elapsed_ms() const
	{
		std::chrono::duration<double, std::milli> ms = clock_::now() - _beg;
		return ms.count();
	}

private:
	typedef std::chrono::high_resolution_clock clock_;
	typedef std::chrono::duration<double, std::ratio<1> > second_;
	std::chrono::time_point<clock_> _beg;
};
