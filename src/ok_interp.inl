// These constants from paper:

// "Polynomial Interpolators for
// High-Quality Resampling of
// Oversampled Audio"

// Olli Niemitalo in October 2001.

// only save the HQ one. Also changed types to double
inline double sample_16x_6p_5z(double* buf, int mask, int index, double frac)
{
	double y0 = buf[(index + 0) & mask];
	double y1 = buf[(index + 1) & mask];
	double y2 = buf[(index + 2) & mask];
	double y3 = buf[(index + 3) & mask];
	double y4 = buf[(index + 4) & mask];
	double y5 = buf[(index + 5) & mask];

	// Optimal 16x (6-point, 5th-order) (z-form)
	double z = frac - 0.5;
	double even1 = y3 + y2, odd1 = y3 - y2;
	double even2 = y4 + y1, odd2 = y4 - y1;
	double even3 = y5 + y0, odd3 = y5 - y0;

	double c0 = even1 * 0.41809989254549901 + even2 * 0.08049339946273310 + even3 * 0.00140670799165932;
	double c1 = odd1 * 0.32767596257424964 + odd2 * 0.20978189376640677 + odd3 * 0.00859567104974701;
	double c2 = even1 * -0.206944618112960001 + even2 * 0.18541689550861262 + even3 * 0.02152772260740132;
	double c3 = odd1 * -0.21686095413034051 + odd2 * 0.02509557922091643 + odd3 * 0.02831484751363800;
	double c4 = even1 * 0.04163046817137675 + even2 * -0.06244556931623735 + even3 * 0.02081510113314315;
	double c5 = odd1 * 0.07990500783668089 + odd2 * -0.03994519162531633 + odd3 * 0.00798609327859495;

	return ((((c5 * z + c4) * z + c3) * z + c2) * z + c1) * z + c0;
}

inline double sample_32x_6p_5z(double* buf, int mask, int index, double frac)
{
	double y0 = buf[(index + 0) & mask];
	double y1 = buf[(index + 1) & mask];
	double y2 = buf[(index + 2) & mask];
	double y3 = buf[(index + 3) & mask];
	double y4 = buf[(index + 4) & mask];
	double y5 = buf[(index + 5) & mask];

	// Optimal 32x (6-point, 5th-order) (z-form)
	double z = frac - 0.5;
	double even1 = y3 + y2, odd1 = y3 - y2;
	double even2 = y4 + y1, odd2 = y4 - y1;
	double even3 = y5 + y0, odd3 = y5 - y0;

	double c0 = even1 * 0.42685983409379380 + even2 * 0.07238123511170030 + even3 * 0.00075893079450573;
	double c1 = odd1 * 0.35831772348893259 + odd2 * 0.20451644554758297 + odd3 * 0.00562658797241955;
	double c2 = even1 * -0.217009177221292431 + even2 * 0.20051376594086157 + even3 * 0.01649541128040211;
	double c3 = odd1 * -0.25112715343740988 + odd2 * 0.04223025992200458 + odd3 * 0.02488727472995134;
	double c4 = even1 * 0.04166946673533273 + even2 * -0.06250420114356986 + even3 * 0.02083473440841799;
	double c5 = odd1 * 0.08349799235675044 + odd2 * -0.04174912841630993 + odd3 * 0.00834987866042734;
	return ((((c5 * z + c4) * z + c3) * z + c2) * z + c1) * z + c0;
}





