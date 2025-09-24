// These constants from paper:

// "Polynomial Interpolators for
// High-Quality Resampling of
// Oversampled Audio"

// Olli Niemitalo in October 2001.

#if 0
inline float sample_2x_6p_5z(float* buf, int mask, int index, float frac)
{
	float y0 = buf[(index + 0) & mask];
	float y1 = buf[(index + 1) & mask];
	float y2 = buf[(index + 2) & mask];
	float y3 = buf[(index + 3) & mask];
	float y4 = buf[(index + 4) & mask];
	float y5 = buf[(index + 5) & mask];

	// Optimal 2x (6-point, 5th-order) (z-form)
	float z = frac - 0.5f;
	float even1 = y3 + y2, odd1 = y3 - y2;
	float even2 = y4 + y1, odd2 = y4 - y1;
	float even3 = y5 + y0, odd3 = y5 - y0;

	float c0 = even1 * 0.40513396007145713f + even2 * 0.09251794438424393f		+ even3 * 0.00234806603570670f;
	float c1 = odd1 * 0.28342806338906690f + odd2 * 0.21703277024054901f		+ odd3 * 0.01309294748731515f;
	float c2 = even1 * -0.191337682540351941f + even2 * 0.16187844487943592f		+ even3 * 0.02946017143111912f;
	float c3 = odd1 * -0.16471626190554542f + odd2 * -0.00154547203542499f		+ odd3 * 0.03399271444851909f;
	float c4 = even1 * 0.03845798729588149f + even2 * -0.05712936104242644f		+ even3 * 0.01866750929921070f;
	float c5 = odd1 * 0.04317950185225609f + odd2 * -0.01802814255926417f		+ odd3 * 0.00152170021558204f;
	return ((((c5 * z + c4) * z + c3) * z + c2) * z + c1) * z + c0;
}
#endif

#if 0
inline float sample_4x_6p_5z(float* buf, int mask, int index, float frac)
{
	float y0 = buf[(index + 0) & mask];
	float y1 = buf[(index + 1) & mask];
	float y2 = buf[(index + 2) & mask];
	float y3 = buf[(index + 3) & mask];
	float y4 = buf[(index + 4) & mask];
	float y5 = buf[(index + 5) & mask];

	// Optimal 4x (6-point, 5th-order) (z-form)
	float z = frac - 0.5f;
	float even1 = y3 + y2, odd1 = y3 - y2;
	float even2 = y4 + y1, odd2 = y4 - y1;
	float even3 = y5 + y0, odd3 = y5 - y0;

	float c0 = even1 * 0.414969029592408940f + even2 * 0.08343081932889224f + even3 * 0.00160015038681571f;
	float c1 = odd1 * 0.316255150048597830f + odd2 * 0.21197848565176958f + odd3 * 0.00956166668408054f;
	float c2 = even1 * -0.203271896548875371f + even2 * 0.17989908432249280f + even3 * 0.02337283412161328f;
	float c3 = odd1 * -0.202092410698357320f + odd2 * 0.01760734419526000f + odd3 * 0.02985927012435252f;
	float c4 = even1 * 0.041009488587619100f + even2 * -0.06147760875085254f + even3 * 0.02046802954581191f;
	float c5 = odd1 * 0.066077478644169240f + odd2 * -0.03255079211953620f + odd3 * 0.00628989632244913f;

	return ((((c5 * z + c4) * z + c3) * z + c2) * z + c1) * z + c0;
}
#endif

#if 0
inline float sample_8x_6p_5z(float* buf, int mask, int index, float frac)
{
	float y0 = buf[(index + 0) & mask];
	float y1 = buf[(index + 1) & mask];
	float y2 = buf[(index + 2) & mask];
	float y3 = buf[(index + 3) & mask];
	float y4 = buf[(index + 4) & mask];
	float y5 = buf[(index + 5) & mask];

	// Optimal 4x (6-point, 5th-order) (z-form)
	float z = frac - 0.5f;
	float even1 = y3 + y2, odd1 = y3 - y2;
	float even2 = y4 + y1, odd2 = y4 - y1;
	float even3 = y5 + y0, odd3 = y5 - y0;

	float c0 = even1 * 0.41660797292569773f + even2 * 0.08188468587188069f + even3 * 0.00150734119050266f;
	float c1 = odd1 * 0.32232780822726981f + odd2 * 0.21076321997422021f + odd3 * 0.00907649978070957f;
	float c2 = even1 * -0.205219993961471501f + even2 * 0.18282942057327367f + even3 * 0.02239057377093268f;
	float c3 = odd1 * -0.21022298520246224f + odd2 * 0.02176417471349534f + odd3 * 0.02898626924395209f;
	float c4 = even1 * 0.04149963966704384f + even2 * -0.06224707096203808f + even3 * 0.02074742969707599f;
	float c5 = odd1 * 0.07517133281176167f + odd2 * -0.03751837438141215f + odd3 * 0.00747588873055296f;

	return ((((c5 * z + c4) * z + c3) * z + c2) * z + c1) * z + c0;
}
#endif

#if 1
inline float sample_16x_6p_5z(float* buf, int mask, int index, float frac)
{
	float y0 = buf[(index + 0) & mask];
	float y1 = buf[(index + 1) & mask];
	float y2 = buf[(index + 2) & mask];
	float y3 = buf[(index + 3) & mask];
	float y4 = buf[(index + 4) & mask];
	float y5 = buf[(index + 5) & mask];

	// Optimal 16x (6-point, 5th-order) (z-form)
	float z = frac - 0.5f;
	float even1 = y3 + y2, odd1 = y3 - y2;
	float even2 = y4 + y1, odd2 = y4 - y1;
	float even3 = y5 + y0, odd3 = y5 - y0;

	float c0 = even1 * 0.41809989254549901f + even2 * 0.08049339946273310f + even3 * 0.00140670799165932f;
	float c1 = odd1 * 0.32767596257424964f + odd2 * 0.20978189376640677f + odd3 * 0.00859567104974701f;
	float c2 = even1 * -0.206944618112960001f + even2 * 0.18541689550861262f + even3 * 0.02152772260740132f;
	float c3 = odd1 * -0.21686095413034051f + odd2 * 0.02509557922091643f + odd3 * 0.02831484751363800f;
	float c4 = even1 * 0.04163046817137675f + even2 * -0.06244556931623735f + even3 * 0.02081510113314315f;
	float c5 = odd1 * 0.07990500783668089f + odd2 * -0.03994519162531633f + odd3 * 0.00798609327859495f;

	return ((((c5 * z + c4) * z + c3) * z + c2) * z + c1) * z + c0;
}
#endif


