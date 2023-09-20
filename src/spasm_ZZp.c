#include <assert.h>

#include "spasm.h"

void spasm_field_init(i64 p, spasm_field F)
{
	F->p = p;
	if (p < 0)
		return;
	assert(2 <= p);
	assert(p <= 0xfffffffbLL);
	F->halfp = p / 2;
	F->mhalfp = p / 2 - p + 1;
	F->dinvp = 1. / ((double) p);
}

static inline spasm_ZZp NORMALISE(const spasm_field F, i64 x)
{
	if (x < F->mhalfp)
		x += F->p;
	else if (x > F->halfp)
		x -= F->p;
	return x;
}

inline spasm_ZZp spasm_ZZp_init(const spasm_field F, i64 x)
{
	i64 p = F->p;
	return NORMALISE(F, x % p);
}

inline spasm_ZZp spasm_ZZp_add(const spasm_field F, spasm_ZZp a, spasm_ZZp b)
{
	return NORMALISE(F, (i64) a + (i64) b);
}

inline spasm_ZZp spasm_ZZp_sub(const spasm_field F, spasm_ZZp a, spasm_ZZp b)
{
	return NORMALISE(F, (i64) a - (i64) b);
}

inline spasm_ZZp spasm_ZZp_mul(const spasm_field F, spasm_ZZp a, spasm_ZZp b)
{
	i64 q = ((double) a) * ((double) b) * F->dinvp;
	return NORMALISE(F, (i64) a * (i64) b - q * F->p);
}

/* compute bezout relation u*a + v*p == 1; returns u */
static i64 gcdext(i64 a, i64 p)
{
	assert(a >= 0);
	i64 t = 0, u = 1;
	i64 r = p, s = a;
	while (s != 0) {
		i64 q = r / s;
		i64 foo = u;
		u = t - q * u;
		t = foo;

		i64 bar = s;
		s = r - q * s;
		r = bar;
	}
	return t;
}

spasm_ZZp spasm_ZZp_inverse(const spasm_field F, spasm_ZZp a)
{
	i64 aa = a;
	if (aa < 0)
		aa += F->p;
	i64 inva = gcdext(aa, F->p);
	return NORMALISE(F, inva);
}


inline spasm_ZZp spasm_ZZp_axpy(const spasm_field F, spasm_ZZp a, spasm_ZZp x, spasm_ZZp y)
{
	i64 q = (((((double) a) * ((double) x)) + (double) y) * F->dinvp);
	i64 aa = a;
	i64 xx = x;
	i64 yy = y;
	return NORMALISE(F, aa * xx + yy - q * F->p);
}
