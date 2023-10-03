#include <stdlib.h>
#include <assert.h>

#include "spasm.h"
#include "test_tools.h"

void check_all(i64 prime)
{
	spasm_field F;
	spasm_field_init(prime, F);  
	for (i64 i = 1; i < prime; i++) {
		spasm_ZZp x = spasm_ZZp_init(F, i);
		spasm_ZZp y = spasm_ZZp_inverse(F, x);
		spasm_ZZp k = spasm_ZZp_mul(F, x, y);
		// fprintf(stderr, "%d * %d mod == %d mod %" PRId64 "\n", x, y, k, prime);
		assert(k == 1);
		assert(y <= prime / 2);
		assert(y >= -prime / 2);
	}
	printf("ok inversion mod %" PRId64"\n", prime);	
}

void check_some(i64 prime)
{
	spasm_field F;
	spasm_field_init(prime, F);
	spasm_prng_ctx ctx;  
	spasm_prng_seed_simple(prime, 0, 0, &ctx);
	for (i64 k = 1; k < 10000; k++) {
		spasm_ZZp x = spasm_prng_ZZp(&ctx);
		assert(x <= prime / 2);
		assert(x >= -prime / 2);

		spasm_ZZp y = spasm_ZZp_inverse(F, x);
		assert(y <= prime / 2);
		assert(y >= -prime / 2);

		spasm_ZZp k = spasm_ZZp_mul(F, x, y);
		assert(k <= prime / 2);
		assert(k >= -prime / 2);

		// fprintf(stderr, "%d * %d mod == %d mod %" PRId64 "\n", x, y, k, prime);
		assert(k == 1);
	}
	printf("ok inversion mod %" PRId64"\n", prime);	
	for (i64 k = 1; k < 10000; k++) {
		spasm_ZZp x = spasm_prng_ZZp(&ctx);
		spasm_ZZp y = spasm_prng_ZZp(&ctx);
		spasm_ZZp z = spasm_prng_ZZp(&ctx);
		spasm_ZZp zz = spasm_ZZp_axpy(F, x, y, z);          // zz == x*y + z
		spasm_ZZp o = spasm_ZZp_axpy(F, -x, y, zz);          // o == zz - x*y == z
		assert(o == z);
	}
	printf("ok axpy mod %" PRId64"\n", prime);	
}

int main()
{
	check_all(2);
	check_all(3);
	check_all(257);
	check_all(65537);

	check_some(67108859);             /* largest 26-bit prime */
	check_some(189812507);            /* largest prime that works for ModularBalanced<double> */

	check_some(0x7fffffff);            /* largest 31-bit prime */
	check_some(3037000493);            /* largest prime s.t. a*x+y fits in 63 bits */
	check_some(0xfffffffb);            /* largest 32-bit prime */

	exit(EXIT_SUCCESS);
}