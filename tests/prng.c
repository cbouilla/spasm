#include <stdlib.h>
#include <assert.h>

#include "spasm.h"
#include "test_tools.h"

void check(i64 prime, u64 seed, u32 sequence)
{
	spasm_prng_ctx ctx;
	spasm_prng_seed_simple(prime, seed, sequence, &ctx);
	printf("prime=%" PRId64 ", seed=%016" PRIx64 ", seq=%08x, out=", prime, seed, sequence);
	for (int i = 0; i < 10; i++) {
		if (i > 0)
			printf(", ");
		printf("%6d", spasm_prng_ZZp(&ctx));
	}
	printf("\n");
}

int main()
{
	check(257, 0, 0);
	check(257, 0, 1);
	check(257, 1, 0);
	check(257, 1, 1);

	check(65537, 0xdead00000000beefull, 0);
	
	exit(EXIT_SUCCESS);
}