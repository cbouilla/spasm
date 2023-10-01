#include <stdlib.h>
#include <assert.h>

#include "spasm.h"
#include "test_tools.h"

void check(u64 seed, u64 sequence)
{
	spasm_prng_seed(seed, sequence);
	printf("seed=%016" PRIx64 ", seq=%016" PRIx64 ", out=", seed, sequence);
	for (int i = 0; i < 4; i++) {
		if (i > 0)
			printf(", ");
		printf("%016" PRIx64, spasm_prng_next());
	}
	printf("\n");
}

int main()
{
	check(0, 0);
	check(0, 1);
	check(1, 0);
	check(1, 1);

	check(0xff00000000000000ull, 0xdead00000000beefull);
	
	exit(EXIT_SUCCESS);
}