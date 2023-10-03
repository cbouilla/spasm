#include <string.h>
#include <stdlib.h>

#include "spasm.h"

void check(const char *msg)
{
	spasm_sha256_ctx ctx;
	u8 hash[32];
	spasm_SHA256_init(&ctx);
	spasm_SHA256_update(&ctx, msg, strlen(msg));
	spasm_SHA256_final(hash, &ctx);
	for (int i = 0; i < 32; i++)
		printf("%02x", hash[i]);
	printf("\n");
}

int main()
{
	check("");
	check("X");
	check("Hello World");
	check("abcdefghijklmnopqrstuvwxyz0123456789ABCDEFGHIJKLMNOPQRSTUVWXYZ+-*/=");
		
	exit(EXIT_SUCCESS);
}