#include <arpa/inet.h>           // htonl

#include "spasm.h"

/* This PRNG is SHA256 in counter mode */

static void rehash(spasm_prng_ctx *ctx)
{
        spasm_sha256_ctx hctx;
        spasm_SHA256_init(&hctx);
        spasm_SHA256_update(&hctx, ctx->block, 44);
        spasm_SHA256_final((u8 *) ctx->hash, &hctx);
        ctx->counter += 1;
        ctx->block[9] = htonl(ctx->counter);
        ctx->i = 0;
}

/*
 * Return a uniformly random 32-bit integer
 */
u32 spasm_prng_u32(spasm_prng_ctx *ctx)
{
        if (ctx->i == 8)
                rehash(ctx);
        u32 res = ctx->hash[ctx->i];
        ctx->i += 1;
        return htonl(res);
}

/*
 * Return a uniformly integer modulo prime (rejection sampling)
 */
spasm_ZZp spasm_prng_ZZp(spasm_prng_ctx *ctx)
{
        for (;;) {
                u32 x = spasm_prng_u32(ctx) & ctx->mask;
                if (x < ctx->prime)
                        return spasm_ZZp_init(ctx->field, x);
        }
}

/*
 * Seed must be 32 bytes.
 */
void spasm_prng_seed(const u8 *seed, i64 prime, u32 seq, spasm_prng_ctx *ctx)
{
        u8 *block8 = (u8 *) ctx->block;
        for (int i = 0; i < 32; i++)
                block8[i] = seed[i];
        ctx->prime = prime;
        i64 mask = 1;
        while (mask < prime)
                mask <<= 1;
        ctx->mask = mask - 1;
        ctx->block[8] = htonl(prime);
        ctx->block[9] = 0;
        ctx->block[10] = htonl(seq);
        ctx->counter = 0;
        spasm_field_init(prime, ctx->field);
        rehash(ctx);
}

/*
 * In case where a 32-byte seed (i.e. a SHA256 digest) is not available
 */
void spasm_prng_seed_simple(i64 prime, u64 seed, u32 seq, spasm_prng_ctx *ctx)
{
        u32 block[8];
        block[0] = htonl(seed & 0xffffffff);
        block[1] = htonl(seed >> 32);
        for (int i = 2; i < 8; i++)
                block[i] = 0;
        spasm_prng_seed((u8 *) block, prime, seq, ctx);
}