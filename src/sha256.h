#include <stddef.h>
#include "spasm.h"

typedef struct {
    u32 h[8];
    u32 Nl, Nh;
    u32 data[16];
    u32 num, md_len;
} SHA256_CTX;

void spasm_SHA256_init(SHA256_CTX *c);
void spasm_SHA256_update(SHA256_CTX *c, const void *data, size_t len);
void spasm_SHA256_final(u8 *md, SHA256_CTX *c);
