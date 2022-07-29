#include "bsalign/filereader.h"

struct KMER_T
{
	u8i s;
};

struct KMER_V
{
        u8i size;
	u8i cap;
	struct KMER_T *buffer;
};

void kmer_push(struct KMER_V *v, struct KMER_T t)
{
	if(v->size == v->cap)
	{
		v->cap = v->cap ? v->cap << 1 : 2;
		v->buffer = (struct KMER_T *)realloc(v->buffer, sizeof(struct KMER_T) * v->cap);
	}
	v->buffer[v->size++] = t;
}
