#include "SRY-C.h"

int usage_zhangya()
{
	fprintf(stdout,
        "usage: SRY-C zhangya [options]\n"
	" -k <int>    kmer size, 1 <= k <= 32, [21]\n"
	" -r <string> input kmer file, [NULL]\n"
	" -i <string> input read file, [NULL]\n"
	" -o <string> output file, [NULL]\n"
        "\n"
        "for example:\n"
        "SRY-C zhangya -k 21 -r kmer.txt -i read.fa -o output.txt\n"
        );
        return 1;
}

int usage_zy()
{
        fprintf(stdout,
        "usage: SRY-C zy [options]\n"
        " -k <int>    kmer size, 33 <= k <= 64, [53]\n"
        " -r <string> input kmer file, [NULL]\n"
        " -i <string> input read file, [NULL]\n"
        " -o <string> output file, [NULL]\n"
        "\n"
        "for example:\n"
        "SRY-C zy -k 53 -r kmer.txt -i read.fa -o output.txt\n"
        );
        return 1;
}

int usage_syn()
{
        fprintf(stdout,
        "usage: SRY-C syn [options]\n"
        " -k <int>    kmer size, 65 <= k <= 96, [68]\n"
        " -r <string> input kmer file, [NULL]\n"
        " -i <string> input read file, [NULL]\n"
        " -o <string> output file, [NULL]\n"
        "\n"
        "for example:\n"
        "SRY-C syn -k 68 -r kmer.txt -i read.fa -o output.txt\n"
        );
        return 1;
}


int usage_songyanni()
{
        fprintf(stdout,
        "usage: SRY-C songyanni [options]\n"
        " -k <int>    kmer size, 97 <= k <= 128, [100]\n"
        " -r <string> input kmer file, [NULL]\n"
        " -i <string> input read file, [NULL]\n"
        " -o <string> output file, [NULL]\n"
        "\n"
        "for example:\n"
        "SRY-C songyanni -k 100 -r kmer.txt -i read.fa -o output.txt\n"
        );
        return 1;
}

int main_zhangya(int argc, char **argv)
{
        FileReader *fr;
        BioSequence *seq;
        FILE *kmerfp, *oputfp;
        char *kmerfn, *oputfn;
        char **readfn;
        b4i c, ksize, bi, bj, bk, bl;
        u1i nucleobase[256], nt[4];
        u8i kmer, mask, ui, uj;
        struct KMER_T  kmer_t;
        static struct KMER_V *kmer_v[67108864];
	//
        kmerfn = oputfn = NULL;
        readfn = (char **)malloc(sizeof(char *));
        readfn[0] = NULL;
        ksize = 21;
        while((c = getopt(argc, argv, "hk:r:i:o:")) != -1)
        {
                switch(c)
                {
                        case 'k':
                                ksize = atoi(optarg);
                                break;
                        case 'r':
                                kmerfn = optarg;
                                break;
                        case 'i':
                                readfn[0] = optarg;
                                break;
                        case 'o':
                                oputfn = optarg;
                                break;
                        default:
                                return usage_zhangya();
                }
        }
	if(!(ksize >= 1 && ksize <= 32))
        {
                return usage_zhangya();
        }
        for(bi = 0; bi < 256; bi++)
        {
                nucleobase[bi] = 4;
        }
        nucleobase[65] = nucleobase[97] = 0;
        nucleobase[67] = nucleobase[99] = 1;
        nucleobase[71] = nucleobase[103] = 2;
        nucleobase[84] = nucleobase[85] = nucleobase[116] = nucleobase[117] = 3;
        nt[0] = 65;
        nt[1] = 67;
        nt[2] = 71;
        nt[3] = 84;
        mask = (1ULL << (ksize * 2)) - 1;
        for(ui = 0; ui < 67108864; ui++)
        {
                kmer_v[ui] = (struct KMER_V *)malloc(sizeof(struct KMER_V));
                kmer_v[ui]->size = kmer_v[ui]->cap = 0;
        }
	//
        if((optind == argc && optind == 1) || (optind != argc))
        {
                return usage_zhangya();
        }
	//
        kmerfp = fopen(kmerfn, "r");
        while((c = fgetc(kmerfp)) != -1)
        {
                kmer = nucleobase[c];
                for(bi = 1; bi < ksize; bi++)
                {
                        c = fgetc(kmerfp);
                        kmer = (kmer << 2) | nucleobase[c];
                }
                kmer_t.s = kmer;
                kmer_push(kmer_v[kmer & 67108863], kmer_t);
                c = fgetc(kmerfp);
        }
        fclose(kmerfp);
	//
        fr = open_all_filereader(1, readfn, 0);
        seq = init_biosequence();
        oputfp = fopen(oputfn, "w");
        while(readseq_filereader(fr, seq))
        {
                fprintf(oputfp, "%s,%d", seq->tag->string, seq->seq->size);
                kmer = 0;
                bk = 0;
                for(bj = 0; bj < seq->seq->size; bj++)
                {
                        c = nucleobase[(int)seq->seq->string[bj]];
                        if(c < 4)
                        {
                                kmer = ((kmer << 2) & mask) | c;
                                bk++;
                                if(bk >= ksize)
                                {
                                        uj = kmer & 67108863;
                                        for(ui = 0; ui < kmer_v[uj]->size; ui++)
                                        {
                                                if(kmer_v[uj]->buffer[ui].s == kmer)
                                                {
                                                        fprintf(oputfp, " ");
                                                        for(bl = (ksize - 1) * 2; bl >= 0; bl = bl - 2)
                                                        {
                                                                fprintf(oputfp, "%c", nt[(kmer >> bl) & 3ULL]);
                                                        }
                                                        fprintf(oputfp, ",%d", bj - ksize + 1);
                                                        break;
                                                }
                                        }
                                }
                        }
                        else
                        {
                                bk = 0;
                        }
                }
                fprintf(oputfp, "\n");
        }
        free_biosequence(seq);
        close_filereader(fr);
        fclose(oputfp);
	//
        free(readfn);
        for(ui = 0; ui < 67108864; ui++)
        {
                free(kmer_v[ui]->buffer);
                free(kmer_v[ui]);
        }
	//
        return 0;
}

int main_zy(int argc, char **argv)
{
        FileReader *fr;
        BioSequence *seq;
        FILE *kmerfp, *oputfp;
        char *kmerfn, *oputfn;
        char **readfn;
        b4i c, cc, ksize, bi, bj, bk, bl;
        u1i nucleobase[256], nt[4];
        u8i kmer[4], mask[2], ui, uj;
        struct KMER_T  kmer_t;
        static struct KMER_V *kmer_v[67108864];
	//
        kmerfn = oputfn = NULL;
        readfn = (char **)malloc(sizeof(char *));
        readfn[0] = NULL;
        ksize = 53;
        while((c = getopt(argc, argv, "hk:r:i:o:")) != -1)
        {
                switch(c)
                {
                        case 'k':
                                ksize = atoi(optarg);
                                break;
                        case 'r':
                                kmerfn = optarg;
                                break;
                        case 'i':
                                readfn[0] = optarg;
                                break;
                        case 'o':
                                oputfn = optarg;
                                break;
                        default:
                                return usage_zy();
                }
        }
        if(!(ksize >= 33 && ksize <= 64))
        {
                return usage_zy();
        }
        for(bi = 0; bi < 256; bi++)
        {
                nucleobase[bi] = 4;
        }
        nucleobase[65] = nucleobase[97] = 0;
        nucleobase[67] = nucleobase[99] = 1;
        nucleobase[71] = nucleobase[103] = 2;
        nucleobase[84] = nucleobase[85] = nucleobase[116] = nucleobase[117] = 3;
        nt[0] = 65;
        nt[1] = 67;
        nt[2] = 71;
        nt[3] = 84;
        mask[0] = (1ULL << ((ksize - 32) * 2)) - 1;
        mask[1] = 3ULL << 62;
        for(ui = 0; ui < 67108864; ui++)
        {
                kmer_v[ui] = (struct KMER_V *)malloc(sizeof(struct KMER_V));
                kmer_v[ui]->size = kmer_v[ui]->cap = 0;
        }
	//
        if((optind == argc && optind == 1) || (optind != argc))
        {
                return usage_zy();
        }
	//
        kmerfp = fopen(kmerfn, "r");
        while((c = fgetc(kmerfp)) != -1)
        {
                kmer[0] = nucleobase[c];
                for(bi = 1; bi < ksize - 32; bi++)
                {
                        c = fgetc(kmerfp);
                        kmer[0] = (kmer[0] << 2) | nucleobase[c];
                }
                c = fgetc(kmerfp);
                kmer[1] = nucleobase[c];
                for(bi = 1; bi < 32; bi++)
                {
                        c = fgetc(kmerfp);
                        kmer[1] = (kmer[1] << 2) | nucleobase[c];
                }
                kmer_t.s = kmer[0];
                kmer_push(kmer_v[kmer[1] & 67108863], kmer_t);
                kmer_t.s = kmer[1];
                kmer_push(kmer_v[kmer[1] & 67108863], kmer_t);
                c = fgetc(kmerfp);
        }
        fclose(kmerfp);
	//
        fr = open_all_filereader(1, readfn, 0);
        seq = init_biosequence();
        oputfp = fopen(oputfn, "w");
        while(readseq_filereader(fr, seq))
        {
                fprintf(oputfp, "%s,%d", seq->tag->string, seq->seq->size);
                kmer[0] = kmer[1] = 0;
                bk = 0;
                for(bj = 0; bj < seq->seq->size; bj++)
                {
                        c = nucleobase[(int)seq->seq->string[bj]];
                        if(c < 4)
                        {
                                cc = (kmer[1] & mask[1]) >> 62;
                                kmer[1] = (kmer[1] << 2) | c;
                                c = cc;
                                kmer[0] = ((kmer[0] << 2) & mask[0]) | c;
                                bk++;
                                if(bk >= ksize)
                                {
                                        uj = kmer[1] & 67108863;
                                        for(ui = 0; ui < kmer_v[uj]->size; ui = ui + 2)
                                        {
						if(kmer_v[uj]->buffer[ui + 0].s == kmer[0] && kmer_v[uj]->buffer[ui + 1].s == kmer[1])
                                                {
                                                        fprintf(oputfp, " ");
                                                        for(bl = (ksize - 33) * 2; bl >= 0; bl = bl - 2)
                                                        {
                                                                fprintf(oputfp, "%c", nt[(kmer[0] >> bl) & 3ULL]);
                                                        }
                                                        for(bl = (32 - 1) * 2; bl >= 0; bl = bl - 2)
                                                        {
                                                                fprintf(oputfp, "%c", nt[(kmer[1] >> bl) & 3ULL]);
                                                        }
                                                        fprintf(oputfp, ",%d", bj - ksize + 1);
                                                        break;
                                                }
                                        }
                                }
                        }
                        else
                        {
                                bk = 0;
                        }
                }
                fprintf(oputfp, "\n");
        }
        free_biosequence(seq);
        close_filereader(fr);
        fclose(oputfp);
	//
        free(readfn);
        for(ui = 0; ui < 67108864; ui++)
        {
                free(kmer_v[ui]->buffer);
                free(kmer_v[ui]);
        }
	//
        return 0;
}

int main_syn(int argc, char **argv)
{
        FileReader *fr;
        BioSequence *seq;
        FILE *kmerfp, *oputfp;
        char *kmerfn, *oputfn;
        char **readfn;
        b4i c, cc, ksize, bi, bj, bk, bl;
        u1i nucleobase[256], nt[4];
        u8i kmer[4], mask[2], ui, uj;
        struct KMER_T  kmer_t;
        static struct KMER_V *kmer_v[67108864];
	//
        kmerfn = oputfn = NULL;
        readfn = (char **)malloc(sizeof(char *));
        readfn[0] = NULL;
        ksize = 68;
        while((c = getopt(argc, argv, "hk:r:i:o:")) != -1)
        {
                switch(c)
                {
                        case 'k':
                                ksize = atoi(optarg);
                                break;
                        case 'r':
                                kmerfn = optarg;
                                break;
                        case 'i':
                                readfn[0] = optarg;
                                break;
                        case 'o':
                                oputfn = optarg;
                                break;
                        default:
                                return usage_syn();
                }
        }
        if(!(ksize >= 65 && ksize <= 96))
        {
                return usage_syn();
        }
        for(bi = 0; bi < 256; bi++)
        {
                nucleobase[bi] = 4;
        }
        nucleobase[65] = nucleobase[97] = 0;
        nucleobase[67] = nucleobase[99] = 1;
        nucleobase[71] = nucleobase[103] = 2;
        nucleobase[84] = nucleobase[85] = nucleobase[116] = nucleobase[117] = 3;
        nt[0] = 65;
        nt[1] = 67;
        nt[2] = 71;
        nt[3] = 84;
        mask[0] = (1ULL << ((ksize - 64) * 2)) - 1;
        mask[1] = 3ULL << 62;
        for(ui = 0; ui < 67108864; ui++)
        {
                kmer_v[ui] = (struct KMER_V *)malloc(sizeof(struct KMER_V));
                kmer_v[ui]->size = kmer_v[ui]->cap = 0;
        }
	//
        if((optind == argc && optind == 1) || (optind != argc))
        {
                return usage_syn();
        }
	//
        kmerfp = fopen(kmerfn, "r");
        while((c = fgetc(kmerfp)) != -1)
        {
                kmer[0] = nucleobase[c];
                for(bi = 1; bi < ksize - 64; bi++)
                {
                        c = fgetc(kmerfp);
                        kmer[0] = (kmer[0] << 2) | nucleobase[c];
                }
                c = fgetc(kmerfp);
                kmer[1] = nucleobase[c];
                for(bi = 1; bi < 32; bi++)
                {
                        c = fgetc(kmerfp);
                        kmer[1] = (kmer[1] << 2) | nucleobase[c];
                }
                c = fgetc(kmerfp);
                kmer[2] = nucleobase[c];
                for(bi = 1; bi < 32; bi++)
                {
                        c = fgetc(kmerfp);
                        kmer[2] = (kmer[2] << 2) | nucleobase[c];
                }
                kmer_t.s = kmer[0];
                kmer_push(kmer_v[kmer[2] & 67108863], kmer_t);
                kmer_t.s = kmer[1];
                kmer_push(kmer_v[kmer[2] & 67108863], kmer_t);
                kmer_t.s = kmer[2];
                kmer_push(kmer_v[kmer[2] & 67108863], kmer_t);
                c = fgetc(kmerfp);
        }
        fclose(kmerfp);
	//
        fr = open_all_filereader(1, readfn, 0);
        seq = init_biosequence();
        oputfp = fopen(oputfn, "w");
        while(readseq_filereader(fr, seq))
        {
                fprintf(oputfp, "%s,%d", seq->tag->string, seq->seq->size);
                kmer[0] = kmer[1] = kmer[2] = 0;
                bk = 0;
                for(bj = 0; bj < seq->seq->size; bj++)
                {
                        c = nucleobase[(int)seq->seq->string[bj]];
                        if(c < 4)
                        {
                                cc = (kmer[2] & mask[1]) >> 62;
                                kmer[2] = (kmer[2] << 2) | c;
                                c = cc;
                                cc = (kmer[1] & mask[1]) >> 62;
                                kmer[1] = (kmer[1] << 2) | c;
                                c = cc;
                                kmer[0] = ((kmer[0] << 2) & mask[0]) | c;
                                bk++;
                                if(bk >= ksize)
                                {
                                        uj = kmer[2] & 67108863;
                                        for(ui = 0; ui < kmer_v[uj]->size; ui = ui + 3)
                                        {
						if(kmer_v[uj]->buffer[ui + 0].s == kmer[0] && kmer_v[uj]->buffer[ui + 1].s == kmer[1] && kmer_v[uj]->buffer[ui + 2].s == kmer[2])
                                                {
                                                        fprintf(oputfp, " ");
                                                        for(bl = (ksize - 65) * 2; bl >= 0; bl = bl - 2)
                                                        {
                                                                fprintf(oputfp, "%c", nt[(kmer[0] >> bl) & 3ULL]);
                                                        }
                                                        for(bl = (32 - 1) * 2; bl >= 0; bl = bl - 2)
                                                        {
                                                                fprintf(oputfp, "%c", nt[(kmer[1] >> bl) & 3ULL]);
                                                        }
                                                        for(bl = (32 - 1) * 2; bl >= 0; bl = bl - 2)
                                                        {
                                                                fprintf(oputfp, "%c", nt[(kmer[2] >> bl) & 3ULL]);
                                                        }
                                                        fprintf(oputfp, ",%d", bj - ksize + 1);
                                                        break;
                                                }
                                        }
                                }
                        }
                        else
                        {
                                bk = 0;
                        }
                }
                fprintf(oputfp, "\n");
        }
        free_biosequence(seq);
        close_filereader(fr);
        fclose(oputfp);
	//
        free(readfn);
        for(ui = 0; ui < 67108864; ui++)
        {
                free(kmer_v[ui]->buffer);
                free(kmer_v[ui]);
        }
	//
        return 0;
}

int main_songyanni(int argc, char **argv)
{
	FileReader *fr;
        BioSequence *seq;
	FILE *kmerfp, *oputfp;
	char *kmerfn, *oputfn;
	char **readfn;
	b4i c, cc, ksize, bi, bj, bk, bl;
	u1i nucleobase[256], nt[4];
        u8i kmer[4], mask[2], ui, uj;
	struct KMER_T  kmer_t;
	static struct KMER_V *kmer_v[67108864];
	//
	kmerfn = oputfn = NULL;
	readfn = (char **)malloc(sizeof(char *));
	readfn[0] = NULL;
	ksize = 100;
        while((c = getopt(argc, argv, "hk:r:i:o:")) != -1)
        {
                switch(c)
                {
			case 'k':
                                ksize = atoi(optarg);
                                break;
			case 'r':
                                kmerfn = optarg;
                                break;
			case 'i':
                                readfn[0] = optarg;
                                break;
			case 'o':
                                oputfn = optarg;
                                break;
                        default:
                                return usage_songyanni();
                }
        }
	if(!(ksize >= 97 && ksize <= 128))
	{
		return usage_songyanni();
	}
	for(bi = 0; bi < 256; bi++)
	{
		nucleobase[bi] = 4;
	}
	nucleobase[65] = nucleobase[97] = 0;
	nucleobase[67] = nucleobase[99] = 1;
        nucleobase[71] = nucleobase[103] = 2;
        nucleobase[84] = nucleobase[85] = nucleobase[116] = nucleobase[117] = 3;
	nt[0] = 65;
	nt[1] = 67;
	nt[2] = 71;
	nt[3] = 84;
	mask[0] = (1ULL << ((ksize - 96) * 2)) - 1;
	mask[1] = 3ULL << 62;
	for(ui = 0; ui < 67108864; ui++)
	{
		kmer_v[ui] = (struct KMER_V *)malloc(sizeof(struct KMER_V));
		kmer_v[ui]->size = kmer_v[ui]->cap = 0;
	}
	//
	if((optind == argc && optind == 1) || (optind != argc))
        {
		return usage_songyanni();
        }
	//
	kmerfp = fopen(kmerfn, "r");
	while((c = fgetc(kmerfp)) != -1)
        {
		kmer[0] = nucleobase[c];
		for(bi = 1; bi < ksize - 96; bi++)
		{
			c = fgetc(kmerfp);
			kmer[0] = (kmer[0] << 2) | nucleobase[c];
		}
		c = fgetc(kmerfp);
		kmer[1] = nucleobase[c];
		for(bi = 1; bi < 32; bi++)
                {
                        c = fgetc(kmerfp);
                        kmer[1] = (kmer[1] << 2) | nucleobase[c];
                }
		c = fgetc(kmerfp);
                kmer[2] = nucleobase[c];
                for(bi = 1; bi < 32; bi++)
                {
                        c = fgetc(kmerfp);
                        kmer[2] = (kmer[2] << 2) | nucleobase[c];
                }
		c = fgetc(kmerfp);
                kmer[3] = nucleobase[c];
                for(bi = 1; bi < 32; bi++)
                {
                        c = fgetc(kmerfp);
                        kmer[3] = (kmer[3] << 2) | nucleobase[c];
                }
		kmer_t.s = kmer[0];
		kmer_push(kmer_v[kmer[3] & 67108863], kmer_t);
		kmer_t.s = kmer[1];
                kmer_push(kmer_v[kmer[3] & 67108863], kmer_t);
		kmer_t.s = kmer[2];
                kmer_push(kmer_v[kmer[3] & 67108863], kmer_t);
		kmer_t.s = kmer[3];
                kmer_push(kmer_v[kmer[3] & 67108863], kmer_t);
		c = fgetc(kmerfp);
        }
	fclose(kmerfp);
	//
	fr = open_all_filereader(1, readfn, 0);
	seq = init_biosequence();
	oputfp = fopen(oputfn, "w");
        while(readseq_filereader(fr, seq))
	{
		fprintf(oputfp, "%s,%d", seq->tag->string, seq->seq->size);
		kmer[0] = kmer[1] = kmer[2] = kmer[3] = 0;
		bk = 0;
		for(bj = 0; bj < seq->seq->size; bj++)
		{
			c = nucleobase[(int)seq->seq->string[bj]];
			if(c < 4)
			{
				cc = (kmer[3] & mask[1]) >> 62;
				kmer[3] = (kmer[3] << 2) | c;
				c = cc;
				cc = (kmer[2] & mask[1]) >> 62;
				kmer[2] = (kmer[2] << 2) | c;
				c = cc;
				cc = (kmer[1] & mask[1]) >> 62;
				kmer[1] = (kmer[1] << 2) | c;
				c = cc;
				kmer[0] = ((kmer[0] << 2) & mask[0]) | c;
				bk++;
				if(bk >= ksize)
				{
					uj = kmer[3] & 67108863;
					for(ui = 0; ui < kmer_v[uj]->size; ui = ui + 4)
					{
						if(kmer_v[uj]->buffer[ui + 0].s == kmer[0] && kmer_v[uj]->buffer[ui + 1].s == kmer[1] && kmer_v[uj]->buffer[ui + 2].s == kmer[2] && kmer_v[uj]->buffer[ui + 3].s == kmer[3])
						{
							fprintf(oputfp, " ");
							for(bl = (ksize - 97) * 2; bl >= 0; bl = bl - 2)
							{
								fprintf(oputfp, "%c", nt[(kmer[0] >> bl) & 3ULL]);
							}
							for(bl = (32 - 1) * 2; bl >= 0; bl = bl - 2)
                                                        {
                                                                fprintf(oputfp, "%c", nt[(kmer[1] >> bl) & 3ULL]);
                                                        }
							for(bl = (32 - 1) * 2; bl >= 0; bl = bl - 2)
                                                        {
                                                                fprintf(oputfp, "%c", nt[(kmer[2] >> bl) & 3ULL]);
                                                        }
							for(bl = (32 - 1) * 2; bl >= 0; bl = bl - 2)
                                                        {
                                                                fprintf(oputfp, "%c", nt[(kmer[3] >> bl) & 3ULL]);
                                                        }
							fprintf(oputfp, ",%d", bj - ksize + 1);
							break;
						}
					}
				}
			}
			else
			{
				bk = 0;
			}
		}
		fprintf(oputfp, "\n");
	}
	free_biosequence(seq);
        close_filereader(fr);
	fclose(oputfp);
	//
	free(readfn);
	for(ui = 0; ui < 67108864; ui++)
	{
		free(kmer_v[ui]->buffer);
        	free(kmer_v[ui]);
	}
	//
        return 0;
}

int usage()
{
        fprintf(stdout,
	"Program: SRY-C\n"
	"Version: %s\n"
	"Author : Shigang Wu <wushigang@caas.cn>\n"
	"Usage  : SRY-C <cmd> [options]\n"
	"\n"
	"commands:\n"
        " zhangya       kmer size in  1 to  32\n"
        " zy            kmer size in 33 to  64\n"
        " syn           kmer size in 65 to  96\n"
        " songyanni     kmer size in 97 to 128\n"
	"example:\n"
	"# To run kmer size in  1 to  32\n"
	"  SRY-C zhangya\n"
	"# To run kmer size in 33 to  64\n"
        "  SRY-C zy\n"
	"# To run kmer size in 65 to  96\n"
        "  SRY-C syn\n"
	"# To run kmer size in 97 to 128\n"
        "  SRY-C songyanni\n"
	"\n", TOSTR(VERSION)
        );
        return 1;
}

int main(int argc, char **argv)
{
	if(argc < 2)
        {
                return usage();
        }
	if(strcasecmp("zhangya", argv[1]) == 0) return main_zhangya(argc - 1, argv + 1);
        if(strcasecmp("zy", argv[1]) == 0) return main_zy(argc - 1, argv + 1);
        if(strcasecmp("syn", argv[1]) == 0) return main_syn(argc - 1, argv + 1);
        if(strcasecmp("songyanni", argv[1]) == 0) return main_songyanni(argc - 1, argv + 1);
	if(strcasecmp("-h", argv[1]) == 0) return usage();
        if(strcasecmp("--help", argv[1]) == 0) return usage();
	fprintf(stderr, " -- unknown command '%s' -- \n", argv[1]);
	return 1;
}
