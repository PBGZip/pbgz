#include "../actg.h"

void test_actg_stretch_mapping_xor()
{
    int32_t n;
    uint8_t base[32], refe[32];
    n = 1; base[n - 1] = 'C';n = 2; base[n - 1] = 'A';
    n = 3; base[n - 1] = 'G';n = 4; base[n - 1] = 'G';
    n = 5; base[n - 1] = 'G';n = 6; base[n - 1] = 'A';
    n = 7; base[n - 1] = 'C';n = 8; base[n - 1] = 'T';

    n = 9; base[n - 1] = 'G';n = 10; base[n - 1] = 'G';
    n = 11; base[n - 1] = 'A';n = 12; base[n - 1] = 'T';
    n = 13; base[n - 1] = 'G';n = 14; base[n - 1] = 'G';
    n = 15; base[n - 1] = 'G';n = 16; base[n - 1] = 'G';

    n = 17; base[n - 1] = 'G';n = 18; base[n - 1] = 'T';
    n = 19; base[n - 1] = 'A';n = 20; base[n - 1] = 'T';
    n = 21; base[n - 1] = 'G';n = 22; base[n - 1] = 'C';
    n = 23; base[n - 1] = 'T';n = 24; base[n - 1] = 'A';

    n = 25; base[n - 1] = 'G';n = 26; base[n - 1] = 'G';
    n = 27; base[n - 1] = 'G';n = 28; base[n - 1] = 'A';
    n = 29; base[n - 1] = 'C';n = 30; base[n - 1] = 'T';
    n = 31; base[n - 1] = 'T';n = 32; base[n - 1] = 'G';
    memcpy(refe, base, 32);
    n = 1; refe[n - 1] = 'T'; 
    n = 2; refe[n - 1] = 'T';
    n = 8; refe[n - 1] = 'A';
    n = 12; refe[n - 1] = 'A';
    n = 16; refe[n - 1] = 'A';
    n = 20; refe[n - 1] = 'A';
    n = 24; refe[n - 1] = 'T';
    n = 28; refe[n - 1] = 'T';
    n = 32; refe[n - 1] = 'T';

    uint8_t squash_base[8], squash_refe[8];
    uint8_t out[32];
    actg_squash(base, 32, squash_base);
    actg_squash(refe, 32, squash_refe);

    actg_stretch_mapping_xor(squash_base, squash_refe, 8, out);

    printf("out:\n");
    for (n = 0; n < 32; n++) {
        printf("%3u,", out[n]);
        if ((n + 1) % 16 == 0)
            printf("\n");
    }
}

int main()
{

    test_actg_stretch_mapping_xor();
    return 0;
}