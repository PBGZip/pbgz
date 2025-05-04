#include <stdlib.h>
#include <memory.h>

#include "fc_model.h"
#include "fc_header.h"

fcModel1 gfcModel1;
fcModel2 gfcModel2;

void fcmemset(void * dst, int size, short v)
{
    for (int i = 0; i < size / 2; ++i) ((short *)dst)[i] = v;
}

int fcinit_model()
{
    for (int mixer = 0; mixer < CHAR_SIZE; ++mixer)
    {
        gfcModel1.tid_mixer[mixer].Init();
        gfcModel1.tid_mixerout_of_range[mixer].Init();
        gfcModel1.run_mixer[mixer].Init();
    }
    for (int bit = 0; bit < 8; ++bit)
    {
        gfcModel1.tid_mixerbits_value[bit].Init();
        for (int context = 0; context < 8; ++context)
            gfcModel1.tid_mixerbit_width[context][bit].Init();
    }
    for (int bit = 0; bit < 32; ++bit)
    {
        gfcModel1.run_mixerbits_value[bit].Init();
        for (int context = 0; context < 32; ++context)
            gfcModel1.run_mixerbit_width[context][bit].Init();
    }

    fcmemset(&gfcModel1.tid_t, sizeof(gfcModel1.tid_t), 2048);
    fcmemset(&gfcModel1.run_t, sizeof(gfcModel1.run_t), 2048);

    fcmemset(&gfcModel2.tid_t, sizeof(gfcModel2.tid_t), 4096);
    fcmemset(&gfcModel2.run_t, sizeof(gfcModel2.run_t), 1024);

    return FC_OK;
}

void fcinit_model(fcModel1 * model)
{
    memcpy(model, &gfcModel1, sizeof(fcModel1));
}

void fcinit_model(fcModel2 * model)
{
    memcpy(model, &gfcModel2, sizeof(fcModel2));
}
