/*
Supported platforms are currently:

  * Linux
  * Mac OS X
  * Window
*/
#include <stdio.h>
#include "cfgpath.h"

int main()
{

    char cfgdir[MAX_PATH];
    get_user_config_file(cfgdir, sizeof(cfgdir), "myapp");
    if (cfgdir[0] == 0)
    {
        printf("Unable to find home directory.\n");
        return 1;
    }
    printf("Saving config file to %s\n", cfgdir);

    get_user_config_folder(cfgdir, sizeof(cfgdir), "myapp");
    if (cfgdir[0] == 0)
    {
        printf("Unable to find home directory.\n");
        return 1;
    }
    printf("Saving  config folder to %s\n", cfgdir);

    get_user_data_folder(cfgdir, sizeof(cfgdir), "myapp");
    if (cfgdir[0] == 0)
    {
        printf("Unable to find home directory.\n");
        return 1;
    }
    printf("Saving data folder to %s\n", cfgdir);

    get_user_cache_folder(cfgdir, sizeof(cfgdir), "myapp");
    if (cfgdir[0] == 0)
    {
        printf("Unable to find home directory.\n");
        return 1;
    }
    printf("Saving cache folder to %s\n", cfgdir);
}
