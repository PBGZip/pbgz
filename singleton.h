#ifndef _SINGLETON_H_
#define _SINGLETON_H_

#include <iostream>
#include <stdlib.h>
#include <pthread.h>
#include <unistd.h>
#include <sys/types.h>
#include <stdio.h>
#include <errno.h>
#include <string.h>

using namespace std;

template <typename T>
class singleton
{

public:
    static T &instance()
    {
        Init_();
        return *instance_;
    }

private:
    static void Init_()
    {
        if (instance_ == 0)
        {

            pthread_mutex_lock(&mutex);
            if (instance_ == 0)
            {
                instance_ = new T;
                atexit(Destroy); 
            }
            pthread_mutex_unlock(&mutex);
        }
    }

    static void Destroy()
    {
        delete instance_;
    }

    singleton(const singleton &other);
    singleton &operator=(const singleton &other);
    singleton();
    ~singleton();

    static T *volatile instance_;
    static pthread_mutex_t mutex;
};

template <typename T>
T *volatile singleton<T>::instance_ = 0;

template <typename T>
pthread_mutex_t singleton<T>::mutex = PTHREAD_MUTEX_INITIALIZER;

#endif