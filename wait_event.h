/*
 * wait-event.h
 *
 */

#ifndef WAIT_EVENT_H_
#define WAIT_EVENT_H_

#include <pthread.h>
#include <stdint.h>

//Comment out the #define below to disable WFMO support if not used (recommended)
//Compiling with WFMO support will add some overhead to all event objects
#define WFMO 1

namespace neosmart {
//Type declarations
struct neosmart_event_t_;
typedef neosmart_event_t_ * neosmart_event_t;

//WIN32-style functions
neosmart_event_t CreateEvent(bool manualReset = false,
        bool initialState = false);
int DestroyEvent(neosmart_event_t event);
int WaitForEvent(neosmart_event_t event, uint64_t milliseconds = -1);
int SetEvent(neosmart_event_t event);
int ResetEvent(neosmart_event_t event);
#ifdef WFMO
int WaitForMultipleEvents(neosmart_event_t *events, int count, bool waitAll,
        uint64_t milliseconds);
int WaitForMultipleEvents(neosmart_event_t *events, int count, bool waitAll,
        uint64_t milliseconds, int &index);
#endif

//POSIX-style functions
//TBD
}

class wait_event {
    public:
        neosmart::neosmart_event_t evt; // simulate Windows Event Object

        wait_event();
        void reset();
        void wait();
        void wakeup();

        ~wait_event();
};

#endif /* WAIT_EVENT_H_ */
