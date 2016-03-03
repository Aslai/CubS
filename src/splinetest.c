#include <stdio.h>

#include "cubs/spline.h"
#include <math.h>
#include <time.h>

void stopwatch_start( struct timespec *t ){
    clock_gettime(CLOCK_THREAD_CPUTIME_ID, t);
}
double stopwatch_time( struct timespec *t ){
    struct timespec tr, tn;
    clock_getres(CLOCK_THREAD_CPUTIME_ID, &tr);
    clock_gettime(CLOCK_THREAD_CPUTIME_ID, &tn);
    double s = (tn.tv_sec - t->tv_sec);
    s *= 1000000000;
    s += tn.tv_nsec;
    s -= t->tv_nsec;
    s /= 1000000000;
    return s;
}



int main(){
    #define tau 3.1415926539 * 2
    struct timespec timer;

    cubs_spline_t spline = cubs_create_spline( 2, 0, 10 );
    cubs_cache_t cache = cubs_create_cache(1000000, 1);
    double t = 0;
    for( t = 0; t <= tau; t += tau / 5000.0 ){
        float point[2] = {cos(t), sin(t)};
        cubs_spline_add(spline, point);
    }

    stopwatch_start(&timer);
    cubs_spline_compile( spline, cache );
    printf("Initial compile: % 15.10f seconds\n", stopwatch_time(&timer) );
    //cubs_spline_bake( spline );
    spline = cubs_reference_spline(spline);
    cubs_destroy_spline(spline);
    stopwatch_start(&timer);
    cubs_spline_compile( spline, cache );
    printf("      Recompile: % 15.10f seconds\n", stopwatch_time(&timer) );

    stopwatch_start(&timer);
    for( t = -10; t <= 20; t += 0.2 ){
        float point[2];
        cubs_spline_evaluate(spline, t, point );
        printf("(%f, %f)\n", point[0], point[1] );
    }
    printf("     Evaluation: % 15.10f seconds\n", stopwatch_time(&timer) );

    cubs_destroy_spline(spline);
    cubs_destroy_cache(cache);
}
