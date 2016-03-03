#include <stdio.h>

#include "cubs/spline.h"
#include <math.h>


int main(){
    #define tau 3.1415926539 * 2

    cubs_spline_t spline = cubs_create_spline( 2, 0, 10 );
    cubs_cache_t cache = cubs_create_cache(1000000, 1);
    double t = 0;
    for( t = 0; t <= tau; t += tau / 200.0 ){
        float point[2] = {cos(t), sin(t)};
        cubs_spline_add(spline, point);
    }

    cubs_spline_compile( spline, cache );
    //cubs_spline_bake( spline );
    spline = cubs_reference_spline(spline);
    cubs_destroy_spline(spline);
    printf("recompiling\n");
    cubs_spline_compile( spline, cache );
    for( t = -10; t <= 20; t += 0.2 ){
        float point[2];
        cubs_spline_evaluate(spline, t, point );
        printf("(%f, %f)\n", point[0], point[1] );
    }
    cubs_destroy_spline(spline);
    cubs_destroy_cache(cache);
}
