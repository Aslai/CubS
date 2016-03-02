#include "cubs/spline.h"
#include <stdlib.h>
#include <string.h>
#include <math.h>

typedef struct spline_cubic_func {
    cubs_number_t a, b, c, d;
} cubs_func_t;

struct spline_data {
    cubs_func_t *functions;
    cubs_index_t func_size;
    cubs_number_t minimum, maximum;
    cubs_index_t dimensions, references;
    char compiled;
    cubs_number_t *point_data;
    cubs_index_t size, reserve;
};

typedef struct spline_matrix {
    cubs_number_t *inv_mat;
    cubs_index_t length, dist;
    cubs_number_t divisor;
} *cubs_matrix_t;

struct spline_cache {
    struct spline_matrix *matrices; //sorted on length
    cubs_index_t length, references, max_size;
};


static cubs_index_t cubs_spline_mem_calc( cubs_index_t count, cubs_index_t dimensions ) {
    return count * dimensions * sizeof( cubs_number_t );
}

cubs_spline_t cubs_create_spline( cubs_index_t dimensions, cubs_number_t minimum, cubs_number_t maximum ) {
    cubs_spline_t spline = malloc(sizeof(struct spline_data));
    spline->size = 0;
    spline->reserve = 100;
    spline->minimum = minimum;
    spline->maximum = maximum;
    spline->dimensions = dimensions;
    spline->references = 1;
    spline->point_data = malloc(cubs_spline_mem_calc(spline->reserve, spline->dimensions));
    spline->compiled = 0;
    spline->functions = NULL;
    spline->func_size = 0;
    return spline;
}

cubs_spline_t cubs_reference_spline( cubs_spline_t other ) {
    if( !other ) {
        return other;
    }
    other->references ++;
    return other;
}

//Reduce reference count, and destroy if references reaches zero.
void cubs_destroy_spline( cubs_spline_t other ) {
    if( !other ) {
        return;
    }
    other->references --;
    if( other->references == 0 ) {
        if( other->point_data ) {
            free( other->point_data );
        }
        free( other );
    }
}

//Does not do any allocation, expects preallocated pointers
struct spline_matrix *cubs_create_matrix( struct spline_matrix *matrix, cubs_index_t size ) {
    cubs_index_t i, div;
    matrix->inv_mat = malloc(sizeof(cubs_number_t) * size + 2);
    matrix->length = size;
    matrix->divisor = 1;
    matrix->inv_mat[0] = 0;
    matrix->inv_mat[1] = 1;
    //Generate first row of [A]^-1
    for( i = 2; i < size; ++i ) {
        matrix->inv_mat[i] = matrix->inv_mat[i - 1] * 4 - matrix->inv_mat[i - 2];
        //Keep floats from overflowing to infinity
        if( matrix->inv_mat[i] > 1e8 ) {
            for( div = 0; div <= i; ++div ) {
                matrix->inv_mat[div] /= matrix->inv_mat[i];
            }
        }
    }
    for( div = 0; div <= size; ++div ) {
        matrix->inv_mat[div] /= matrix->inv_mat[size - 1];
    }
    matrix->inv_mat[size] = 0;
    matrix->dist = size;
    for( div = 0; div <= size; ++div ) {
        if( matrix->inv_mat[size - 2 - div] == 0.0 ) {
            matrix->dist = sqrt(div * div + div * div) / 2;
            matrix->dist += 1;
            break;
        }
    }
    return matrix;
}

//Assumes external free for matrix itself
void cubs_destroy_matrix( struct spline_matrix *matrix ) {
    if( matrix && matrix->inv_mat ) {
        free( matrix->inv_mat );
    }
}

cubs_cache_t cubs_create_cache( cubs_index_t max_size ) {
    cubs_cache_t cache = malloc(sizeof(struct spline_cache));
    cache->matrices = NULL;
    cache->length = 0;
    cache->max_size = max_size;
    cache->references = 1;
}

cubs_cache_t cubs_reference_cache( cubs_cache_t other ) {
    if( !other ) {
        return other;
    }
    other->references ++;
    return other;
}

//Reduce reference count, and destroy if references reaches zero.
void cubs_destroy_cache( cubs_cache_t other ) {
    cubs_index_t i;
    if( !other ) {
        return;
    }
    other->references --;
    if( other->references == 0 ) {
        if( other->matrices ) {
            for( i = 0; i < other->length; ++i ) {
                cubs_destroy_matrix(other->matrices + i );
            }
            free( other->matrices );
        }
        free( other );
    }
}

cubs_number_t cubs_matrix_get( cubs_matrix_t matrix, cubs_index_t i, cubs_index_t j ) {
    cubs_index_t pos = abs((int)i - (int)j);
    if( pos > matrix->dist ) {
        return 0;
    }
    cubs_number_t value = 0;
    cubs_index_t m = matrix->length - 2;

    cubs_index_t effectiveSize = m / 2 + m % 2;
    cubs_index_t effectiveI = i;
    cubs_index_t effectiveJ = j;
    if (i >= effectiveSize) {
        effectiveI = m - 1 - i;
    }
    if (j >= effectiveSize) {
        effectiveJ = m - 1 - j;
    }
    cubs_index_t terms = effectiveI;
    if (effectiveJ < terms) {
        terms = effectiveJ;
    }

    //Repeatedly add terms from the sequence
    cubs_index_t t;
    for (t = 0; t <= terms; ++t) {
        value += matrix->inv_mat[matrix->length - 2 - (pos + 2 * t)];
    }
    return value * ((i + j) % 2 ? -1.0 : 1.0);
}

struct spline_matrix *cubs_cache_add( cubs_cache_t cache, cubs_index_t length ) {
    cubs_index_t idx = 0;

    if( cache->length != 0 ) {
        //Do a binary search to find a cached version of this matrix
        cubs_index_t low, high, mid;
        low = 0;
        high = cache->length - 1;
        mid = high / 2;
        while( low != high ) {
            if( length > cache->matrices[mid].length ) {
                if( mid == cache->length - 1 ) {
                    mid++;
                    break;
                }
                low = mid + 1;
            } else if( length < cache->matrices[mid].length ) {
                if( mid == 0 ) {
                    break;
                }
                high = mid - 1;
            } else {
                return &cache->matrices[mid];
            }
            mid = low + (high - low) / 2;
        }
        idx = mid;
    }
    if( cache->matrices == NULL ) {
        cache->matrices = malloc( sizeof( struct spline_matrix ) );
        cache->length = 1;
        return cubs_create_matrix( cache->matrices, length );
    }
    //We need to realloc and generate the new matrix at location idx
    cache->matrices = realloc( cache->matrices, (cache->length + 1) * sizeof( struct spline_matrix ) );
    memmove( cache->matrices + idx + 1, cache->matrices + idx, cache->length - idx );
    cache->length ++;
    return cubs_create_matrix( cache->matrices + idx, length );
}

//Returns index of added point
cubs_index_t cubs_spline_add( cubs_spline_t spline, cubs_point_t point ) {
    spline->compiled = 0;

    if( spline->size + 1 >= spline->reserve ) {
        spline->reserve *= 2;
        spline->point_data = realloc(   spline->point_data,
                                        cubs_spline_mem_calc(spline->reserve, spline->dimensions) );
    }
    memcpy(     spline->point_data + spline->size * spline->dimensions,
                point,
                spline->dimensions * sizeof( cubs_number_t ) );
    return spline->size ++;
}

//Returns index of added point
cubs_index_t cubs_spline_insert( cubs_spline_t spline, cubs_index_t index, cubs_point_t point ) {
    spline->compiled = 0;

    if( spline->size + 1 >= spline->reserve ) {
        spline->reserve *= 2;
        spline->point_data = realloc(   spline->point_data,
                                        cubs_spline_mem_calc(spline->reserve, spline->dimensions) );
    }

    memmove(    spline->point_data + (index + 1) * spline->dimensions,
                spline->point_data + index * spline->dimensions,
                spline->dimensions * sizeof( cubs_number_t ) * (spline->size - index) );

    memcpy(     spline->point_data + index * spline->dimensions,
                point,
                spline->dimensions * sizeof( cubs_number_t ) );

    spline->size ++;
    return index;
}

void cubs_spline_erase( cubs_spline_t spline, cubs_index_t index ) {
    spline->compiled = 0;

    memmove(    spline->point_data + index * spline->dimensions,
                spline->point_data + (index + 1) * spline->dimensions,
                spline->dimensions * sizeof( cubs_number_t ) * (spline->size - index - 1) );
    spline->size --;
}

cubs_index_t cubs_spline_size( cubs_spline_t spline ) {
    return spline->size;
}

//Get mutable point
cubs_point_t cubs_spline_get( cubs_spline_t spline, cubs_index_t index ) {
    spline->compiled = 0;
    return spline->point_data + index * spline->dimensions;
}

//Change point at index. Alternative to mutating point obtained from _get
void cubs_spline_mutate( cubs_spline_t spline, cubs_index_t index, cubs_point_t point ) {
    memcpy(     spline->point_data + index * spline->dimensions,
                point,
                spline->dimensions * sizeof( cubs_number_t ) );
}

//Prepare spline for evaluation
//Should be split up a bit more.
void cubs_spline_compile( cubs_spline_t spline, cubs_cache_t cache OPTIONAL ) {
    cubs_index_t n = spline->size;
    cubs_index_t m = n - 2;
    cubs_index_t i, j, k;
    if( spline->size <= 1 ) {
        return;
    }

    cubs_matrix_t matrix;
    struct spline_matrix non_cached;
    if( cache ) {
        matrix = cubs_cache_add( cache, spline->size );
    } else {
        matrix = cubs_create_matrix( &non_cached, spline->size );
    }

    cubs_number_t *C = malloc( (n - 2) * sizeof(cubs_number_t) * spline->dimensions );
    cubs_number_t *C2 = malloc( n * sizeof(cubs_number_t) * spline->dimensions );
    cubs_number_t h = (spline->maximum - spline->minimum) / (n - 1);

    for (i = 1; i < n - 1; ++i) {
        for (j = 0; j < spline->dimensions; ++j) {
            const cubs_number_t p1 = spline->point_data[(i - 1) * spline->dimensions + j];
            const cubs_number_t p2 = spline->point_data[(i) * spline->dimensions + j];
            const cubs_number_t p3 = spline->point_data[(i + 1) * spline->dimensions + j];
            const cubs_number_t constant = 1; // 6 / (h * h);

            C[(i - 1) * spline->dimensions + j] = constant * (p1 - 2 * p2 + p3);
        }
    }

    for( i = 0; i < spline->dimensions * n; ++i ) {
        C2[i] = 0;
    }
    for ( i = 0; i < m; ++i) {
        cubs_number_t *value = C2 + (i+1) * spline->dimensions;
        for ( j = 0; j < m; ++j) {
            cubs_number_t mat = cubs_matrix_get( matrix, i, j );
            for (k = 0; k < spline->dimensions; ++k) {
                value[k] += mat * C[j * spline->dimensions + k];
            }
        }
    }
    free(C);
    C = C2;
    for( i = 0; i < n * spline->dimensions; ++i ) {
        C[i] *= 6.0 / (h * h);
    }

    if( spline->functions ) {
        free(spline->functions);
    }
    spline->func_size = n - 1;
    spline->functions = malloc( spline->func_size * spline->dimensions * sizeof(cubs_func_t));

    //For all the points, generate the functions with appropriate coefficients for the
    //piecewise function.
    for ( i = 0; i < n - 1; ++i ) {
        for ( j = 0; j < spline->dimensions; ++j ) {
            cubs_func_t *f = spline->functions + i * spline->dimensions + j;
            const cubs_index_t first = i * spline->dimensions + j;
            const cubs_index_t second = (i+1) * spline->dimensions + j;
            const cubs_number_t first_value = spline->point_data[first];
            const cubs_number_t second_value = spline->point_data[second];


            f->a = (C[second] - C[first]) / (6.0 * h);
            f->b = C[first] / 2;
            f->c = (second_value - first_value) / h - h * (C[second] + 2.0 * C[first]) / 6.0;
            f->d = first_value;
        }
    }
    free( C );

    if( !cache ) {
        cubs_destroy_matrix( matrix );
    }
}

static void cubs_evaluate_params( cubs_spline_t spline, cubs_number_t *t IN_OUT, cubs_index_t *idx OUT ) {
    cubs_number_t range = spline->maximum - spline->minimum;
    cubs_number_t frac = *t / range;
    if( frac < 0 ) {
        frac = 0;
    } else if( frac > 1 ) {
        frac = 1;
    }
    *idx = spline->func_size * frac;
    if( *idx >= spline->func_size ) {
        *idx = spline->func_size - 1;
    }

    *t -= spline->minimum;
    *t -= (cubs_number_t)*idx / (cubs_number_t)spline->func_size * range;
}

//Evaluate spline at time t
void cubs_spline_evaluate( cubs_spline_t spline, cubs_number_t t, cubs_point_t point OUT ) {
    cubs_index_t idx;
    cubs_evaluate_params( spline, &t, &idx );
    cubs_func_t *funcs = spline->functions + idx * spline->dimensions;
    for( idx = 0; idx < spline->dimensions; ++idx ) {
        point[idx] = funcs[idx].a * t * t * t + funcs[idx].b * t * t + funcs[idx].c * t + funcs[idx].d;
    }
}

//Evaluate first derivative of spline at time t
void cubs_spline_evaluate_dv1( cubs_spline_t spline, cubs_number_t t, cubs_point_t point OUT ) {
    cubs_index_t idx;
    cubs_evaluate_params( spline, &t, &idx );
    cubs_func_t *funcs = spline->functions + idx * spline->dimensions;
    for( idx = 0; idx < spline->dimensions; ++idx ) {
        point[idx] = 3 * funcs[idx].a * t * t + 2 * funcs[idx].b * t + funcs[idx].c;
    }
}

//Evaluate second derivative of spline at time t
void cubs_spline_evaluate_dv2( cubs_spline_t spline, cubs_number_t t, cubs_point_t point OUT ) {
    cubs_index_t idx;
    cubs_evaluate_params( spline, &t, &idx );
    cubs_func_t *funcs = spline->functions + idx * spline->dimensions;
    for( idx = 0; idx < spline->dimensions; ++idx ) {
        point[idx] = 6 * funcs[idx].a * t + 2 * funcs[idx].b;
    }
}

//Reduce a compiled spline to minimum size (cannot be edited)
void cubs_spline_bake( cubs_spline_t spline ) {
    if( spline->point_data != NULL ) {
        free( spline->point_data );
        spline->point_data = NULL;
    }
    spline->size = spline->reserve = 0;
}
