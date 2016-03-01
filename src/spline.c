#include "cubs/spline.h"
#include <stdlib.h>

struct spline_data{
	cubs_number_t *point_data;
	cubs_index_t size, reserve, references;
	cubs_number_t minimum, maximum;
	cubs_number_t dimensions;
};

struct spline_matrix{
	cubs_number_t *inv_mat;
	cubs_index_t length;
	cubs_number_t divisor;
};

struct spline_cache{
	struct spline_matrix *matrices; //sorted on length
	cubs_index_t length, references, max_size;
};

static cubs_index_t cubs_spline_mem_calc(cubs_index_t count, cubs_index_t dimensions){
	return count * dimensions * sizeof( cubs_number_t );
}

cubs_spline_t cubs_create_spline(
		cubs_index_t dimensions, 
		cubs_number_t minimum, //Time value of first element
		cubs_number_t maximum //Time value of last element
)
{
	cubs_spline_t spline = malloc(sizeof(struct spline_data));
	spline->size = 0;
	spline->reserve = 100;
	spline->minimum = minimum;
	spline->maximum = maximum;
	spline->dimensions = dimensions;
	spline->references = 1;
	spline->point_data = malloc(cubs_spline_mem_calc(spline->reserve, spline->dimensions));
	return spline;
}

cubs_spline_t cubs_reference_spline(
	cubs_spline_t other
)
{
	if( !other ){
		return other;
	}
	other->references ++;
	return other;
}

//Reduce reference count, and destroy if references reaches zero.
void cubs_destroy_spline(
	cubs_spline_t other
)
{
	if( !other ){
		return;
	}
	other->references --;
	if( other->references == 0 ){
		if( other->point_data ){
			free( other->point_data );
		}
		free( other );
	}
}


struct spline_matrix *cubs_create_matrix(
	struct spline_matrix *base, 
	cubs_index_t size
)
{
	
}

void cubs_destroy_matrix(
	struct spline_matrix *matrix
)
{
	free( matrix );
}

cubs_cache_t cubs_create_cache(
	cubs_index_t max_size
)
{
	cubs_cache_t cache = malloc(sizeof(struct spline_cache));
	cache->matrices = NULL;
	cache->length = 0;
	cache->max_size = max_size;
	cache->references = 1;
}

cubs_cache_t cubs_reference_cache(
	cubs_cache_t other
)
{
	if( !other ){
		return other;
	}
	other->references ++;
	return other;
}

//Reduce reference count, and destroy if references reaches zero.
void cubs_destroy_cache(
	cubs_cache_t other
)
{
	if( !other ){
		return;
	}
	other->references --;
	if( other->references == 0 ){
		if( other->matrices ){
			free( other->matrices );
		}
		free( other );
	}
}


//Returns index of added point
cubs_index_t cubs_spline_add(
	cubs_spline_t spline, 
	cubs_point_t point
);

//Returns index of added point
cubs_index_t cubs_spline_insert(
	cubs_spline_t spline, 
	cubs_index_t index,
	cubs_point_t point
);

void cubs_spline_erase(
	cubs_spline_t spline, 
	cubs_index_t index
);

cubs_index_t cubs_spline_size(
	cubs_spline_t spline
);

//Get mutable point
cubs_point_t cubs_spline_get(
	cubs_spline_t spline, 
	cubs_index_t index
);

//Change point at index. Alternative to mutating point obtained from _get
void cubs_spline_mutate(
	cubs_spline_t spline, 
	cubs_index_t index,
	cubs_point_t point
);

//Prepare spline for evaluation
void cubs_spline_compile(
	cubs_spline_t spline,
	cubs_cache_t cache OPTIONAL
);

//Evaluate spline at time t
void cubs_spline_evaluate(
	cubs_spline_t spline,
	cubs_number_t t,
	cubs_point_t point OUT
);

//Evaluate first derivative of spline at time t
void cubs_spline_evaluate_dv1(
	cubs_spline_t spline,
	cubs_number_t t,
	cubs_point_t point OUT
);

//Evaluate second derivative of spline at time t
void cubs_spline_evaluate_dv2(
	cubs_spline_t spline,
	cubs_number_t t,
	cubs_point_t point OUT
);

