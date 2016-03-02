#ifndef INC_CUBS_SPLINE_H
#define INC_CUBS_SPLINE_H

#ifndef CUBS_DATA_T
#define CUBS_DATA_T float
#endif

#define OPTIONAL
#define OUT
#define IN_OUT

typedef CUBS_DATA_T cubs_number_t;
typedef struct spline_data* cubs_spline_t;
typedef cubs_number_t* cubs_point_t;
typedef unsigned int cubs_index_t;
typedef struct spline_cache* cubs_cache_t;

cubs_spline_t cubs_create_spline(
		cubs_index_t dimensions,
		cubs_number_t minimum, //Time value of first element
		cubs_number_t maximum //Time value of last element
);

cubs_spline_t cubs_reference_spline(
	cubs_spline_t other
);

//Reduce reference count, and destroy if references reaches zero.
void cubs_destroy_spline(
	cubs_spline_t other
);

cubs_cache_t cubs_create_cache(
	cubs_index_t max_size
);

cubs_cache_t cubs_reference_cache(
	cubs_cache_t other
);

//Reduce reference count, and destroy if references reaches zero.
void cubs_destroy_cache(
	cubs_cache_t other
);


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

//Get dimensionality of spline
cubs_index_t cubs_spline_dimensions(
	cubs_spline_t spline
);

//Get mutable point
//Modifying a mutable point after altering the internal layout of the spline
//(compiling, adding, deleting nodes) will cause bad things to happen...
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

//Reduce a compiled spline to minimum size (cannot be edited)
void cubs_spline_bake(
	cubs_spline_t spline
);

#endif

