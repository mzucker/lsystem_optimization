//////////////////////////////////////////////////////////////////////
//
// lsystems_v7.c
// Matt Zucker
//
//////////////////////////////////////////////////////////////////////
//
// Based on documentation in https://en.wikipedia.org/wiki/L-system and
// http://paulbourke.net/fractals/lsys/

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <ctype.h>
#include <math.h>
#include <time.h>

#ifdef _OPENMP
#include <omp.h>
#endif

//////////////////////////////////////////////////////////////////////
// for benchmarking

double get_time_as_double();

//////////////////////////////////////////////////////////////////////
// dynamic array

// dynamic array data type
typedef struct darray {

    size_t elem_size;
    size_t capacity;
    size_t count;
    
    unsigned char* data;

} darray_t;

// dynamic array functions
void darray_create(darray_t* darray, size_t elem_size, size_t capacity);
void darray_resize(darray_t* darray, size_t new_count);
void darray_extend(darray_t* darray, const void* elements, size_t count);
void darray_push_back(darray_t* darray, const void* elem);
void darray_pop_back(darray_t* darray, void* elem);
void* darray_elem_ptr(darray_t* darray, size_t idx);
const void* darray_const_elem_ptr(const darray_t* darray, size_t idx);
void darray_get(const darray_t* darray, size_t idx, void* dst);
void darray_set(darray_t* darray, size_t idx, const void* src);
void darray_clear(darray_t* darray);
void darray_destroy(darray_t* darray);

//////////////////////////////////////////////////////////////////////
// geometry utils

// 2D point
typedef struct point2d {
    float x, y;
} point2d_t;

// 2D rotation
typedef struct rot2d {
    float c, s;
} rot2d_t;

// 2D transformation
typedef struct xform {
    point2d_t pos;
    rot2d_t   rot;
    float     angle;
} xform_t;

// identity transform
static const xform_t IDENTITY_XFORM = {
    { 0.f, 0.f, }, { 1.f, 0.f }, 0.f
};

// 2D geometry functions
int positive_mod(int ticks, int divisor);
point2d_t rotate_point(const rot2d_t R, const point2d_t p);
rot2d_t rotate_compose(const rot2d_t R2, const rot2d_t R1);
point2d_t translate_point(const point2d_t p, const point2d_t q);
xform_t xform_inverse(xform_t xform);
xform_t xform_compose(xform_t xform2, xform_t xform1);
point2d_t xform_transform_point(xform_t xform, point2d_t p);

//////////////////////////////////////////////////////////////////////
// L-System types/functions

// misc enums
enum {
    LSYS_MAX_RULES = 128,
    LSYS_MAX_CYCLE_LENGTH = 256,
    LSYS_INIT_STRING_CAPACITY = 4096,
    LSYS_INIT_STATES_CAPACITY = 64,
    LSYS_INIT_SEGMENTS_CAPACITY = 1024
};

// line segment is 2 points
typedef struct lsys_segment {
    point2d_t p0, p1;
} lsys_segment_t;

// rule tagged with string length for string replacement
typedef struct lsys_sized_string {
    const char* replacement;
    size_t length;
} lsys_sized_string_t;

// L-System datatype
typedef struct lsystem {

    const char*         name;
    const char*         start;
    lsys_sized_string_t rules[LSYS_MAX_RULES];
    unsigned char       draw_chars[LSYS_MAX_RULES];
    double              turn_angle_rad;
    int                 rotation_cycle_length;
    rot2d_t             rotations[LSYS_MAX_CYCLE_LENGTH];
    
} lsys_t;

// lsystem character + replacement pair, for defining L-Systems
typedef struct lsys_rule_def {
    char symbol;
    const char* replacement;
} lsys_rule_def_t;

// datatype for memoizing a single L-System rule
typedef struct lsys_memo {

    size_t memo_depth;
    
    size_t segment_start;
    size_t segment_count;
    
    xform_t init_inverse;
    xform_t delta_xform;
    
} lsys_memo_t;

// datatype for memoizing an entire L-system
typedef struct lsys_memo_set {

    size_t max_depth;
    size_t min_memo_segments;
    size_t min_parallel_segments;
    
    lsys_memo_t* memos[LSYS_MAX_RULES];

} lsys_memo_set_t;

void lsys_create(lsys_t* lsys,
                 const char* name,
                 const char* start,
                 lsys_rule_def_t const rules[],
                 double turn_angle_deg,
                 const char* draw_chars);

void lsys_print(const lsys_t* lsys);

char* lsys_build_string(const lsys_t* lsys, size_t max_depth);

darray_t* lsys_segments_from_string(const lsys_t* lsys,
                                    const char* lstring);

darray_t* lsys_segments_recursive(const lsys_t* lsys,
                                  size_t max_depth,
                                  size_t min_memo_segments,
                                  size_t min_parallel_segments);

//////////////////////////////////////////////////////////////////////
// set up some known L-Systems

enum {
    LSYS_SIERPINSKI_ARROWHEAD = 0, // depth 17 takes ~5 sec
    LSYS_SIERPINSKI_TRIANGLE,      // depth 16 takes ~3 sec
    LSYS_DRAGON_CURVE,             // depth 26 takes ~6 sec
    LSYS_BARNSLEY_FERN,            // depth 13 takes ~7 sec
    LSYS_STICKS,                   // depth 16 takes ~4 sec
    LSYS_HILBERT,                  // depth 13 takes ~3 sec
    LSYS_PENTAPLEXITY,             // depth 9 takes ~2 sec
    NUM_KNOWN_LSYSTEMS
};

lsys_t KNOWN_LSYSTEMS[NUM_KNOWN_LSYSTEMS];

void initialize_known_lsystems();

//////////////////////////////////////////////////////////////////////
// options for running this program

typedef enum lsys_method {
    LSYS_METHOD_RECURSION,
    LSYS_METHOD_STRING
} lsys_method_t;

typedef struct options {
    lsys_t*       lsys;
    size_t        max_depth;
    size_t        max_segments;
    lsys_method_t method;
    size_t        min_memo_segments;
    size_t        min_parallel_segments;
} options_t;

void parse_options(int argc, char** argv, options_t* opts);

//////////////////////////////////////////////////////////////////////

double get_time_as_double() {

    struct timespec tp;

    clock_gettime(CLOCK_REALTIME, &tp);

    return (double)tp.tv_sec + (double)tp.tv_nsec * 1e-9;

}

//////////////////////////////////////////////////////////////////////

int positive_mod(int ticks, int divisor) {

    int rval = ticks % divisor;
    if (ticks < 0) {
        rval += divisor;
    }

    return rval;

}

point2d_t rotate_point(const rot2d_t R, const point2d_t p) {

    return (point2d_t){
        R.c * p.x - R.s * p.y,
        R.s * p.x + R.c * p.y
    };

}

rot2d_t rotate_compose(const rot2d_t R2, const rot2d_t R1) {

    return (rot2d_t) {
        R2.c * R1.c - R2.s * R1.s,
        R2.s * R1.c + R2.c * R1.s,
    };

}
 
point2d_t translate_point(const point2d_t p, const point2d_t q) {
    
    return (point2d_t) { p.x + q.x, p.y + q.y };

}

xform_t xform_inverse(xform_t xform) {

    rot2d_t rinv = { xform.rot.c, -xform.rot.s };

    return (xform_t) {
        .pos = rotate_point(rinv, (point2d_t){ -xform.pos.x, -xform.pos.y }),
        .rot = rinv,
        .angle = -xform.angle
    };

}

xform_t xform_compose(xform_t xform2,
                      xform_t xform1) {

    return (xform_t) {

        .pos = translate_point(xform2.pos,
                               rotate_point(xform2.rot, xform1.pos)),

        .rot = rotate_compose(xform2.rot, xform1.rot),

        .angle = xform2.angle + xform1.angle

    };
    
}

point2d_t xform_transform_point(xform_t xform,
                            point2d_t p) {

    return translate_point(rotate_point(xform.rot, p), xform.pos);

}
 
//////////////////////////////////////////////////////////////////////

void darray_create(darray_t* darray, size_t elem_size, size_t capacity) {
    
    size_t alloc_size = elem_size * capacity;
    
    darray->elem_size = elem_size;
    darray->count = 0;
    darray->capacity = capacity;
    darray->data = malloc(alloc_size);
    
}

void darray_resize(darray_t* darray, size_t new_count) {

    if (new_count > darray->capacity) {
        
        size_t new_capacity = darray->capacity;
        
        while (new_capacity <= new_count) {
            new_capacity *= 2;
        }

        size_t alloc_size = darray->elem_size * new_capacity;

        darray->data = realloc(darray->data, alloc_size);
        darray->capacity = new_capacity;
        
    }

    darray->count = new_count;
        
}

void darray_extend(darray_t* darray, const void* elements, size_t count) {

    size_t offset = darray->elem_size * darray->count;

    darray_resize(darray, darray->count + count);

    memcpy(darray->data + offset, elements, count*darray->elem_size);

}

void darray_push_back(darray_t* darray, const void* elem) {
    darray_extend(darray, elem, 1);
}

void* darray_elem_ptr(darray_t* darray, size_t idx) {
    return darray->data + idx*darray->elem_size;
}

const void* darray_const_elem_ptr(const darray_t* darray, size_t idx) {
    return darray->data + idx*darray->elem_size;
}

void darray_get(const darray_t* darray, size_t idx, void* dst) {
    memcpy(dst, darray_const_elem_ptr(darray, idx), darray->elem_size);
}

void darray_set(darray_t* darray, size_t idx, const void* src) {
    memcpy(darray_elem_ptr(darray, idx), src, darray->elem_size);
}
 
void darray_pop_back(darray_t* darray, void* dst) {
    darray->count -= 1;
    size_t offset = darray->count * darray->elem_size;
    memcpy(dst, (const void*)darray->data + offset, darray->elem_size);
}

void darray_clear(darray_t* darray) {
    darray->count = 0;
}

void darray_destroy(darray_t* darray) {
    free(darray->data);
    memset(darray, 0, sizeof(darray_t));
}

//////////////////////////////////////////////////////////////////////

void lsys_create(lsys_t* lsys,
                    const char* name,
                    const char* start,
                    lsys_rule_def_t const rules[],
                    double turn_angle_deg,
                    const char* draw_chars) {

    lsys->name = name;
    lsys->start = start;

    memset(lsys->rules, 0, sizeof(lsys->rules));
    memset(lsys->draw_chars, 0, sizeof(lsys->draw_chars));

    for (const lsys_rule_def_t* src_rule=rules; src_rule->symbol; ++src_rule) {
        lsys_sized_string_t* dst_rule = lsys->rules + (int)src_rule->symbol;
        dst_rule->replacement = src_rule->replacement;
        dst_rule->length = strlen(src_rule->replacement);
    }

    lsys->turn_angle_rad = turn_angle_deg * M_PI / 180.f;

    for (int i=0; i<=LSYS_MAX_CYCLE_LENGTH; ++i) {
        
        if (i > 0 && fmod(turn_angle_deg*i, 360.) == 0) {
            lsys->rotation_cycle_length = i;
            break;
        }
        
        float theta = lsys->turn_angle_rad * i;
        
        lsys->rotations[i].c = cosf(theta);
        lsys->rotations[i].s = sinf(theta);
        
    }

    if (draw_chars) {
        lsys->draw_chars[0] = 1;
        for (const char* c=draw_chars; *c; ++c) {
            lsys->draw_chars[(int)*c] = 1;
        }
    }

}

void lsys_print(const lsys_t* lsys) {

    printf("%s:\n", lsys->name);
    printf("  start: %s\n", lsys->start);
    printf("  rules:\n");
    for (int i=0; i<LSYS_MAX_RULES; ++i) {
        if (lsys->rules[i].replacement) {
            printf("    %c -> %s\n", i, lsys->rules[i].replacement);
        }
    }
    printf("  turn_angle_deg: %g\n", lsys->turn_angle_rad * 180.f / M_PI);
    if (lsys->rotation_cycle_length) {
        printf("  rotation_cycle_length: %d\n", lsys->rotation_cycle_length);
    }
    printf("\n");

}

char* lsys_build_string(const lsys_t* lsys, size_t max_depth) {
    
    darray_t string_darrays[2];
    
    for (int i=0; i<2; ++i) {
        darray_create(string_darrays + i, sizeof(char),
                        LSYS_INIT_STRING_CAPACITY);
    }

    int cur_idx = 0;
    
    darray_extend(string_darrays + cur_idx,
                  lsys->start,
                  strlen(lsys->start));

    for (int i=0; i<max_depth; ++i) {

        int next_idx = 1 - cur_idx;

        darray_t* src_darray = string_darrays + cur_idx;
        darray_t* dst_darray = string_darrays + next_idx;

        darray_clear(dst_darray);
        
        const char* start = (const char*)src_darray->data;
        const char* end = start + src_darray->count;

        for (const char* c=start; c!=end; ++c) {

            const lsys_sized_string_t* rule = lsys->rules + (int)*c;
            
            if (rule->replacement) {
                
                darray_extend(dst_darray,
                              rule->replacement,
                              rule->length);

            } else {

                darray_push_back(dst_darray, c);

            }

        }

        cur_idx = next_idx;
        
    }

    const char nul = '\0';
    darray_push_back(string_darrays + cur_idx, &nul);

    darray_destroy(string_darrays + (1 - cur_idx));

    return (char*)string_darrays[cur_idx].data;

}

void _lsys_execute_symbol(const lsys_t* lsys,
                          const char symbol,
                          darray_t* segments,
                          xform_t* state,
                          darray_t* xform_stack) {

    if (isalpha(symbol)) {

        if (lsys->draw_chars[0] && !lsys->draw_chars[(int)symbol]) {
            return;
        }
            
        float xnew = state->pos.x + state->rot.c;
        float ynew = state->pos.y + state->rot.s;

        lsys_segment_t seg = { { state->pos.x, state->pos.y},
                               { xnew, ynew } };

        darray_push_back(segments, &seg);

        state->pos.x = xnew;
        state->pos.y = ynew;

    } else if (symbol == '+' || symbol == '-') {

        if (lsys->rotation_cycle_length) {

            int delta = (symbol == '+') ? 1 : -1;

            int t = state->angle;

            t = positive_mod(t + delta,
                             lsys->rotation_cycle_length);

            const rot2d_t* r = lsys->rotations + t;

            state->angle = t;
            state->rot = *r;
            
        } else {

            float delta = ( (symbol == '+') ?
                            lsys->turn_angle_rad : -lsys->turn_angle_rad );
            
            state->angle += delta;

            state->rot.c = cosf(state->angle);
            state->rot.s = sinf(state->angle);

        }

    } else if (symbol == '[') {

        darray_push_back(xform_stack, state);

    } else if (symbol == ']') {

        darray_pop_back(xform_stack, state);

    } else {

        fprintf(stderr, "invalid character in string: %c\n", symbol);
        exit(1);

    }

}

darray_t* lsys_segments_from_string(const lsys_t* lsys,
                                    const char* lstring) {

    darray_t* segments = malloc(sizeof(darray_t));

    darray_create(segments, sizeof(lsys_segment_t),
                    LSYS_INIT_SEGMENTS_CAPACITY);

    darray_t xform_stack;
    darray_create(&xform_stack, sizeof(xform_t),
                    LSYS_INIT_STATES_CAPACITY);

    xform_t cur_state = IDENTITY_XFORM;

    for (const char* psymbol=lstring; *psymbol; ++psymbol) {

        _lsys_execute_symbol(lsys, *psymbol, segments,
                             &cur_state, &xform_stack);

    }

    darray_destroy(&xform_stack);

    return segments;

}

void _lsys_segments_r(const lsys_t* lsys,
                      const char* lstring, 
                      size_t remaining_depth,
                      darray_t* segments,
                      xform_t* cur_state,
                      darray_t* xform_stack,
                      lsys_memo_set_t* mset) {

    for (const char* psymbol=lstring; *psymbol; ++psymbol) {

        int symbol = *psymbol;
        
        const lsys_sized_string_t* rule = lsys->rules + symbol;

        if (remaining_depth && rule->replacement) {

            size_t segment_start = segments->count;
            xform_t xform_start = *cur_state;

            if (mset) {

                lsys_memo_t* memo = mset->memos[symbol];

                if (memo && memo->memo_depth == remaining_depth) {

                    //printf("playback %c %d\n", symbol, (int)memo->segment_count);

                    xform_t update_xform =
                        xform_compose(*cur_state, memo->init_inverse);

                    darray_resize(segments, segment_start + memo->segment_count);

                    const lsys_segment_t* src =
                        darray_const_elem_ptr(segments, memo->segment_start);

                    lsys_segment_t* dst =
                        darray_elem_ptr(segments, segment_start);

#ifdef _OPENMP                    
                    int do_parallelize =
                        (memo->segment_count >= mset->min_parallel_segments);
#endif                    
                    
                    #pragma omp parallel for if (do_parallelize)
                    for (size_t i=0; i<memo->segment_count; ++i) {
                        lsys_segment_t newsrc = {
                            xform_transform_point(update_xform, src[i].p0),
                            xform_transform_point(update_xform, src[i].p1)
                        };
                        dst[i] = newsrc;
                    }
                    
                    *cur_state = xform_compose(*cur_state, memo->delta_xform);

                    continue;

                }

            }

            _lsys_segments_r(lsys, rule->replacement,
                             remaining_depth-1,
                             segments, cur_state, xform_stack,
                             mset);

            if (mset && !mset->memos[symbol]) {

                size_t segment_count = segments->count - segment_start;

                if (segment_count > mset->min_memo_segments ||
                    remaining_depth == mset->max_depth - 2) {

                    lsys_memo_t* new_memo = malloc(sizeof(lsys_memo_t));

                    new_memo->memo_depth = remaining_depth;
                    new_memo->segment_start = segment_start;
                    new_memo->segment_count = segment_count;
                    new_memo->init_inverse = xform_inverse(xform_start);

                    new_memo->delta_xform = xform_compose(new_memo->init_inverse,
                                                          *cur_state);

                    mset->memos[symbol] = new_memo;

                }

            }

        } else {

            _lsys_execute_symbol(lsys, *psymbol, segments,
                                 cur_state, xform_stack);

        }
        
    }

}

darray_t* lsys_segments_recursive(const lsys_t* lsys,
                                  size_t max_depth,
                                  size_t min_memo_segments,
                                  size_t min_parallel_segments) {

    darray_t* segments = malloc(sizeof(darray_t));
    
    darray_create(segments, sizeof(lsys_segment_t),
                    LSYS_INIT_SEGMENTS_CAPACITY);

    darray_t xform_stack;
    darray_create(&xform_stack, sizeof(xform_t),
                    LSYS_INIT_STATES_CAPACITY);

    xform_t cur_state = IDENTITY_XFORM;

    lsys_memo_set_t mset;
    memset(&mset, 0, sizeof(mset));

    mset.max_depth = max_depth;
    mset.min_memo_segments = min_memo_segments;
    mset.min_parallel_segments = min_parallel_segments;

    _lsys_segments_r(lsys, lsys->start, 
                     max_depth, segments,
                     &cur_state, &xform_stack,
                     &mset);

    for (int i=0; i<LSYS_MAX_RULES; ++i) {
        if (mset.memos[i]) {
            free(mset.memos[i]);
        }
    }

    darray_destroy(&xform_stack);

    return segments;

}

//////////////////////////////////////////////////////////////////////
// definitions for L-Systems from websites listed at top of file

void initialize_known_lsystems() {

    lsys_create(KNOWN_LSYSTEMS + LSYS_SIERPINSKI_TRIANGLE,
                   "sierpinski_triangle", "F-G-G",
                   (lsys_rule_def_t[]){
                     { 'F', "F-G+F+G-F" },
                     { 'G', "GG" }, 
                     { 0, 0 }
                   }, 120, NULL);

    lsys_create(KNOWN_LSYSTEMS + LSYS_SIERPINSKI_ARROWHEAD,
                   "sierpinski_arrowhead", "A",
                   (lsys_rule_def_t[]){
                     { 'A', "B-A-B" }, 
                     { 'B', "A+B+A" }, 
                     { 0, 0 }
                   }, 60, NULL);
    
    lsys_create(KNOWN_LSYSTEMS + LSYS_DRAGON_CURVE,
                   "dragon_curve", "FX",
                   (lsys_rule_def_t[]){
                     { 'X', "X+YF+" },
                     { 'Y', "-FX-Y" },
                   }, 90, NULL);
    
    lsys_create(KNOWN_LSYSTEMS + LSYS_BARNSLEY_FERN,
                   "barnsley_fern", "X",
                   (lsys_rule_def_t[]){
                     { 'X', "F+[[X]-X]-F[-FX]+X" },
                     { 'F', "FF" },
                     { 0, 0 }
                   }, 25, NULL);

    lsys_create(KNOWN_LSYSTEMS + LSYS_STICKS,
                   "sticks", "X",
                   (lsys_rule_def_t[]){
                     { 'X', "F[+X]F[-X]+X" },
                     { 'F', "FF" },
                     { 0, 0 }
                   }, 20, "F");

    lsys_create(KNOWN_LSYSTEMS + LSYS_HILBERT,
                   "hilbert", "L",
                   (lsys_rule_def_t[]){
                     { 'L', "+RF-LFL-FR+" },
                     { 'R', "-LF+RFR+FL-" },
                     { 0, 0 }
                   }, 90, "F");

    lsys_create(KNOWN_LSYSTEMS + LSYS_PENTAPLEXITY,
                   "pentaplexity", "F++F++F++F++F",
                   (lsys_rule_def_t[]){
                     { 'F', "F++F++F+++++F-F++F" },
                     { 0, 0 }
                   }, 36, NULL);

}

//////////////////////////////////////////////////////////////////////

void parse_options(int argc, char** argv, options_t* opts) {

    int ok = 1;

    int disable_precomputed_rotation = 0;
    int disable_memoization = 0;
    int disable_parallelization = 0;

    memset(opts, 0, sizeof(options_t));
    opts->max_segments = 100000;

    int i=1;
    int required_count = 0;

    for (; i<argc; ++i) {

        const char* arg = argv[i];

        if (*arg && arg[0] != '-') {

            if (required_count == 0) {

                ok = 0;

                for (int j=0; j<NUM_KNOWN_LSYSTEMS; ++j) {
                    if (!strcmp(arg, KNOWN_LSYSTEMS[j].name)) {
                        ok = 1;
                        opts->lsys = KNOWN_LSYSTEMS + j;
                        break;
                    }
                }

                if (!ok) { break; }

            } else if (required_count == 1) {

                int d;
            
                if (sscanf(arg, "%d", &d) != 1 || d <= 0) {
                    ok = 0;
                    break;
                }

                opts->max_depth = d;

            } else {

                ok = 0;
                break;

            }

            ++required_count;
            
        } else if (!strcmp(arg, "-s")) {
            
            opts->method = LSYS_METHOD_STRING;
            
        } else if (!strcmp(arg, "-r")) {

            opts->method = LSYS_METHOD_RECURSION;
            
        } else if (!strcmp(arg, "-x")) {
            
            if (++i == argc) { ok = 0; break; }
            
            int d;
            
            if (sscanf(argv[i], "%d", &d) != 1) {
                ok = 0; break;
            }
            
            if (d >= -1) {
                opts->max_segments = (size_t)d;
            } else {
                ok = 0;
                break;
            }

        } else if (!strcmp(arg, "-M")) {

            disable_memoization = 1;
            
#ifdef _OPENMP            
        } else if (!strcmp(arg, "-P")) {
            
            disable_parallelization = 1;
#endif            
        } else if (!strcmp(arg, "-R")) {
            
            disable_precomputed_rotation = 1;

        } else {

            fprintf(stderr, "error: unrecognized option %s\n\n", arg);

            ok = 0;
            break;

        }

    }

    if (!ok || !opts->lsys || !opts->max_depth) {
        printf("usage: %s [options] LSYSTEM MAXDEPTH\n"
               "\n"
               "where LSYSTEM is one of:\n", argv[0]);
        for (int j=0; j<NUM_KNOWN_LSYSTEMS; ++j) {
            printf("  * %s\n", KNOWN_LSYSTEMS[j].name);
        }
        printf("\n");
        printf("options:\n");
        printf("  -x MAXSEGMENTS maximum number of segments for output\n"
               "  -s             use string building method\n"
               "  -r             use recursive method (default)\n"
               "  -M             disable memoization for recursive method\n"
#ifdef _OPENMP               
               "  -P             disable parallelization for memoization\n"
#endif               
               "  -R             don't precompute rotations\n"
               "\n");
        exit(1);
    }

    printf("using %s method\n",
           opts->method == LSYS_METHOD_STRING ? "string" : "recursion");
    
    if (opts->method == LSYS_METHOD_STRING) {
        if (disable_memoization) {
            printf("warning: disabling memoization has no effect for string method!\n");
        }
    }

    int have_memo = (opts->method == LSYS_METHOD_RECURSION) && !disable_memoization;
    
    if (!have_memo) {
        if (disable_parallelization) {
            printf("warning: disabling parallelization has no effect when not memoizing\n");
        }
        opts->min_memo_segments = (size_t)-1;
        opts->min_parallel_segments = (size_t)-1;
        if (opts->method == LSYS_METHOD_RECURSION) {
            printf("memoization is disabled\n");
        } 
    } else if (disable_parallelization) {
        opts->min_memo_segments = 10000;
        opts->min_parallel_segments = (size_t)-1;
        printf("memoizing runs with > %d segments\n",
               (int)opts->min_memo_segments);
    } else {
        opts->min_memo_segments = 100000;
        opts->min_parallel_segments = 5000;
        printf("memoizing runs with > %d segments and parallelizing when > %d segments\n",
               (int)opts->min_memo_segments, (int)opts->min_parallel_segments);
    }

    if (disable_precomputed_rotation) {
        printf("disabling precomputed rotation!\n");
        opts->lsys->rotation_cycle_length = 0;
    }

    printf("\n");

    lsys_print(opts->lsys);
    
}

//////////////////////////////////////////////////////////////////////
// main program

int main(int argc, char** argv) {

    // initialize lsystems
    initialize_known_lsystems();

    // parse command-line options
    options_t opts;
    parse_options(argc, argv, &opts);

    ////////////////////////////////////////////////////////////
    // now get the segments

    double start = get_time_as_double();

    darray_t* segments;

    if (opts.method == LSYS_METHOD_STRING) {

        char* lstring = lsys_build_string(opts.lsys,
                                             opts.max_depth);

        segments = lsys_segments_from_string(opts.lsys, lstring);

        free(lstring);

    } else {

        segments = lsys_segments_recursive(opts.lsys, opts.max_depth,
                                           opts.min_memo_segments,
                                           opts.min_parallel_segments);

    }

    double elapsed = get_time_as_double() - start;
    printf("generated %d segments in %.6f s (%.3f ns/segment).\n",
           (int)segments->count, elapsed, 1e9 * elapsed/segments->count);

    ////////////////////////////////////////////////////////////
    // either output or not

    if (segments->count > opts.max_segments) {
        
        printf("...maximum of %d segments exceeded, skipping output!\n",
               (int)opts.max_segments);

    } else {

        const lsys_segment_t* segment = (const lsys_segment_t*)segments->data;
        const lsys_segment_t* end = segment + segments->count;

        FILE* outfile = fopen("segments.txt", "w");
        
        for ( ; segment != end; ++segment) {
            fprintf(outfile, "%g %g %g %g\n",
                    segment->p0.x, segment->p0.y,
                    segment->p1.x, segment->p1.y);
        }

        fclose(outfile);

        printf("wrote segments.txt\n");

    }

    ////////////////////////////////////////////////////////////
    // clean up

    darray_destroy(segments);
    free(segments);

    return 0;

}
