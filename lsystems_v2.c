//////////////////////////////////////////////////////////////////////
//
// lsystems_v2.c
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

//////////////////////////////////////////////////////////////////////
// for benchmarking

double get_time_as_double(void);

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
    size_t      length;
} lsys_sized_string_t;

// L-System datatype
typedef struct lsystem {

    const char*         name;
    const char*         start;
    lsys_sized_string_t rules[LSYS_MAX_RULES];
    unsigned char       draw_chars[LSYS_MAX_RULES];
    double              turn_angle_rad;

} lsys_t;

// lsystem character + replacement pair, for defining L-Systems
typedef struct lsys_rule_def {
    char        symbol;
    const char* replacement;
} lsys_rule_def_t;

typedef struct lsys_state {
    point2d_t pos;
    float     angle;
} lsys_state_t;

static const lsys_state_t LSYS_START_STATE = { { 0.f, 0.f, }, 0.f };

void lsys_create(lsys_t* lsys,
                 const char* name,
                 const char* start,
                 lsys_rule_def_t const rules[],
                 double turn_angle_deg,
                 const char* draw_chars);

void lsys_print(const lsys_t* lsys);

char* lsys_build_string(const lsys_t* lsys, size_t total_iterations);

darray_t* lsys_segments_from_string(const lsys_t* lsys,
                                    const char* lstring);

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

void initialize_known_lsystems(void);

//////////////////////////////////////////////////////////////////////
// options for running this program

typedef struct options {
    lsys_t* lsys;
    size_t  total_iterations;
    size_t  max_segments;
} options_t;

void parse_options(int argc, char** argv, options_t* opts);

//////////////////////////////////////////////////////////////////////

double get_time_as_double(void) {

    struct timespec tp;

    clock_gettime(CLOCK_REALTIME, &tp);

    return (double)tp.tv_sec + (double)tp.tv_nsec * 1e-9;

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
    printf("\n");

}

char* lsys_build_string(const lsys_t* lsys, size_t total_iterations) {

    darray_t string_darrays[2];

    for (int i=0; i<2; ++i) {
        darray_create(string_darrays + i, sizeof(char),
                      LSYS_INIT_STRING_CAPACITY);
    }

    int cur_idx = 0;

    darray_extend(string_darrays + cur_idx,
                  lsys->start,
                  strlen(lsys->start));

    for (int i=0; i<total_iterations; ++i) {

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
                          lsys_state_t* state,
                          darray_t* state_stack) {

    if (isalpha(symbol)) {

        if (lsys->draw_chars[0] && !lsys->draw_chars[(int)symbol]) {
            return;
        }

        float c = cosf(state->angle);
        float s = sinf(state->angle);

        float xnew = state->pos.x + c;
        float ynew = state->pos.y + s;

        lsys_segment_t seg = { { state->pos.x, state->pos.y},
                               { xnew, ynew } };

        darray_push_back(segments, &seg);

        state->pos.x = xnew;
        state->pos.y = ynew;

    } else if (symbol == '+' || symbol == '-') {

        float delta = ( (symbol == '+') ?
                        lsys->turn_angle_rad : -lsys->turn_angle_rad );

        state->angle += delta;

    } else if (symbol == '[') {

        darray_push_back(state_stack, state);

    } else if (symbol == ']') {

        darray_pop_back(state_stack, state);

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

    darray_t state_stack;
    darray_create(&state_stack, sizeof(lsys_state_t),
                  LSYS_INIT_STATES_CAPACITY);

    lsys_state_t cur_state = LSYS_START_STATE;

    for (const char* psymbol=lstring; *psymbol; ++psymbol) {

        _lsys_execute_symbol(lsys, *psymbol, segments,
                             &cur_state, &state_stack);

    }

    darray_destroy(&state_stack);

    return segments;

}

//////////////////////////////////////////////////////////////////////
// definitions for L-Systems from websites listed at top of file

void initialize_known_lsystems(void) {

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

                opts->total_iterations = d;

            } else {

                ok = 0;
                break;

            }

            ++required_count;

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

        } else {

            fprintf(stderr, "error: unrecognized option %s\n\n", arg);

            ok = 0;
            break;

        }

    }

    if (!ok || !opts->lsys || !opts->total_iterations) {
        printf("usage: %s [options] LSYSTEM ITERATIONS\n"
               "\n"
               "where LSYSTEM is one of:\n", argv[0]);
        for (int j=0; j<NUM_KNOWN_LSYSTEMS; ++j) {
            printf("  * %s\n", KNOWN_LSYSTEMS[j].name);
        }
        printf("\n");
        printf("options:\n");
        printf("  -x MAXSEGMENTS maximum number of segments for output\n"
               "\n");
        exit(1);
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

    char* lstring = lsys_build_string(opts.lsys,
                                      opts.total_iterations);

    segments = lsys_segments_from_string(opts.lsys, lstring);

    free(lstring);

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
