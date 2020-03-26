#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <ctype.h>
#include <math.h>
#include <time.h>

typedef struct buffer {

    size_t elem_size;
    size_t capacity;
    size_t count;
    
    unsigned char* data;

} buffer_t;

//////////////////////////////////////////////////////////////////////

enum {
    MAX_RULES = 128,
    MAX_CYCLE_LENGTH = 256,
    INIT_STRING_CAPACITY = 4096,
    INIT_STATES_CAPACITY = 64,
    INIT_SEGMENTS_CAPACITY = 1024
};

typedef struct vec2d {
    float x, y;
} vec2d_t;

typedef struct rot2d {
    float c, s;
} rot2d_t;

typedef struct xform {
    vec2d_t pos;
    rot2d_t rot;
    float   turn;
} xform_t;

static const xform_t IDENTITY_XFORM = {
    { 0.f, 0.f, }, { 1.f, 0.f }, 0.f
};

typedef struct lsystem_segment {
    vec2d_t p0, p1;
} lsystem_segment_t;

typedef struct lsystem_rule_tagged {
    
    const char* replacement;
    size_t length;
    
} lsystem_rule_tagged_t;

typedef struct lsystem {

    const char*           name;
    const char*           start;
    lsystem_rule_tagged_t rules[MAX_RULES];
    unsigned char         draw_chars[MAX_RULES];
    double                turn_angle_rad;
    int                   rotation_cycle_length;
    rot2d_t               rotations[MAX_CYCLE_LENGTH];
    
} lsystem_t;

typedef struct lsystem_rule {
    
    char variable;
    const char* replacement;
    
} lsystem_rule_t;

typedef struct lsystem_memo {
    
    size_t segment_start;
    size_t segment_count;
    
    xform_t init_inverse;
    xform_t delta_xform;
    
} lsystem_memo_t;

typedef struct lsystem_memo_set {

    size_t memo_depth;
    lsystem_memo_t* memos[MAX_RULES];

} lsystem_memo_set_t;


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

lsystem_t known_lsystems[NUM_KNOWN_LSYSTEMS];

//////////////////////////////////////////////////////////////////////

typedef enum method {
    METHOD_STRING,
    METHOD_RECURSION
} method_t;

typedef struct options {
    lsystem_t* lsys;
    size_t     max_depth;
    size_t     memo_depth;
    method_t   method;
    size_t     max_segments;
} options_t;

//////////////////////////////////////////////////////////////////////

int positive_mod(int ticks, int divisor) {

    int rval = ticks % divisor;
    if (ticks < 0) {
        rval += divisor;
    }

    return rval;

}

vec2d_t rotate(const rot2d_t R, const vec2d_t p) {

    vec2d_t result = {
        R.c * p.x - R.s * p.y,
        R.s * p.x + R.c * p.y
    };

    return result;

}

rot2d_t rotateR(const rot2d_t R2, const rot2d_t R1) {

    rot2d_t result = {
        R2.c * R1.c - R2.s * R1.s,
        R2.s * R1.c + R2.c * R1.s,
    };

    return result;

}

vec2d_t translate(const vec2d_t p, const vec2d_t q) {
    vec2d_t result = { p.x + q.x, p.y + q.y };
    return result;
}

xform_t xform_inverse(xform_t xform) {

    xform_t rval = {
        .pos = { -xform.pos.x, -xform.pos.y },
        .rot = {  xform.rot.c, -xform.rot.s },
        .turn = -xform.turn
    };

    rval.pos = rotate(rval.rot, rval.pos);

    return rval;

}

xform_t xform_compose(xform_t xform2,
                      xform_t xform1) {

    xform_t rval = {

        .pos = translate(xform2.pos,
                         rotate(xform2.rot, xform1.pos)),

        .rot = rotateR(xform2.rot, xform1.rot),

        .turn = xform2.turn + xform1.turn

    };
    
    return rval;

}

vec2d_t xform_transform_point(xform_t xform,
                             vec2d_t p) {

    return translate(rotate(xform.rot, p), xform.pos);

}

void xform_print(const char* name, xform_t xform) {

    printf("%s x=%.2f, y=%.2f, c=%.2f, s=%.2f, turn=%.2f\n",
           name,
           xform.pos.x, xform.pos.y,
           xform.rot.c, xform.rot.s,
           xform.turn);

}

//////////////////////////////////////////////////////////////////////

void buffer_create(buffer_t* buffer, size_t elem_size, size_t capacity) {
    
    size_t alloc_size = elem_size * capacity;
    
    buffer->elem_size = elem_size;
    buffer->count = 0;
    buffer->capacity = capacity;
    buffer->data = malloc(alloc_size);
    
}

void buffer_resize(buffer_t* buffer, size_t new_count) {

    if (new_count > buffer->capacity) {
        
        size_t new_capacity = buffer->capacity;
        
        while (new_capacity <= new_count) {
            new_capacity *= 2;
        }

        size_t alloc_size = buffer->elem_size * new_capacity;

        buffer->data = realloc(buffer->data, alloc_size);
        buffer->capacity = new_capacity;
        
    }

    buffer->count = new_count;
        
}

void buffer_extend(buffer_t* buffer, const void* elements, size_t count) {

    size_t offset = buffer->elem_size * buffer->count;

    buffer_resize(buffer, buffer->count + count);

    memcpy(buffer->data + offset, elements, count*buffer->elem_size);

}

const void* buffer_read(const buffer_t* buffer, size_t idx) {
    return buffer->data + idx*buffer->elem_size;
}
 
void buffer_append(buffer_t* buffer, const void* elem) {
    buffer_extend(buffer, elem, 1);
}

void buffer_pop(buffer_t* buffer, void* dst) {
    buffer->count -= 1;
    size_t offset = buffer->count * buffer->elem_size;
    memcpy(dst, (const void*)buffer->data + offset, buffer->elem_size);
}

void buffer_destroy(buffer_t* buffer) {
    free(buffer->data);
    memset(buffer, 0, sizeof(buffer_t));
}

//////////////////////////////////////////////////////////////////////

void lsystem_create(lsystem_t* lsys,
                    const char* name,
                    const char* start,
                    lsystem_rule_t const rules[],
                    double turn_angle_deg,
                    const char* draw_chars) {

    lsys->name = name;
    lsys->start = start;

    memset(lsys->rules, 0, sizeof(lsys->rules));
    memset(lsys->draw_chars, 0, sizeof(lsys->draw_chars));

    for (const lsystem_rule_t* src_rule=rules; src_rule->variable; ++src_rule) {
        lsystem_rule_tagged_t* dst_rule = lsys->rules + (int)src_rule->variable;
        dst_rule->replacement = src_rule->replacement;
        dst_rule->length = strlen(src_rule->replacement);
    }

    lsys->turn_angle_rad = turn_angle_deg * M_PI / 180.f;

    for (int i=0; i<=MAX_CYCLE_LENGTH; ++i) {
        
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


void lsystem_print(const lsystem_t* lsys) {

    printf("%s:\n", lsys->name);
    printf("  start: %s\n", lsys->start);
    printf("  rules:\n");
    for (int i=0; i<MAX_RULES; ++i) {
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

char* lsystem_build_string(const lsystem_t* lsys,
                           size_t max_depth) {
    
    buffer_t string_buffers[2];
    
    for (int i=0; i<2; ++i) {
        buffer_create(string_buffers + i, sizeof(char), INIT_STRING_CAPACITY);
    }

    int cur_idx = 0;
    
    buffer_extend(string_buffers + cur_idx,
                  lsys->start,
                  strlen(lsys->start));

    for (int i=0; i<max_depth; ++i) {

        int next_idx = 1 - cur_idx;

        buffer_t* src_buffer = string_buffers + cur_idx;
        buffer_t* dst_buffer = string_buffers + next_idx;

        buffer_resize(dst_buffer, 0);
        
        const char* start = (const char*)src_buffer->data;
        const char* end = start + src_buffer->count;

        for (const char* c=start; c!=end; ++c) {

            const lsystem_rule_tagged_t* rule = lsys->rules + (int)*c;
            
            if (rule->replacement) {
                
                buffer_extend(dst_buffer,
                              rule->replacement,
                              rule->length);

            } else {

                buffer_append(dst_buffer, c);

            }

        }

        cur_idx = next_idx;
        
    }

    const char nul = '\0';
    buffer_append(string_buffers + cur_idx, &nul);

    buffer_destroy(string_buffers + (1 - cur_idx));

    return (char*)string_buffers[cur_idx].data;

}

void lsystem_execute_symbol(const lsystem_t* lsys,
                            const char symbol,
                            buffer_t* segments,
                            xform_t* state,
                            buffer_t* xform_stack) {

    if (isalpha(symbol)) {

        if (lsys->draw_chars[0] && !lsys->draw_chars[(int)symbol]) {
            return;
        }
            
        float xnew = state->pos.x + state->rot.c;
        float ynew = state->pos.y + state->rot.s;

        lsystem_segment_t seg = { { state->pos.x, state->pos.y},
                                  { xnew, ynew } };

        buffer_append(segments, &seg);

        state->pos.x = xnew;
        state->pos.y = ynew;

    } else if (symbol == '+' || symbol == '-') {

        if (lsys->rotation_cycle_length) {

            int delta = (symbol == '+') ? 1 : -1;

            int t = state->turn;

            t = positive_mod(t + delta,
                             lsys->rotation_cycle_length);

            const rot2d_t* r = lsys->rotations + t;

            state->turn = t;
            state->rot = *r;
            
        } else {

            float delta = ( (symbol == '+') ?
                            lsys->turn_angle_rad : -lsys->turn_angle_rad );
            
            state->turn += delta;

            state->rot.c = cosf(state->turn);
            state->rot.s = sinf(state->turn);

        }

    } else if (symbol == '[') {

        buffer_append(xform_stack, state);

    } else if (symbol == ']') {

        buffer_pop(xform_stack, state);

    } else {

        fprintf(stderr, "invalid character in string: %c\n", symbol);
        exit(1);

    }

}

void lsystem_segments_r(const lsystem_t* lsys,
                        const char* lstring, 
                        size_t max_depth,
                        buffer_t* segments,
                        xform_t* cur_state,
                        buffer_t* xform_stack,
                        lsystem_memo_set_t* mset) {

    for (const char* psymbol=lstring; *psymbol; ++psymbol) {

        const lsystem_rule_tagged_t* rule = lsys->rules + (int)*psymbol;

        if (max_depth && rule->replacement) {

            lsystem_memo_t* active_memo = 0;

            if (mset && max_depth == mset->memo_depth) {

                lsystem_memo_t* memo = mset->memos[(int)*psymbol];

                if (memo) {

                    // cur_state * init_inverse
                    xform_t update_xform =
                        xform_compose(*cur_state, memo->init_inverse);

                    size_t init_count = segments->count;
                    buffer_resize(segments, init_count + memo->segment_count);

                    const lsystem_segment_t* start =
                        buffer_read(segments, memo->segment_start);

                    const lsystem_segment_t* end = start + memo->segment_count;

                    lsystem_segment_t* dst =
                        (lsystem_segment_t*)(segments->data +
                                             init_count * segments->elem_size);

                    for (const lsystem_segment_t* seg=start; seg != end; ++seg) {

                        lsystem_segment_t newseg = {
                            xform_transform_point(update_xform, seg->p0),
                            xform_transform_point(update_xform, seg->p1)
                        };

                        *dst++ = newseg;

                    }

                    *cur_state = xform_compose(*cur_state,
                                               memo->delta_xform);

                    continue;

                } else {
                    
                    active_memo = malloc(sizeof(lsystem_memo_t));
                    
                    active_memo->segment_start = segments->count;
                    active_memo->segment_count = segments->count;
                    active_memo->init_inverse = xform_inverse(*cur_state);
                    
                    mset->memos[(int)*psymbol] = active_memo;

                }

            }

            lsystem_segments_r(lsys, rule->replacement,
                               max_depth-1,
                               segments, cur_state, xform_stack,
                               mset);

            if (active_memo) {
                
                active_memo->segment_count =
                    segments->count - active_memo->segment_start;
                
                active_memo->delta_xform =
                    xform_compose(active_memo->init_inverse,
                                  *cur_state);

            }

        } else {

            lsystem_execute_symbol(lsys, *psymbol, segments,
                                   cur_state, xform_stack);

        }
        
    }

}

buffer_t* lsystem_segments_recursive(const lsystem_t* lsys,
                                     size_t max_depth,
                                     size_t memo_depth) {

    buffer_t* segments = malloc(sizeof(buffer_t));
    
    buffer_create(segments, sizeof(lsystem_segment_t), INIT_SEGMENTS_CAPACITY);

    buffer_t xform_stack;
    buffer_create(&xform_stack, sizeof(xform_t), INIT_STATES_CAPACITY);

    xform_t cur_state = IDENTITY_XFORM;

    lsystem_memo_set_t mset;
    memset(&mset, 0, sizeof(mset));

    mset.memo_depth = memo_depth;

    lsystem_segments_r(lsys, lsys->start, 
                       max_depth, segments,
                       &cur_state, &xform_stack,
                       memo_depth ? &mset : NULL);

    for (int i=0; i<MAX_RULES; ++i) {
        if (mset.memos[i]) {
            free(mset.memos[i]);
        }
    }

    buffer_destroy(&xform_stack);

    return segments;

}

buffer_t* lsystem_segments_from_string(const lsystem_t* lsys,
                                       const char* lstring) {

    buffer_t* segments = malloc(sizeof(buffer_t));

    buffer_create(segments, sizeof(lsystem_segment_t), INIT_SEGMENTS_CAPACITY);

    buffer_t xform_stack;
    buffer_create(&xform_stack, sizeof(xform_t), INIT_STATES_CAPACITY);

    xform_t cur_state = IDENTITY_XFORM;

    for (const char* psymbol=lstring; *psymbol; ++psymbol) {

        lsystem_execute_symbol(lsys, *psymbol, segments,
                               &cur_state, &xform_stack);

    }

    buffer_destroy(&xform_stack);

    return segments;

}

void initialize_known_lsystems() {

    lsystem_create(known_lsystems + LSYS_SIERPINSKI_TRIANGLE,
                   "sierpinski_triangle", "F-G-G",
                   (lsystem_rule_t[]){
                     { 'F', "F-G+F+G-F" },
                     { 'G', "GG" }, 
                     { 0, 0 }
                   }, 120, NULL);

    lsystem_create(known_lsystems + LSYS_SIERPINSKI_ARROWHEAD,
                   "sierpinski_arrowhead", "A",
                   (lsystem_rule_t[]){
                     { 'A', "B-A-B" }, 
                     { 'B', "A+B+A" }, 
                     { 0, 0 }
                   }, 60, NULL);
    
    lsystem_create(known_lsystems + LSYS_DRAGON_CURVE,
                   "dragon_curve", "FX",
                   (lsystem_rule_t[]){
                     { 'X', "X+YF+" },
                     { 'Y', "-FX-Y" },
                   }, 90, NULL);
    
    lsystem_create(known_lsystems + LSYS_BARNSLEY_FERN,
                   "barnsley_fern", "X",
                   (lsystem_rule_t[]){
                     { 'X', "F+[[X]-X]-F[-FX]+X" },
                     { 'F', "FF" },
                     { 0, 0 }
                   }, 25, NULL);

    lsystem_create(known_lsystems + LSYS_STICKS,
                   "sticks", "X",
                   (lsystem_rule_t[]){
                     { 'X', "F[+X]F[-X]+X" },
                     { 'F', "FF" },
                     { 0, 0 }
                   }, 20, "F");

    lsystem_create(known_lsystems + LSYS_HILBERT,
                   "hilbert", "L",
                   (lsystem_rule_t[]){
                     { 'L', "+RF-LFL-FR+" },
                     { 'R', "-LF+RFR+FL-" },
                     { 0, 0 }
                   }, 90, "F");

    lsystem_create(known_lsystems + LSYS_PENTAPLEXITY,
                   "pentaplexity", "F++F++F++F++F",
                   (lsystem_rule_t[]){
                     { 'F', "F++F++F+++++F-F++F" },
                     { 0, 0 }
                   }, 36, NULL);
                   
}

//////////////////////////////////////////////////////////////////////

double get_time_as_double() {

    struct timespec tp;

    clock_gettime(CLOCK_REALTIME, &tp);

    return (double)tp.tv_sec + (double)tp.tv_nsec * 1e-9;

}

//////////////////////////////////////////////////////////////////////

void parse_options(int argc, char** argv, options_t* opts) {

    int ok = 1;

    int allow_precomputed_rotation = 1;

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
                    if (!strcmp(arg, known_lsystems[j].name)) {
                        ok = 1;
                        opts->lsys = known_lsystems + j;
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
            
            opts->method = METHOD_STRING;
            
        } else if (!strcmp(arg, "-r")) {

            opts->method = METHOD_RECURSION;
            
        } else if (!strcmp(arg, "-m")) {
            
            if (++i == argc) { ok = 0; break; }
            
            int d;
            
            if (sscanf(argv[i], "%d", &d) != 1) {
                ok = 0; break;
            }
            
            if (d >= 0) {
                opts->memo_depth = d;
            } else {
                ok = 0;
                break;
            }

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
            
        } else if (!strcmp(arg, "-P")) {
            
            allow_precomputed_rotation = 0;

        } else {

            fprintf(stderr, "error: unrecognized option %s\n\n", arg);

            ok = 0;
            break;

        }

    }

    if (opts->memo_depth && opts->method == METHOD_STRING) {
        printf("warning: ignoring memo depth for string method!\n");
    }
    
    if (!ok || !opts->lsys || !opts->max_depth) {
        printf("usage: %s [options] LSYSTEM MAXDEPTH\n"
               "\n"
               "where LSYSTEM is one of:\n", argv[0]);
        for (int j=0; j<NUM_KNOWN_LSYSTEMS; ++j) {
            printf("  * %s\n", known_lsystems[j].name);
        }
        printf("\n");
        printf("options:\n");
        printf("  -s             use string building method\n"
               "  -r             use recursive method\n"
               "  -m MEMODEPTH   enable memoization for recursive method\n"
               "  -x MAXSEGMENTS maximum number of segments for output\n"
               "  -P             don't precompute rotations\n"
               "\n");
        exit(1);
    }

    if (!allow_precomputed_rotation) {
        printf("disabling precomputed rotation!\n");
        opts->lsys->rotation_cycle_length = 0;
    }

}

//////////////////////////////////////////////////////////////////////

int main(int argc, char** argv) {

    initialize_known_lsystems();

    options_t opts;

    parse_options(argc, argv, &opts);

    lsystem_print(opts.lsys);

    printf("using %s method\n",
           opts.method == METHOD_STRING ? "string" : "recursion");

    if (opts.method == METHOD_RECURSION && opts.memo_depth) {
        printf("memo depth is %d\n", (int)opts.memo_depth);
    }

    double start = get_time_as_double();

    buffer_t* segments;

    if (opts.method == METHOD_STRING) {

        char* lstring = lsystem_build_string(opts.lsys,
                                             opts.max_depth);

        segments = lsystem_segments_from_string(opts.lsys, lstring);

        free(lstring);

    } else {

        segments = lsystem_segments_recursive(opts.lsys, opts.max_depth,
                                              opts.memo_depth);

    }

    double elapsed = get_time_as_double() - start;
    printf("got %d segments in %.6f sec.\n", (int)segments->count, elapsed);

    if (segments->count > opts.max_segments) {
        
        printf("...maximum of %d segments exceeded, skipping output!\n",
               (int)opts.max_segments);

    } else {

        const lsystem_segment_t* segment = (const lsystem_segment_t*)segments->data;
        const lsystem_segment_t* end = segment + segments->count;

        FILE* outfile = fopen("segments.txt", "w");
        
        for ( ; segment != end; ++segment) {
            fprintf(outfile, "%g %g %g %g\n",
                    segment->p0.x, segment->p0.y,
                    segment->p1.x, segment->p1.y);
        }

        fclose(outfile);

        printf("wrote segments.txt\n");

    }

    buffer_destroy(segments);
    free(segments);

    return 0;

}
