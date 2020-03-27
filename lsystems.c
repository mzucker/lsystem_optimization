#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <ctype.h>
#include <math.h>
#include <time.h>
#include <unistd.h>

#include <assert.h>
#ifdef NDEBUG
#undef NDEBUG
#endif

#ifdef _OPENMP
#include <omp.h>
#endif

#include <pthread.h>

typedef struct dynarray {

    size_t elem_size;
    size_t capacity;
    size_t count;
    
    unsigned char* data;

} dynarray_t;

//////////////////////////////////////////////////////////////////////

#ifndef _OPENMP

typedef void (thread_pool_func_t)(void*);

typedef struct thread_pool_task {

    thread_pool_func_t* func;
    void* data;
    
} thread_pool_task_t;

enum {
    INIT_TASKS_CAPACITY = 16
};

typedef struct thread_pool {

    size_t thread_count;
    pthread_t* threads;

    dynarray_t tasks;
    size_t next_task_idx;
    size_t finished_count;

    pthread_mutex_t lock;
    pthread_cond_t start_cond;
    pthread_cond_t end_cond;

} thread_pool_t;

#endif

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
    size_t min_parallel_segments;

#ifndef _OPENMP    
    size_t num_tasks;
    thread_pool_t pool;
#endif
    
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

lsystem_t KNOWN_LSYSTEMS[NUM_KNOWN_LSYSTEMS];

//////////////////////////////////////////////////////////////////////

typedef enum method {
    METHOD_STRING,
    METHOD_RECURSION
} method_t;

typedef struct options {
    lsystem_t* lsys;
    size_t     max_depth;
    size_t     memo_depth;
    size_t     min_parallel_segments;
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

    return (vec2d_t){
        R.c * p.x - R.s * p.y,
        R.s * p.x + R.c * p.y
    };

}

rot2d_t rotateR(const rot2d_t R2, const rot2d_t R1) {

    return (rot2d_t) {
        R2.c * R1.c - R2.s * R1.s,
        R2.s * R1.c + R2.c * R1.s,
    };

}

vec2d_t translate(const vec2d_t p, const vec2d_t q) {
    
    return (vec2d_t) { p.x + q.x, p.y + q.y };

}

xform_t xform_inverse(xform_t xform) {

    rot2d_t rinv = { xform.rot.c, -xform.rot.s };

    return (xform_t) {
        .pos = rotate(rinv, (vec2d_t){ -xform.pos.x, -xform.pos.y }),
        .rot = rinv,
        .turn = -xform.turn
    };

}

xform_t xform_compose(xform_t xform2,
                      xform_t xform1) {

    return (xform_t) {

        .pos = translate(xform2.pos,
                         rotate(xform2.rot, xform1.pos)),

        .rot = rotateR(xform2.rot, xform1.rot),

        .turn = xform2.turn + xform1.turn

    };
    
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

void dynarray_create(dynarray_t* dynarray, size_t elem_size, size_t capacity) {
    
    size_t alloc_size = elem_size * capacity;
    
    dynarray->elem_size = elem_size;
    dynarray->count = 0;
    dynarray->capacity = capacity;
    dynarray->data = malloc(alloc_size);
    
}

void dynarray_resize(dynarray_t* dynarray, size_t new_count) {

    if (new_count > dynarray->capacity) {
        
        size_t new_capacity = dynarray->capacity;
        
        while (new_capacity <= new_count) {
            new_capacity *= 2;
        }

        size_t alloc_size = dynarray->elem_size * new_capacity;

        dynarray->data = realloc(dynarray->data, alloc_size);
        dynarray->capacity = new_capacity;
        
    }

    dynarray->count = new_count;
        
}

void dynarray_extend(dynarray_t* dynarray, const void* elements, size_t count) {

    size_t offset = dynarray->elem_size * dynarray->count;

    dynarray_resize(dynarray, dynarray->count + count);

    memcpy(dynarray->data + offset, elements, count*dynarray->elem_size);

}

void* dynarray_elem_ptr(dynarray_t* dynarray, size_t idx) {
    return dynarray->data + idx*dynarray->elem_size;
}

const void* dynarray_const_elem_ptr(const dynarray_t* dynarray, size_t idx) {
    return dynarray->data + idx*dynarray->elem_size;
}

void dynarray_get(const dynarray_t* dynarray, size_t idx, void* dst) {
    memcpy(dst, dynarray_const_elem_ptr(dynarray, idx), dynarray->elem_size);
}

void dynarray_set(dynarray_t* dynarray, size_t idx, const void* src) {
    memcpy(dynarray_elem_ptr(dynarray, idx), src, dynarray->elem_size);
}
 
void dynarray_append(dynarray_t* dynarray, const void* elem) {
    dynarray_extend(dynarray, elem, 1);
}

void dynarray_pop(dynarray_t* dynarray, void* dst) {
    dynarray->count -= 1;
    size_t offset = dynarray->count * dynarray->elem_size;
    memcpy(dst, (const void*)dynarray->data + offset, dynarray->elem_size);
}

void dynarray_clear(dynarray_t* dynarray) {
    dynarray->count = 0;
}

void dynarray_destroy(dynarray_t* dynarray) {
    free(dynarray->data);
    memset(dynarray, 0, sizeof(dynarray_t));
}


//////////////////////////////////////////////////////////////////////

#ifndef _OPENMP

void* thread_pool_start(void* p) {
    
    thread_pool_t* pool = p;

    while (1) {

        pthread_mutex_lock(&pool->lock);

        while (pool->next_task_idx == pool->tasks.count) {
            ++pool->finished_count;
            pthread_cond_broadcast(&pool->end_cond);
            //printf("finished! count is now %d\n", (int)pool->finished_count);
            pthread_cond_wait(&pool->start_cond, &pool->lock);
        }

        thread_pool_task_t task;
        dynarray_get(&pool->tasks, pool->next_task_idx, &task);
        ++pool->next_task_idx;

        //printf("got task %d/%d with func %p\n",
        //(int)pool->next_task_idx, (int)pool->tasks.count, task.func);

        pthread_mutex_unlock(&pool->lock);
        //pthread_yield_np();

        if (!task.func) {
            return 0;
        }

        task.func(task.data);

    }

    return 0;

}

void thread_pool_create(thread_pool_t* pool, size_t nthreads) {

    memset(pool, 0, sizeof(thread_pool_t));

    dynarray_create(&pool->tasks, sizeof(thread_pool_task_t), INIT_TASKS_CAPACITY);

    if (nthreads == 1) { return; }
    
    const pthread_mutexattr_t* m_attr = NULL;
    pthread_mutex_init(&pool->lock, m_attr);

    const pthread_condattr_t* c_attr = NULL;
    pthread_cond_init(&pool->start_cond, c_attr);
    pthread_cond_init(&pool->end_cond, c_attr);

    pool->thread_count = nthreads;
    pool->threads = (pthread_t*) malloc(nthreads * sizeof(pthread_t));

    pthread_mutex_lock(&pool->lock);

    for (size_t i=0; i<nthreads; ++i) {
        const pthread_attr_t* t_attr = NULL;
        pthread_create(pool->threads + i, t_attr, thread_pool_start, pool);
    }

    //printf("waiting for threads to start up...\n");

    while (pool->finished_count < nthreads) {
        pthread_cond_wait(&pool->end_cond, &pool->lock);
    }

    pool->finished_count = 0;

    //printf("all threads started!\n");
    pthread_mutex_unlock(&pool->lock);

}

void thread_pool_add_task(thread_pool_t* pool, thread_pool_func_t* func, void* data) {

    //printf("added a task with func %p\n", func);
    
    thread_pool_task_t task = { func, data };
    dynarray_append(&pool->tasks, &task);
    
}


void thread_pool_run_in_main_thread(thread_pool_t* pool) {

    for (size_t i=0; i<pool->tasks.count; ++i) {
        thread_pool_task_t task;
        dynarray_get(&pool->tasks, i, &task);
        task.func(task.data);
    }

    dynarray_clear(&pool->tasks);

}

void thread_pool_run(thread_pool_t* pool) {

    if (!pool->thread_count) {
        
        thread_pool_run_in_main_thread(pool);
        
    } else {

        pool->finished_count = 0;
        pool->next_task_idx = 0;
        
        pthread_mutex_lock(&pool->lock);
        pthread_cond_broadcast(&pool->start_cond);

        while (pool->finished_count < pool->thread_count) {
            pthread_cond_wait(&pool->end_cond, &pool->lock);
        }

        pthread_mutex_unlock(&pool->lock);
        pool->next_task_idx = 0;
        dynarray_clear(&pool->tasks);

    }
    
}

void thread_pool_destroy(thread_pool_t* pool) {

    if (pool->thread_count) {
        
        for (size_t i=0; i<pool->thread_count; ++i) {
            thread_pool_add_task(pool, NULL, NULL);
        }

        pthread_cond_broadcast(&pool->start_cond);
        
        for (size_t i=0; i<pool->thread_count; ++i) {
            pthread_join(pool->threads[i], NULL);
        }
        
        free(pool->threads);

        pthread_mutex_destroy(&pool->lock);
        pthread_cond_destroy(&pool->start_cond);
        pthread_cond_destroy(&pool->end_cond);
        
    }

    dynarray_destroy(&pool->tasks);
    
    
}

#endif

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

void initialize_known_lsystems() {

    lsystem_create(KNOWN_LSYSTEMS + LSYS_SIERPINSKI_TRIANGLE,
                   "sierpinski_triangle", "F-G-G",
                   (lsystem_rule_t[]){
                     { 'F', "F-G+F+G-F" },
                     { 'G', "GG" }, 
                     { 0, 0 }
                   }, 120, NULL);

    lsystem_create(KNOWN_LSYSTEMS + LSYS_SIERPINSKI_ARROWHEAD,
                   "sierpinski_arrowhead", "A",
                   (lsystem_rule_t[]){
                     { 'A', "B-A-B" }, 
                     { 'B', "A+B+A" }, 
                     { 0, 0 }
                   }, 60, NULL);
    
    lsystem_create(KNOWN_LSYSTEMS + LSYS_DRAGON_CURVE,
                   "dragon_curve", "FX",
                   (lsystem_rule_t[]){
                     { 'X', "X+YF+" },
                     { 'Y', "-FX-Y" },
                   }, 90, NULL);
    
    lsystem_create(KNOWN_LSYSTEMS + LSYS_BARNSLEY_FERN,
                   "barnsley_fern", "X",
                   (lsystem_rule_t[]){
                     { 'X', "F+[[X]-X]-F[-FX]+X" },
                     { 'F', "FF" },
                     { 0, 0 }
                   }, 25, NULL);

    lsystem_create(KNOWN_LSYSTEMS + LSYS_STICKS,
                   "sticks", "X",
                   (lsystem_rule_t[]){
                     { 'X', "F[+X]F[-X]+X" },
                     { 'F', "FF" },
                     { 0, 0 }
                   }, 20, "F");

    lsystem_create(KNOWN_LSYSTEMS + LSYS_HILBERT,
                   "hilbert", "L",
                   (lsystem_rule_t[]){
                     { 'L', "+RF-LFL-FR+" },
                     { 'R', "-LF+RFR+FL-" },
                     { 0, 0 }
                   }, 90, "F");

    lsystem_create(KNOWN_LSYSTEMS + LSYS_PENTAPLEXITY,
                   "pentaplexity", "F++F++F++F++F",
                   (lsystem_rule_t[]){
                     { 'F', "F++F++F+++++F-F++F" },
                     { 0, 0 }
                   }, 36, NULL);
                   
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
    
    dynarray_t string_dynarrays[2];
    
    for (int i=0; i<2; ++i) {
        dynarray_create(string_dynarrays + i, sizeof(char), INIT_STRING_CAPACITY);
    }

    int cur_idx = 0;
    
    dynarray_extend(string_dynarrays + cur_idx,
                  lsys->start,
                  strlen(lsys->start));

    for (int i=0; i<max_depth; ++i) {

        int next_idx = 1 - cur_idx;

        dynarray_t* src_dynarray = string_dynarrays + cur_idx;
        dynarray_t* dst_dynarray = string_dynarrays + next_idx;

        dynarray_resize(dst_dynarray, 0);
        
        const char* start = (const char*)src_dynarray->data;
        const char* end = start + src_dynarray->count;

        for (const char* c=start; c!=end; ++c) {

            const lsystem_rule_tagged_t* rule = lsys->rules + (int)*c;
            
            if (rule->replacement) {
                
                dynarray_extend(dst_dynarray,
                              rule->replacement,
                              rule->length);

            } else {

                dynarray_append(dst_dynarray, c);

            }

        }

        cur_idx = next_idx;
        
    }

    const char nul = '\0';
    dynarray_append(string_dynarrays + cur_idx, &nul);

    dynarray_destroy(string_dynarrays + (1 - cur_idx));

    return (char*)string_dynarrays[cur_idx].data;

}

void lsystem_execute_symbol(const lsystem_t* lsys,
                            const char symbol,
                            dynarray_t* segments,
                            xform_t* state,
                            dynarray_t* xform_stack) {

    if (isalpha(symbol)) {

        if (lsys->draw_chars[0] && !lsys->draw_chars[(int)symbol]) {
            return;
        }
            
        float xnew = state->pos.x + state->rot.c;
        float ynew = state->pos.y + state->rot.s;

        lsystem_segment_t seg = { { state->pos.x, state->pos.y},
                                  { xnew, ynew } };

        dynarray_append(segments, &seg);

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

        dynarray_append(xform_stack, state);

    } else if (symbol == ']') {

        dynarray_pop(xform_stack, state);

    } else {

        fprintf(stderr, "invalid character in string: %c\n", symbol);
        exit(1);

    }

}

typedef struct segment_xform_data {
    
    lsystem_segment_t* dst;
    const lsystem_segment_t* src;
    size_t count;
    const xform_t* update_xform;
       
} segment_xform_data_t;

void lsystem_transform_segments(void* ptr) {

    const segment_xform_data_t* sdata = ptr;

    lsystem_segment_t* dst = sdata->dst;
    const lsystem_segment_t* src = sdata->src;
                    
    for (size_t i=0; i<sdata->count; ++i) {

        lsystem_segment_t newsrc = {
            xform_transform_point(*sdata->update_xform, src->p0),
            xform_transform_point(*sdata->update_xform, src->p1)
        };

        *dst++ = newsrc;
        src++;

    }

}


void lsystem_segments_r(const lsystem_t* lsys,
                        const char* lstring, 
                        size_t max_depth,
                        dynarray_t* segments,
                        xform_t* cur_state,
                        dynarray_t* xform_stack,
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
                    dynarray_resize(segments, init_count + memo->segment_count);

                    int do_parallelize =
                        (memo->segment_count >= mset->min_parallel_segments);
                    
                    const lsystem_segment_t* src =
                        dynarray_const_elem_ptr(segments, memo->segment_start);

                    lsystem_segment_t* dst =
                        dynarray_elem_ptr(segments, init_count);

#ifdef _OPENMP

                    #pragma omp parallel for if (do_parallelize)
                    for (size_t i=0; i<memo->segment_count; ++i) {
                        lsystem_segment_t newsrc = {
                            xform_transform_point(update_xform, src[i].p0),
                            xform_transform_point(update_xform, src[i].p1)
                        };
                        dst[i] = newsrc;
                    }
                    
#else
                    
                    if (do_parallelize) {
                        
                        size_t batch_size =
                            (int)ceil((double)memo->segment_count / mset->num_tasks);

                        segment_xform_data_t sdata[mset->num_tasks];

                        size_t start_idx = 0;

                        for (size_t i=0; i<mset->num_tasks; ++i) {

                            size_t end_idx = start_idx + batch_size;
                            
                            if (end_idx > memo->segment_count) {
                                end_idx = memo->segment_count;
                            }
                            
                            segment_xform_data_t* si = sdata + i;
                            
                            si->dst = dst + start_idx;
                            si->src = src + start_idx;
                            si->count = end_idx - start_idx;
                            si->update_xform = &update_xform;

                            start_idx = end_idx;

                            thread_pool_add_task(&mset->pool,
                                                 lsystem_transform_segments,
                                                 si);

                        }

                        thread_pool_run(&mset->pool);
                        
                    } else {

                        segment_xform_data_t sdata = {
                            dst, src, memo->segment_count, &update_xform
                        };

                        lsystem_transform_segments(&sdata);

                    }

#endif                        
                    *cur_state = xform_compose(*cur_state, memo->delta_xform);

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

dynarray_t* lsystem_segments_recursive(const lsystem_t* lsys,
                                       size_t max_depth,
                                       size_t memo_depth,
                                       size_t min_parallel_segments) {

    dynarray_t* segments = malloc(sizeof(dynarray_t));
    
    dynarray_create(segments, sizeof(lsystem_segment_t), INIT_SEGMENTS_CAPACITY);

    dynarray_t xform_stack;
    dynarray_create(&xform_stack, sizeof(xform_t), INIT_STATES_CAPACITY);

    xform_t cur_state = IDENTITY_XFORM;

    lsystem_memo_set_t mset;
    memset(&mset, 0, sizeof(mset));

    mset.memo_depth = memo_depth;
    mset.min_parallel_segments = min_parallel_segments;


    if (min_parallel_segments != (size_t)-1) {
#ifndef _OPENMP    
        mset.num_tasks = sysconf(_SC_NPROCESSORS_ONLN);
        thread_pool_create(&mset.pool, mset.num_tasks);
        printf("created a thread pool with %d threads!\n", (int)mset.num_tasks);
#else
        printf("will parallelize with OpenMP!\n");
#endif
    }

    lsystem_segments_r(lsys, lsys->start, 
                       max_depth, segments,
                       &cur_state, &xform_stack,
                       memo_depth ? &mset : NULL);

    for (int i=0; i<MAX_RULES; ++i) {
        if (mset.memos[i]) {
            free(mset.memos[i]);
        }
    }

#ifndef _OPENMP    
    thread_pool_destroy(&mset.pool);
#endif    

    dynarray_destroy(&xform_stack);

    return segments;

}

dynarray_t* lsystem_segments_from_string(const lsystem_t* lsys,
                                       const char* lstring) {

    dynarray_t* segments = malloc(sizeof(dynarray_t));

    dynarray_create(segments, sizeof(lsystem_segment_t), INIT_SEGMENTS_CAPACITY);

    dynarray_t xform_stack;
    dynarray_create(&xform_stack, sizeof(xform_t), INIT_STATES_CAPACITY);

    xform_t cur_state = IDENTITY_XFORM;

    for (const char* psymbol=lstring; *psymbol; ++psymbol) {

        lsystem_execute_symbol(lsys, *psymbol, segments,
                               &cur_state, &xform_stack);

    }

    dynarray_destroy(&xform_stack);

    return segments;

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
    opts->min_parallel_segments = -1;

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

        } else if (!strcmp(arg, "-t")) {
            
            if (++i == argc) { ok = 0; break; }
            
            int d;
            
            if (sscanf(argv[i], "%d", &d) != 1) {
                ok = 0; break;
            }
            
            if (d >= -1) {
                opts->min_parallel_segments = (size_t)d;
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
        opts->memo_depth = 0;
    }

#ifdef _OPENMP    
    if (opts->min_parallel_segments != (size_t)-1 && !opts->memo_depth) {
        printf("warning: ignoring min parallel segments without memoization\n");
        opts->min_parallel_segments = -1;
    }
    if (opts->min_parallel_segments != (size_t)-1) {
        printf("will parallelize >= %d segments with %d processors\n",
               (int)opts->min_parallel_segments,
               (int)omp_get_num_procs());
    }
#endif    
    
    if (!ok || !opts->lsys || !opts->max_depth) {
        printf("usage: %s [options] LSYSTEM MAXDEPTH\n"
               "\n"
               "where LSYSTEM is one of:\n", argv[0]);
        for (int j=0; j<NUM_KNOWN_LSYSTEMS; ++j) {
            printf("  * %s\n", KNOWN_LSYSTEMS[j].name);
        }
        printf("\n");
        printf("options:\n");
        printf("  -s             use string building method\n"
               "  -r             use recursive method\n"
               "  -m MEMODEPTH   enable memoization for recursive method\n"
               "  -t SEGMENTS    parallelize memoization for more than SEGMENTS segments\n"
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

    dynarray_t* segments;

    if (opts.method == METHOD_STRING) {

        char* lstring = lsystem_build_string(opts.lsys,
                                             opts.max_depth);

        segments = lsystem_segments_from_string(opts.lsys, lstring);

        free(lstring);

    } else {

        segments = lsystem_segments_recursive(opts.lsys, opts.max_depth,
                                              opts.memo_depth,
                                              opts.min_parallel_segments);

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

    dynarray_destroy(segments);
    free(segments);

    return 0;

}
