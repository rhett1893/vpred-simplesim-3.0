#ifndef _stride_h_
#define _stride_h_

#include "machine.h"
#include "memory.h"
/* constants for entry lookup, update status in stride table: */
/* MISS - not found */
#define MISS	 	-1

/* HIT - found and correct prediction */
#define HIT 	 	1

/* HITFAULT - found and uncorrect prediction */
#define HITFAULT  	2


/* op data type table constants */
#define SINGLE	 	0
#define DOUBLE		1

/* pc stride (in bits) */
#define INST_OFFSET	0

#define is_PRED   ((MD_OP_FLAGS(op) & F_LOAD))
/* Single precision table entry definition */
typedef struct ln_i1{
  unsigned long	pc;     /* entry PC (tag) */
  int last_val;         /* last output value */
  int stride;           /* the stride between the two known last values */
  int last_alloc;       /* allocation sim_cycle time */
  int last_ref;         /* last lookup sim_cycle time */
  char ps;              /* classification fsm present state:
			             0,1 - don't use prediction (don't go)
			             1,2 - use prediction (go)           */
  char valid;           /* valid bit */
  char x;               /* X machanism - for high bandwidth fetch:
			   when using X mechanism predicted value is not:
			   last_val + stride, but: last_val + X*stride.
			   X is incremented with every fetch (lookup) and
			   decremented when the instruction leaves (either
			   by normal terminition (update) or when instruction
			   doesn't exit via commit, i.e. pipe is flushed
			   (lookup_undo) */
  int accessed;
  int pred_correct;
} Hash_i1, *Hash_ptr_i1;

/* Double precision table entry definition */
typedef struct ln_i2{
  unsigned long	pc;
  int last_val[2];
  int stride[2];
  int last_alloc;
  int last_ref;
  char ps;
  char valid;
  char x;
} Hash_i2, *Hash_ptr_i2;

typedef struct value_t_t
{
  int single_p;
  int double_p[2];
}value_t;

typedef struct VAL_TAG_TYPE_t
{
  int8_t fsm_pred;
  value_t value;
}VAL_TAG_TYPE;



extern int                       en_lookup_stride;
// extern tick_t                    sim_cycle;
extern int                       use_stride;
extern int                       use_fsm;
extern int                       use_trace_cache;
extern int                       vp_replace;
extern int                       start_fsm;
extern int                       hash_no_i1;
extern int                       hash_no_i2;
extern int                       hash_asso_i1;
extern int                       hash_asso_i2;


/* function headers : */

void create_stride();

int update(md_addr_t instpc,md_inst_t inst,VAL_TAG_TYPE pred_val,VAL_TAG_TYPE calc_val,struct mem_t *mem, int my_ref);

int lookup_undo(md_addr_t instpc, md_inst_t inst, VAL_TAG_TYPE pred_val);

void allocate(md_addr_t pred_PC,md_inst_t inst,VAL_TAG_TYPE calc_val,int my_ref);

void lookup(md_addr_t pred_PC,md_inst_t inst,VAL_TAG_TYPE *pred_val,struct mem_t *mem);

#endif
