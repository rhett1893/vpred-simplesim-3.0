/*
 * dpred.h - data value predictor interfaces
 *
 * This file is a part of the SimpleScalar tool suite written by
 * Todd M. Austin as a part of the Multiscalar Research Project.
 * The tool suite is currently maintained by Doug Burger and Todd M. Austin.
 * Copyright (C) 1994, 1995, 1996, 1997 by Todd M. Austin
 * INTERNET: dburger@cs.wisc.edu
 * US Mail:  1210 W. Dayton Street, Madison, WI 53706
 *
 * This routine is written by Sang-Jeong Lee, 
 * Soonchunhyang University, Asan, Chungnam, Korea.
 * Copyright (C) 1998 by Sang-Jeong Lee
 * INTERNET: sjlee@asan.sch.ac.kr
 *
 */
             

#ifndef DPRED_H

#define DPRED_H
#define dassert(a) assert(a)
/* define the predicted instruction type */
#define is_PRED 	( ((SS_OP_FLAGS(op) & F_ICOMP) || 	\
                         (SS_OP_FLAGS(op) & F_LOAD))  &&	\
                         !(SS_OP_FLAGS(op) & F_LONGLAT) )

#include <stdio.h>
#include "misc.h"
#include "ss.h"
#include "stats.h"


/*
 * This module implements a number of data value predictor mechanisms.  The
 * following predictors are supported:
 *
 *	DPred2Level:  two level hybrid data value predictor
 *		K.Wang and M.Franklin, "Highly Accurate Data Prediction
 *                 using Hybrid Predictors", Micro-30, 1997.
 *	
 *  	DPredStride:  strid-based data value predictor	
 *	
 *	DPredLast:  last data value predictor
 *   		M.H.Llipasti and J.P.Shen, "Exceeding the Limit via Value
 *		   Prediction", Micro-29, 1996.	
 *
 *
 */

/* data value predictor types */
enum dpred_class {
  DPredHybrid,                  /* hybrid data pred */
  DPred2Level,                  /* 2-level data pred */
  DPredStride,			/* stride-based data pred */
  DPredLast,			/* last data pred */
  DPred_NUM
};

/* states for stride predictor */
enum stride_state {
  Init, 
  Transient,
  Steady
};

/* an entry in a cache table for each data pred */
struct dpred_ent {
  enum dpred_class class;	/* type of predictor */
  SS_ADDR_TYPE addr;		/* address of an instruction being predicted */
  enum ss_opcode op;		/* opcode of branch corresp. to addr */
  
  union {
    union {
      unsigned char ct_counter;	/* prediction state counter in CT */
      int *vpt_data;		/* data values in VPT */
    } last;
    struct {
      enum stride_state state;
      int value;
      int stride;
    } stride;  
    union {
      struct {
      	unsigned char *lru_info;  /* descending from [0] to [n] by the degree 
                                     of least recently used */
      	int *values;		/* data values */
        unsigned int vhp;		/* value history pattern */
      } vht_val;
      unsigned char *pht_val;	/* pht counter */
     } two;
    union {
      struct {
      	enum stride_state state;   /* stride state */
      	int stride;		   /* stride */
      	unsigned char *lru_info;  /* descending from [0] to [n] by the degree 
                                     of least recently used */
      	int *values;		/* data values */
        unsigned int vhp;		/* value history pattern */
      } vht_val;
      unsigned char *pht_val;	/* pht counter */
     } hybrid;
  } ent;
  struct dpred_ent *prev, *next; /* lru chaining pointers */
};  

/* cache table def */
struct dpred_cache {
  int sets;		/* num cache sets */
  int assoc;		/* cache associativity */
  struct dpred_ent *table;	/* prediction state table */
}; 

/* data predictor def */
struct dpred_data {
  enum dpred_class class;	/* type of predictor */
  union {
    struct {
      unsigned int ct_size;	/* number of entries in CT table */
      unsigned int counter_size;/* bits number of a counter in CT */
      unsigned int vpt_size;	/* number of entries in VPT table */
      unsigned int hist;	/* number of history depth in VPT table */
      struct dpred_cache *ct;	/* classification table */
      struct dpred_cache *vpt;	/* value prediction table */
    } last;
    struct {
      unsigned int vht_size;	/* number of entries in VHT(value history table) */
      struct dpred_cache *vht;  /* Value History Table */
    } stride;
    struct {
      unsigned int vht_size;    /* number of entries in VHT table */
      unsigned int threshold;   /* threshold values for each counter in PHT */
      unsigned int pht_size;	/* number of entries in PHT table */
      unsigned int hist;	/* number of history depth in VHT table */     
      unsigned int xor;		/* VHP xor flag */
      struct dpred_cache *vht;	/* value history table */
      struct dpred_cache *pht;	/* pattern history table */
    } two;
    struct {
      unsigned int vht_size;    /* number of entries in VHT table */
      unsigned int threshold;   /* threshold values for each counter in PHT */
      unsigned int pht_size;	/* number of entries in PHT table */
      unsigned int hist;	/* number of history depth in VHT table */     
      unsigned int xor;		/* VHP xor flag */
      struct dpred_cache *vht;	/* value history table */
      struct dpred_cache *pht;	/* pattern history table */
      char is_stride;		/* flag for tarce info. */
    } hybrid;
  } config;
};

/* data predictor def */
struct dpred {
  enum dpred_class class;	/* type of predictor */
  int sh_trace;		/* flag for show instruction trace */ 
  struct {
    struct dpred_data *last;     /* last data predictor */
    struct dpred_data *stride;   /* stride data predictor */
    struct dpred_data *two;	 /* two level data predictor */
    struct dpred_data *hybrid;	 /* hybrid data predictor */
  } datapred;

  /* stats */
  SS_COUNTER_TYPE lookups;	/* num lookups */

  SS_COUNTER_TYPE data_hits;	/* num correct data-predictions */
  SS_COUNTER_TYPE misses;	/* num incorrect predictions */

  SS_COUNTER_TYPE no_hits;	/* num correct no predictions */
  SS_COUNTER_TYPE no_misses;	/* num incorrect no predictions */

  SS_COUNTER_TYPE l1_misses;	/* num l1 table cache miss */
  SS_COUNTER_TYPE l2_misses;	/* num l2 table cache miss */

  SS_COUNTER_TYPE alias;	/* num aliasing in PHT */
  SS_COUNTER_TYPE alias_hits;	/* correct even if aliasing */
  SS_COUNTER_TYPE alias_misses;	/* incorrect and aliasing */
};


/* create a data predictor */
struct dpred *				/* data predictory instance */
dpred_create(enum dpred_class class,	/* type of predictor to create */
      	     unsigned int l1_size,	/* number of entries in l1 table */
	     unsigned int counter,      /* bits number of a counter 
					   or threshold in two level */
	     unsigned int l2_size,	/* number of entries in l2 table */
	     unsigned int hist,         /* number of history depth  */
	     unsigned int xor);		/* VHP xor flag in two level */

/* create a data predictor table */
struct dpred_data *		/* data predictor table instance */
dpred_data_create (enum dpred_class class,	/* type of predictor to create */
      	     	   unsigned int l1_size,	/* number of entries in l1 table */
	           unsigned int counter,        /* bits number of a counter
					           or threshold in two level */
	           unsigned int l2_size,	/* number of entries in l2 table */
	           unsigned int hist,	        /* number of history depth */
	            unsigned int xor);		/* VHP xor flag in two level */

/* print data predictor configuration */
void
dpred_data_config(
  struct dpred_data *pred_data, /* data predictor instance */
  char name[],                  /* predictor name */
  FILE *stream);        	/* output stream */

/* print predictor configuration */
void
dpred_config(struct dpred *pred,	/* data predictor instance */
	     FILE *stream);		/* output stream */

/* print predictor stats */
void
dpred_stats(struct dpred *pred,		/* data predictor instance */
	    FILE *stream);		/* output stream */

/* register data predictor stats */
void
dpred_reg_stats(struct dpred *pred,	/* data predictor instance */
		struct stat_sdb_t *sdb);/* stats database */

/* reset stats after priming, if appropriate */
void dpred_after_priming(struct dpred *dpred);

/* lookup cache table */
struct dpred_ent *				/* pointer to table entry */
dpred_cache_lookup (
  struct dpred_cache *pred_cache,/* cache table instance */
  SS_ADDR_TYPE daddr);            /* PC address */

/* probe a predictor for a data value, the predictor is probed
   with instruction PC address */
SS_WORD_TYPE				/* predicted data value */
dpred_lookup(struct dpred *pred,	/* data predictor instance */
	     SS_ADDR_TYPE pcaddr,	/* instruction PC address */
	     enum ss_opcode,	        /* opcode of instruction */  
	     char *no_pred,		/* predictor state pointer */
	     struct dpred_ent **ptbl1,	/* pointer of table 1 entry */
	     struct dpred_ent **ptbl2);	/* pointer of table 2 entry */

/* lookup cache table and update LRU state */
struct dpred_ent *				/* pointer to table entry */
dpred_cache_lru_update (
  struct dpred_cache *pred_cache,/* cache table instance */
  SS_ADDR_TYPE daddr);            /* PC address */

/* update the data predictor.  Predictor statistics are updated 
   with result of prediction, indicated by CORRECT. 
   Predictor state to be updated is indicated by *DATA_UPDATE_PTR 
 */
void
dpred_update(struct dpred *pred,	/* data predictor instance */
	     SS_ADDR_TYPE daddr,	/* PC address */
	     SS_WORD_TYPE data,		/* resolved data */
	     char no_pred,		/* no pred flag */
	     int correct,		/* was earlier data prediction ok? */
	     enum ss_opcode op);	/* opcode of instruction */

/* print the trace of the table content */
void				
dpred_trace(struct dpred *pred,	/* data predictor instance */
	     SS_ADDR_TYPE pcaddr,	/* instruction PC address */
	     SS_WORD_TYPE pdata,	/* predicted data */
	     SS_WORD_TYPE rdata,	/* resolved data */
	     char *no_pred,		/* predictor state pointer */
	     struct dpred_ent *ptbl1,	/* pointer of table 1 entry */
	     struct dpred_ent *ptbl2);	/* pointer of table 2 entry */
#endif /* DPRED_H */
