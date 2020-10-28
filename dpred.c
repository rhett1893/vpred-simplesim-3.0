/*
 * dpred.c - data predictor routines
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
            
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>

#include "misc.h"
#include "ss.h"
#include "dpred.h"

/* create a data predictor */
struct dpred *				/* data predictory instance */
dpred_create(enum dpred_class class,	/* type of predictor to create */
      	     unsigned int l1_size,	/* number of entries in l1 table */
	     unsigned int counter,      /* bits number of a counter or
					   threshold of two level */
	     unsigned int l2_size,	/* number of entries in l2 table */
	     unsigned int hist,	        /* number of history depth */
	     unsigned int xor)		/* VHP xor flag in two level */
{
  struct dpred *pred;

  if (!(pred = calloc(1, sizeof(struct dpred))))
    fatal("out of virtual memory");

  pred->class = class;

  switch (class) {
  case DPredLast:
    pred->datapred.last = 
      dpred_data_create(class, l1_size, counter, l2_size, hist, 0);
    break;

  case DPredStride: 
     pred->datapred.stride = 
      dpred_data_create(class, l1_size, 0, 0, 0,0);
   break;

  case DPred2Level: 
    pred->datapred.two = 
      dpred_data_create(class, l1_size, counter, l2_size, hist, xor);
    break;

  case DPredHybrid: 
    pred->datapred.hybrid = 
      dpred_data_create(class, l1_size, counter, l2_size, hist, xor);
    break;

  default:
    panic("bogus predictor class");
  }
  return pred;
}


/* create a data predictor table */
struct dpred_data *		/* data predictor table instance */
dpred_data_create (enum dpred_class class,	/* type of predictor to create */
      	     	   unsigned int l1_size,	/* number of entries in l1 table */
	           unsigned int counter,        /* bits number of a counter 
					           counter threshold of two level */
	           unsigned int l2_size,	/* number of entries in l2 table */
	           unsigned int hist,	        /* number of history depth */
	           unsigned int xor)		/* VHP xor flag in two level */
{
  struct dpred_data *pred_data;
  unsigned int i;
  int flipflop;

  if (!(pred_data = calloc(1, sizeof(struct dpred_data))))
    fatal("out of virtual memory");

  pred_data->class = class;

  switch (class) {
  case DPredLast:
    {
      unsigned int ct_size, vpt_size;
      struct dpred_cache *cache_ct, *cache_vpt;

      /* mapping l1 to ct, l2 to vpt */
      ct_size = l1_size;
      vpt_size = l2_size;

      /** dpred_data.config.last **/
      if (!ct_size || (ct_size & (ct_size-1)) != 0)
	fatal("CT size, `%d', must be non-zero and a power of two", 
	      ct_size);
      pred_data->config.last.ct_size = ct_size;
      
      if (!counter)   /* counter size */
	fatal("counter size, `%d', must be non-zero", 
	      counter);
      pred_data->config.last.counter_size = counter;
      
      if (!vpt_size || (vpt_size & (vpt_size-1)) != 0)
	fatal("VPT size, `%d', must be non-zero and a power of two", 
	      vpt_size);
      pred_data->config.last.vpt_size = vpt_size;

      if (!hist)
	fatal("history depth, `%d', must be non-zero", 
	      hist);
      pred_data->config.last.hist = hist;
     
      cache_ct = pred_data->config.last.ct = calloc(1, sizeof(struct dpred_cache));
      if (!pred_data->config.last.ct)
    	fatal("out of virtual memory");

      cache_vpt = pred_data->config.last.vpt = calloc(1, sizeof(struct dpred_cache));
      if (!pred_data->config.last.vpt)
    	fatal("out of virtual memory");

      /* allocate cache table CT */
      cache_ct->sets = ct_size;

      /* ####### supporting only for direct mapping table on current version */
      cache_ct->assoc = 1; 

      cache_ct->table = calloc(cache_ct->sets * cache_ct->assoc, sizeof(struct dpred_ent));
      if (!cache_ct->table)
	fatal("cannot allocate CT table");

      if (cache_ct->assoc > 1)  /* if associative mapping */
	for (i=0; i < (cache_ct->assoc*cache_ct->sets); i++)
	  {
	    if (i % cache_ct->assoc != cache_ct->assoc - 1)
	      cache_ct->table[i].next = &cache_ct->table[i+1];
	    else
	      cache_ct->table[i].next = NULL;
	    
	    if (i % cache_ct->assoc != cache_ct->assoc - 1)
	      cache_ct->table[i+1].prev = &cache_ct->table[i];
	  }

      /* allocate cache table VPT */
      cache_vpt->sets = vpt_size;
      cache_vpt->assoc = 1; 
           /* #########  supporting only for direct mapping table on current version */

      cache_vpt->table = calloc(cache_vpt->sets * cache_vpt->assoc, sizeof(struct dpred_ent));
      if (!cache_vpt->table)
	fatal("cannot allocate VPT table");

      if (cache_vpt->assoc > 1)  /* if associative mapping */
	for (i=0; i < (cache_vpt->assoc*cache_vpt->sets); i++)
	  {
	    /* alocate each VPT table by hitory depth */
            cache_vpt->table[i].ent.last.vpt_data = calloc(hist, sizeof(int));
            if (!cache_vpt->table[i].ent.last.vpt_data)
    		fatal("out of virtual memory");

	    /* lru list */
	    if (i % cache_vpt->assoc != cache_vpt->assoc - 1)
	      cache_vpt->table[i].next = &cache_vpt->table[i+1];
	    else
	      cache_vpt->table[i].next = NULL;
	    
	    if (i % cache_vpt->assoc != cache_vpt->assoc - 1)
	      cache_vpt->table[i+1].prev = &cache_vpt->table[i];
	  }
      else /* direct mapping */
	for (i=0; i < cache_vpt->sets; i++)
	  {
	    /* alocate each VPT table by hitory depth */
            cache_vpt->table[i].ent.last.vpt_data = calloc(hist, sizeof(int));
            if (!cache_vpt->table[i].ent.last.vpt_data)
    		fatal("out of virtual memory");
 	  } 

      /* initialize counters to no prediction */
      flipflop = 0;
      for (i=0; i < (cache_ct->assoc*cache_ct->sets); i++)
        cache_ct->table[i].ent.last.ct_counter = flipflop;
    }
    break;

  case DPredStride: 
    {
      unsigned int vht_size;
      struct dpred_cache *cache_vht;

      /* mapping l1 to vht */
      vht_size = l1_size;

      /** dpred_data.config.last **/
      if (!vht_size || (vht_size & (vht_size-1)) != 0)
	fatal("VHT size, `%d', must be non-zero and a power of two", 
	      vht_size);
      pred_data->config.stride.vht_size = vht_size;

      cache_vht = pred_data->config.stride.vht = calloc(1, sizeof(struct dpred_cache));
      if (!pred_data->config.stride.vht)
    	fatal("out of virtual memory");

      /* allocate cache table VHT */
      cache_vht->sets = vht_size;
      cache_vht->assoc = 1; 
           /* #########  supporting only for direct mapping table on current version */

      cache_vht->table = calloc(cache_vht->sets * cache_vht->assoc, sizeof(struct dpred_ent));
      if (!cache_vht->table)
	fatal("cannot allocate vht table");

      if (cache_vht->assoc > 1)  /* if associative mapping */
	for (i=0; i < (cache_vht->assoc*cache_vht->sets); i++)
	  {
	    /* lru list */
	    if (i % cache_vht->assoc != cache_vht->assoc - 1)
	      cache_vht->table[i].next = &cache_vht->table[i+1];
	    else
	      cache_vht->table[i].next = NULL;
	    
	    if (i % cache_vht->assoc != cache_vht->assoc - 1)
	      cache_vht->table[i+1].prev = &cache_vht->table[i];
	  }
      else /* direct mapping */
	;

      /* initialize state */
      for (i=0; i < (cache_vht->assoc*cache_vht->sets); i++)
        cache_vht->table[i].ent.stride.state = Init;
    }
    break;

  case DPred2Level: 
    { 
      unsigned int vht_size, threshold, pht_size;
      struct dpred_cache *cache_vht, *cache_pht;

      /* mapping l1 to vht, counter to threshpld, l2 to ptn_num */
      vht_size = l1_size;
      threshold = counter;
      pht_size = l2_size;

      /** dpred_data.config.two **/
      if (!vht_size || (vht_size & (vht_size-1)) != 0)
	fatal("VHT size, `%d', must be non-zero and a power of two", 
	      vht_size);
      pred_data->config.two.vht_size = vht_size;
      
      if (!threshold)
	fatal("threshold value of counter, `%d', must be non-zero", 
	      threshold);
      pred_data->config.two.threshold = threshold;     

      if (!pht_size || (pht_size & (pht_size-1)) != 0)
	fatal("PHT size, `%d', must be non-zero and a power of two", 
	      pht_size);
      pred_data->config.two.pht_size = pht_size;

      if (!hist)
	fatal("history depth, `%d', must be non-zero", 
	      hist);
      pred_data->config.two.hist = hist;
      pred_data->config.two.xor = xor;
     

      cache_vht = pred_data->config.two.vht = calloc(1, sizeof(struct dpred_cache));
      if (!pred_data->config.two.vht)
    	fatal("out of virtual memory");

      cache_pht = pred_data->config.two.pht = calloc(1, sizeof(struct dpred_cache));
      if (!pred_data->config.two.pht)
    	fatal("out of virtual memory");

      /** allocate cache table VHT **/
      cache_vht->sets = vht_size;

      /* ####### supporting only for direct mapping table on current version */
      cache_vht->assoc = 1; 

      cache_vht->table = calloc(cache_vht->sets * cache_vht->assoc, sizeof(struct dpred_ent));
      if (!cache_vht->table)
	fatal("cannot allocate VHT table");

      if (cache_vht->assoc > 1)  /* if associative mapping */
	for (i=0; i < (cache_vht->assoc*cache_vht->sets); i++)
	  {
            /* allocate lru_info according to hist */
            cache_vht->table[i].ent.two.vht_val.lru_info = calloc(hist, sizeof(char));
            if (!cache_vht->table[i].ent.two.vht_val.lru_info)
    		fatal("out of virtual memory");

            /* allocate values according to hist */
            cache_vht->table[i].ent.two.vht_val.values = calloc(hist, sizeof(int));
            if (!cache_vht->table[i].ent.two.vht_val.values)
    		fatal("out of virtual memory");

	    /* cache lru list */
	    if (i % cache_vht->assoc != cache_vht->assoc - 1)
	      cache_vht->table[i].next = &cache_vht->table[i+1];
	    else
	      cache_vht->table[i].next = NULL;    
	    if (i % cache_vht->assoc != cache_vht->assoc - 1)
	      cache_vht->table[i+1].prev = &cache_vht->table[i];
	  }
      else /* direct mapping */
	for (i=0; i < cache_vht->sets; i++)
	  {
            /* allocate lru_info according to hist */
            cache_vht->table[i].ent.two.vht_val.lru_info = calloc(hist, sizeof(char));
            if (!cache_vht->table[i].ent.two.vht_val.lru_info)
    		fatal("out of virtual memory");

            /* allocate values according to hist */
            cache_vht->table[i].ent.two.vht_val.values = calloc(hist, sizeof(int));
            if (!cache_vht->table[i].ent.two.vht_val.values)
    		fatal("out of virtual memory");
 	  } 

      /** allocate cache table PHT **/     
      cache_pht->sets = pht_size;
      cache_pht->assoc = 1;          /* only for direct mapping table */

      cache_pht->table = calloc(pht_size, sizeof(struct dpred_ent));
      if (!cache_pht->table)
	fatal("cannot allocate PHT table");

      for (i=0; i < pht_size; i++)
	  {
	    /* allocate each PHT table by hitory depth */
            cache_pht->table[i].ent.two.pht_val = calloc(hist, sizeof(char));
            if (!cache_pht->table[i].ent.two.pht_val)
    		fatal("out of virtual memory");
 	  } 
      break;
    }

  case DPredHybrid: 
    { 
      unsigned int vht_size, threshold, pht_size;
      struct dpred_cache *cache_vht, *cache_pht;

      /* mapping l1 to vht, counter to threshpld, l2 to ptn_num */
      vht_size = l1_size;
      threshold = counter;
      pht_size = l2_size;

      /** dpred_data.config.hybrid **/
      if (!vht_size || (vht_size & (vht_size-1)) != 0)
	fatal("VHT size, `%d', must be non-zero and a power of two", 
	      vht_size);
      pred_data->config.hybrid.vht_size = vht_size;
      
      if (!threshold)
	fatal("threshold value of counter, `%d', must be non-zero", 
	      threshold);
      pred_data->config.hybrid.threshold = threshold;     

      if (!pht_size || (pht_size & (pht_size-1)) != 0)
	fatal("PHT size, `%d', must be non-zero and a power of two", 
	      pht_size);
      pred_data->config.hybrid.pht_size = pht_size;

      if (!hist)
	fatal("history depth, `%d', must be non-zero", 
	      hist);
      pred_data->config.hybrid.hist = hist;
      pred_data->config.hybrid.xor = xor;
     

      cache_vht = pred_data->config.hybrid.vht = calloc(1, sizeof(struct dpred_cache));
      if (!pred_data->config.hybrid.vht)
    	fatal("out of virtual memory");

      cache_pht = pred_data->config.hybrid.pht = calloc(1, sizeof(struct dpred_cache));
      if (!pred_data->config.hybrid.pht)
    	fatal("out of virtual memory");

      /** allocate cache table VHT **/
      cache_vht->sets = vht_size;

      /* ####### supporting only for direct mapping table on current version */
      cache_vht->assoc = 1; 

      cache_vht->table = calloc(cache_vht->sets * cache_vht->assoc, sizeof(struct dpred_ent));
      if (!cache_vht->table)
	fatal("cannot allocate VHT table");

      if (cache_vht->assoc > 1)  /* if associative mapping */
	for (i=0; i < (cache_vht->assoc*cache_vht->sets); i++)
	  {
            /* allocate lru_info according to hist */
            cache_vht->table[i].ent.hybrid.vht_val.lru_info = calloc(hist, sizeof(char));
            if (!cache_vht->table[i].ent.hybrid.vht_val.lru_info)
    		fatal("out of virtual memory");

            /* allocate values according to hist */
            cache_vht->table[i].ent.hybrid.vht_val.values = calloc(hist, sizeof(int));
            if (!cache_vht->table[i].ent.hybrid.vht_val.values)
    		fatal("out of virtual memory");

	    /* cache lru list */
	    if (i % cache_vht->assoc != cache_vht->assoc - 1)
	      cache_vht->table[i].next = &cache_vht->table[i+1];
	    else
	      cache_vht->table[i].next = NULL;    
	    if (i % cache_vht->assoc != cache_vht->assoc - 1)
	      cache_vht->table[i+1].prev = &cache_vht->table[i];
	  }
      else /* direct mapping */
	for (i=0; i < cache_vht->sets; i++)
	  {
            /* allocate lru_info according to hist */
            cache_vht->table[i].ent.hybrid.vht_val.lru_info = calloc(hist, sizeof(char));
            if (!cache_vht->table[i].ent.hybrid.vht_val.lru_info)
    		fatal("out of virtual memory");

            /* allocate values according to hist */
            cache_vht->table[i].ent.hybrid.vht_val.values = calloc(hist, sizeof(int));
            if (!cache_vht->table[i].ent.hybrid.vht_val.values)
    		fatal("out of virtual memory");
 	  } 

      /** allocate cache table PHT **/     
      cache_pht->sets = pht_size;
      cache_pht->assoc = 1;          /* only for direct mapping table */

      cache_pht->table = calloc(pht_size, sizeof(struct dpred_ent));
      if (!cache_pht->table)
	fatal("cannot allocate PHT table");

      for (i=0; i < pht_size; i++)
	  {
	    /* allocate each PHT table by hitory depth */
            cache_pht->table[i].ent.hybrid.pht_val = calloc(hist, sizeof(char));
            if (!cache_pht->table[i].ent.hybrid.pht_val)
    		fatal("out of virtual memory");
 	  } 
      break;
    }

  default:
    panic("bogus data predictor table class");
  }

  return pred_data;
}


/* print data predictor configuration */
void
dpred_data_config(
  struct dpred_data *pred_data, /* data predictor instance */
  char name[],                  /* predictor name */
  FILE *stream)	        	/* output stream */
{
  switch (pred_data->class) {
  case DPredLast:
    fprintf(stream,
      "pred_data: %s: %d ct-sz, %d counter-sz, %d vpt-sz, %d hist-sz\n",
      name, pred_data->config.last.ct_size, pred_data->config.last.counter_size,
      pred_data->config.last.vpt_size, pred_data->config.last.hist);
    fprintf(stderr,
      "pred_data: %s: %d ct-sz, %d counter-sz, %d vpt-sz, %d hist-sz\n",
      name, pred_data->config.last.ct_size, pred_data->config.last.counter_size,
      pred_data->config.last.vpt_size, pred_data->config.last.hist);
    break;

  case DPredStride: 
    fprintf(stream,
      "pred_data: %s: %d vht-sz\n", name, pred_data->config.stride.vht_size);
    fprintf(stderr,
      "pred_data: %s: %d vht-sz\n", name, pred_data->config.stride.vht_size);
    break;

  case DPred2Level: 
    fprintf(stream,
      "pred_data: %s: %d vht-sz, %d threshold, %d pht-sz, %d hist-sz, %d xor\n",
      name, pred_data->config.two.vht_size, pred_data->config.two.threshold,
      pred_data->config.two.pht_size, pred_data->config.two.hist, pred_data->config.two.xor);
    fprintf(stderr,
      "pred_data: %s: %d vht-sz, %d threshold, %d pht-sz, %d hist-sz, %d xor\n",
      name, pred_data->config.two.vht_size, pred_data->config.two.threshold,
      pred_data->config.two.pht_size, pred_data->config.two.hist, pred_data->config.two.xor);
    break;

  case DPredHybrid: 
    fprintf(stream,
      "pred_data: %s: %d vht-sz, %d threshold, %d pht-sz, %d hist-sz, %d xor\n",
      name, pred_data->config.hybrid.vht_size, pred_data->config.hybrid.threshold,
      pred_data->config.hybrid.pht_size, pred_data->config.hybrid.hist, pred_data->config.hybrid.xor);
    fprintf(stderr,
      "pred_data: %s: %d vht-sz, %d threshold, %d pht-sz, %d hist-sz, %d xor\n",
      name, pred_data->config.hybrid.vht_size, pred_data->config.hybrid.threshold,
      pred_data->config.hybrid.pht_size, pred_data->config.hybrid.hist, pred_data->config.hybrid.xor);
    break;
  default:
    panic("bogus data predictor class");
  }
}


/* print predictor configuration */
void
dpred_config(struct dpred *pred,	/* data predictor instance */
	     FILE *stream)		/* output stream */
{
  switch (pred->class) {
  case DPredLast:
    dpred_data_config (pred->datapred.last, "last", stream);
    break;

  case DPredStride: ;
    dpred_data_config (pred->datapred.stride, "stride", stream);
    break;

  case DPred2Level: ;
    dpred_data_config (pred->datapred.two, "2lev", stream);
    break;

  case DPredHybrid: ;
    dpred_data_config (pred->datapred.hybrid, "hybrid", stream);
    break;

  default:
    panic("bogus data predictor class");
  }
}

/* ? print predictor stats */
void
dpred_stats(struct dpred *pred,		/* data predictor instance */
	    FILE *stream)		/* output stream */
{
  fprintf(stream, "pred: total-prediction rate = %f\n",
	  (double)pred->data_hits/(double)pred->lookups);
  fprintf(stream, "pred: data-prediction rate = %f\n",
	  (double)pred->data_hits/(double)(pred->data_hits+pred->misses));
  fprintf(stream, "pred:  no prediction rate = %f\n",
	  (double)pred->no_hits/(double)(pred->no_hits+pred->no_misses));
/*  fprintf(stream, "pred:  aliasing rate = %f\n",
	  (double)pred->alias/(double)pred->lookups);
*/
  fprintf(stream, "pred:  l1 cache miss rate = %f\n",
	  (double)pred->l1_misses/(double)(pred->lookups));
  fprintf(stream, "pred:  l2 cache miss rate = %f\n",
	  (double)pred->l2_misses/(double)(pred->lookups));
}

/* register data predictor stats */
void
dpred_reg_stats(struct dpred *pred,	/* data predictor instance */
		struct stat_sdb_t *sdb)	/* stats database */
{
  char buf[512], buf1[512], *name;

  /* get a name for this predictor */
  switch (pred->class)
    {
    case DPredLast:
      name = "dpred_last";
      break;
    case DPredStride:
      name = "dpred_stride";
      break;
    case DPred2Level:
      name = "dpred_2lev";
      break;
    case DPredHybrid:
      name = "dpred_hybrid";
      break;
    default:
      panic("bogus data predictor class");
    }

  sprintf(buf, "%s.lookups", name);
  stat_reg_counter(sdb, buf, "total number of dpred lookups",
		   &pred->lookups, 0, NULL);

  /* prediction hit ratio */
  sprintf(buf, "%s.dpred_total_hit_rate", name);
  sprintf(buf1, "%s.data_hits / %s.lookups", name, name);
  stat_reg_formula(sdb, buf,
		   "total data prediction hit rate (i.e., data-hits/lookups)",
		   buf1, "%9.4f");

  sprintf(buf, "%s.updates", name);
  sprintf(buf1, "%s.data_hits + %s.misses", name, name);
  stat_reg_formula(sdb, buf, "total number of updates(i.e.,data-hits+data-misses)", 
		   buf1, "%11.0f");

  sprintf(buf, "%s.dpred_pred_rate", name);
  sprintf(buf1, "%s.updates / %s.lookups", name, name);
  stat_reg_formula(sdb, buf,
		   "data-prediction rate (i.e., updates/lookups)",
		   buf1, "%9.4f");

  sprintf(buf, "%s.data_hits", name);
  stat_reg_counter(sdb, buf, "total number of data-predicted hits", 
		   &pred->data_hits, 0, NULL);

  sprintf(buf, "%s.misses", name);
  stat_reg_counter(sdb, buf, "total number of misses", &pred->misses, 0, NULL);

  sprintf(buf, "%s.dpred_data_hit_rate", name);
  sprintf(buf1, "%s.data_hits / %s.updates", name, name);
  stat_reg_formula(sdb, buf,
		   "data-prediction rate (i.e., data-hits/updates)",
		   buf1, "%9.4f");

  /* no prediction hit ratio */
  sprintf(buf, "%s.no_updates", name);
  sprintf(buf1, "%s.no_hits + %s.no_misses", name, name);
  stat_reg_formula(sdb, buf, "total number of no_updates", buf1, "%11.0f");

  sprintf(buf, "%s.no_hits", name);
  stat_reg_counter(sdb, buf, 
		   "total number of no-predicted hits ",
		   &pred->no_hits, 0, NULL);

  sprintf(buf, "%s.no_misses", name);
  stat_reg_counter(sdb, buf, "total number of no-misses", &pred->no_misses, 0, NULL);

  sprintf(buf, "%s.dpred_no_hit_rate", name);
  sprintf(buf1, "%s.no_hits / %s.no_updates", name, name);
  stat_reg_formula(sdb, buf,
		  "data no-prediction rate (i.e., no-hits/no-updates)",
		  buf1, "%9.4f");

  /* aliasing rate */
  if (pred->class == DPred2Level || pred->class == DPredHybrid) {
  	sprintf(buf, "%s.alias", name);
  	stat_reg_counter(sdb, buf, 
			   "total number of aliasing ",
			   &pred->alias, 0, NULL);

  	sprintf(buf, "%s.alias_rate", name);
  	sprintf(buf1, "%s.alias / %s.updates", name, name);
  	stat_reg_formula(sdb, buf,
			   "aliasing rate (i.e., alias/updates)",
			   buf1, "%9.4f");

  	sprintf(buf, "%s.alias_hits", name);
  	stat_reg_counter(sdb, buf, 
			   "hit number on aliasing ",
			   &pred->alias_hits, 0, NULL);

  	sprintf(buf, "%s.alias_misses", name);
  	stat_reg_counter(sdb, buf, 
			   "miss number on aliasing ",
			   &pred->alias_misses, 0, NULL);

  	sprintf(buf, "%s.alias_hit_rate", name);
  	sprintf(buf1, "%s.alias_hits / (%s.alias_hits+%s.alias_misses)", name, name,name);
  	stat_reg_formula(sdb, buf,
			   "data hits on aliasing (i.e., alias_hits/(alias_hits+alias_misses))",
			   buf1, "%9.4f");
  }

  /* cache miss rate */
  sprintf(buf, "%s.l1_misses", name);
  stat_reg_counter(sdb, buf, "total number of l1 cache misses", &pred->l1_misses, 0, NULL);

  sprintf(buf, "%s.l1_miss_rate", name);
  sprintf(buf1, "%s.l1_misses / %s.lookups", name, name);
  stat_reg_formula(sdb, buf,
		   "l1 cache miss rate (i.e., l1_misses/lookups)",
		   buf1, "%9.4f");

  sprintf(buf, "%s.l2_misses", name);
  stat_reg_counter(sdb, buf, "total number of l2 cache misses", &pred->l2_misses, 0, NULL);

  sprintf(buf, "%s.l2_miss_rate", name);
  sprintf(buf1, "%s.l2_misses / %s.lookups", name, name);
  stat_reg_formula(sdb, buf,
		   "l2 cache miss rate (i.e., l2_misses/lookups)",
		   buf1, "%9.4f");
}

/* reset stats after priming, if appropriate */
void dpred_after_priming(struct dpred *dpred)
{
  if (dpred == NULL)
    return;

  dpred->lookups = 0;
  dpred->data_hits = 0;
  dpred->misses = 0;
  dpred->no_hits = 0;
  dpred->no_misses = 0;
  dpred->l1_misses = 0;
  dpred->l2_misses = 0;
  dpred->alias = 0;
  dpred->alias_hits = 0;
  dpred->alias_misses = 0;
}

/* lookup cache table */
struct dpred_ent *				/* pointer to table entry */
dpred_cache_lookup (
  struct dpred_cache *pred_cache,/* cache table instance */
  SS_ADDR_TYPE daddr)            /* PC address */
{
  struct dpred_ent *pent = NULL;
  int index, i;

  /* Get a pointer into the cache */
  index = (daddr >> 3) & (pred_cache->sets - 1);

  if (pred_cache->assoc > 1)  /* associative mapping */
    {
      index *= pred_cache->assoc;

      /* Now we know the set; look for a PC match */
      for (i = index; i < (index+pred_cache->assoc) ; i++)
	if (pred_cache->table[i].addr == daddr)  
	  {
	    /* match */
	    pent = &pred_cache->table[i];
	    break;
	  }
    }	
  else	/* direct mapping */
    {
      pent = &pred_cache->table[index];
      if (pent->addr != daddr)
	pent = NULL;
    }  
  
  return pent;  
}


/* probe a predictor for a data value, the predictor is probed
   with instruction PC address */
SS_WORD_TYPE				/* predicted data value */
dpred_lookup(struct dpred *pred,	/* data predictor instance */
	     SS_ADDR_TYPE pcaddr,	/* instruction PC address */
	     enum ss_opcode op,		/* opcode of onstruction */
	     char *no_pred, 		/* predictor state pointer */
	     struct dpred_ent **pptbl1,	/* pointer of table 1 entry */
	     struct dpred_ent **pptbl2)	/* pointer of table 2 entry */
{
  struct dpred_ent *ptbl1, *ptbl2;
 
  /* if this is not a integer register write instruction, make no predict */
  if (!is_PRED) 
    {
       *no_pred = 4;	/* set it to no prediction */
       return 0;
    }

  pred->lookups++;

  ptbl1 = NULL;
  ptbl2 = NULL;
 
 /* get a pointer to prediction table entry and return predicted data */
  switch (pred->class) {
  case DPredLast:
    /* get CT entry pointer */
    ptbl1 = dpred_cache_lookup (pred->datapred.last->config.last.ct, pcaddr);
    *pptbl1 = ptbl1;
    if (!ptbl1)   /* CT cache miss, no prediction */
       {
	*no_pred = 2;
	return 0;		
       }
    /* get VPT entry pointer */
    ptbl2 = dpred_cache_lookup (pred->datapred.last->config.last.vpt, pcaddr);
    *pptbl2 = ptbl2;
    if (!ptbl2)   /* VPT cache miss, no prediction */
       {
	*no_pred = 3;
	return 0;		
       }

    if (ptbl1->ent.last.ct_counter >= 2) 
      {  /* prediction */
	*no_pred = 0;
      }
    else
      {  /* no prediction by counter value */
	*no_pred = 1;
      }	

    return ptbl2->ent.last.vpt_data[0];
		/* ######## only if hist = 1 */
    break;

  case DPredStride: ;
    /* get VHT entry pointer */
    ptbl1 = dpred_cache_lookup (pred->datapred.stride->config.stride.vht, pcaddr);
    *pptbl1 = ptbl1;
    if (!ptbl1)   /* VHT cache miss, no prediction */
       {
	*no_pred = 2;
	return 0;		
       }

    /* set no_pred by the state */
    if (ptbl1->ent.stride.state == Steady) 
      {  /* prediction */
	*no_pred = 0;
      }
    else
      {  /* no prediction by the state */
	*no_pred = 1;
      }	

    return (ptbl1->ent.stride.value + ptbl1->ent.stride.stride);
    break;

  case DPred2Level: 
  {
    int index, i, j, max_count;

    /* get VHT entry pointer */
    ptbl1 = dpred_cache_lookup (pred->datapred.two->config.two.vht, pcaddr);
    *pptbl1 = ptbl1;
    if (!ptbl1)   /* VHT cache miss, no prediction */
       {
	*no_pred = 2;
	return 0;		
       }

    /* get PHT entry pointer */
    if (!pred->datapred.two->config.two.xor)   
        index = ptbl1->ent.two.vht_val.vhp & (pred->datapred.two->config.two.pht_size - 1);
    else { /* xor */
  	unsigned int pctemp;
	pctemp = pcaddr;
	switch (pred->datapred.two->config.two.xor) {
	  case 6 : pctemp &= 0x3f; break;
	  case 8 : pctemp &= 0xff; break;	  
	  case 12 : pctemp &= 0xfff; break;
	  case 16 : pctemp &= 0xffff; break;
	  default : fatal("illegal xor bit size");
	}
        index = (ptbl1->ent.two.vht_val.vhp ^ pctemp) 
                & (pred->datapred.two->config.two.pht_size - 1);
    }

    ptbl2 = &pred->datapred.two->config.two.pht->table[index];
    *pptbl2 = ptbl2;

    if (!ptbl2)   /* PHT miss, error */
	fatal("PHT table miss");


    /* find max counter in PHT */
    max_count = ptbl2->ent.two.pht_val[0];
    j = 0;
    for (i=1; i < pred->datapred.two->config.two.hist; i++)
	if (max_count < ptbl2->ent.two.pht_val[i]) 
         {
    	     max_count = ptbl2->ent.two.pht_val[i];
	     j = i;
	 }

    if (max_count >= pred->datapred.two->config.two.threshold)
      {  /* prediction */
	*no_pred = 0;
         /* aliasing check */
    	if (ptbl2->addr != pcaddr) {
	 /* aliasing occur, count aliasing */
         	pred->alias++;
    	}
        return ptbl1->ent.two.vht_val.values[j];
      }
    else
      {  /* no prediction by counter value */
	*no_pred = 1;
      }	

    return 0;
    break;
  }

  case DPredHybrid: 
  {
    int index, i, j, max_count;

    /* get VHT entry pointer */
    ptbl1 = dpred_cache_lookup (pred->datapred.hybrid->config.hybrid.vht, pcaddr);
    *pptbl1 = ptbl1;
    if (!ptbl1)   /* VHT cache miss, no prediction */
       {
	*no_pred = 2;
	return 0;		
       }

    /* get PHT entry pointer */
    if (!pred->datapred.hybrid->config.hybrid.xor)   
        index = ptbl1->ent.hybrid.vht_val.vhp & (pred->datapred.hybrid->config.hybrid.pht_size - 1);
    else { /* xor */
  	unsigned int pctemp;
	pctemp = pcaddr;
	switch (pred->datapred.hybrid->config.hybrid.xor) {
	  case 6 : pctemp &= 0x3f; break;
	  case 8 : pctemp &= 0xff; break;	  
	  case 12 : pctemp &= 0xfff; break;
	  case 16 : pctemp &= 0xffff; break;
	  default : fatal("illegal xor bit size");
	}
        index = (ptbl1->ent.hybrid.vht_val.vhp ^ pctemp) 
                & (pred->datapred.hybrid->config.hybrid.pht_size - 1);
    }

    ptbl2 = &pred->datapred.hybrid->config.hybrid.pht->table[index];
    *pptbl2 = ptbl2;

    if (!ptbl2)   /* PHT miss, error */
	fatal("PHT table miss");

    /* find max counter in PHT */
    max_count = ptbl2->ent.hybrid.pht_val[0];
    j = 0;
    for (i=1; i < pred->datapred.hybrid->config.hybrid.hist; i++)
	if (max_count < ptbl2->ent.hybrid.pht_val[i]) 
         {
    	     max_count = ptbl2->ent.hybrid.pht_val[i];
	     j = i;
	 }

    if (max_count >= pred->datapred.hybrid->config.hybrid.threshold)
      {  /* two level prediction */
	*no_pred = 0;
    	/* aliasing check */
    	if (ptbl2->addr != pcaddr) {
		/* aliasing occur, count aliasing */
        	pred->alias++;
    	}
	pred->datapred.hybrid->config.hybrid.is_stride = 0;
    	return ptbl1->ent.hybrid.vht_val.values[j];
      }
    else
      { /* stride prediction  */
	int mru, ht;

	pred->datapred.hybrid->config.hybrid.is_stride = 1;  /* set flag */
    	/* set no_pred by the state */
    	if (ptbl1->ent.hybrid.vht_val.state == Steady) 
           {  	/* prediction */
	 	*no_pred = 0;
		ht = pred->datapred.hybrid->config.hybrid.hist;  /* history */
		mru = ptbl1->ent.hybrid.vht_val.lru_info[ht-1];  /* MRU index */
    		return (ptbl1->ent.hybrid.vht_val.values[mru] +
                           ptbl1->ent.hybrid.vht_val.stride);
           }
    	else
      	   {  /* no prediction by the state */
		*no_pred = 1;
       	   }	
      }	
    return 0;
    break;
  }

  default:
    panic("bogus data predictor class");
  } 
}



/* lookup cache table and update LRU state */
struct dpred_ent *				/* pointer to table entry */
dpred_cache_lru_update (
  struct dpred_cache *pred_cache,/* cache table instance */
  SS_ADDR_TYPE daddr)            /* PC address */
{
  struct dpred_ent *ptbl = NULL;
  struct dpred_ent *lruhead = NULL, *lruitem = NULL;
  int index, i;

  /* Get a pointer into the cache */
  index = (daddr >> 3) & (pred_cache->sets - 1);

  if (pred_cache->assoc > 1)  /* associative mapping */
   {
      index *= pred_cache->assoc;

      /* Now we know the set; look for a PC match; also identify
         MRU and LRU items */
      for (i = index; i < (index+pred_cache->assoc) ; i++)
       {
	if (pred_cache->table[i].addr == daddr)  
	  {
	    /* match */
	    assert(!ptbl);
	    ptbl = &pred_cache->table[i];
	  }
        dassert(pred_cache->table[i].prev 
	        != pred_cache->table[i].next);
        if (pred_cache->table[i].prev == NULL)
	 {
	   /* this is the head of the lru list, ie current MRU item */
	   dassert(lruhead == NULL);
	   lruhead = &pred_cache->table[i];
	 }
	 if (pred_cache->table[i].next == NULL)
	  {
	    /* this is the tail of the lru list, ie the LRU item */
	    dassert(lruitem == NULL);
	    lruitem = &pred_cache->table[i];
	   }
       }	/* for */
       dassert(lruhead && lruitem);

       if (!ptbl)
    	/* missed in cache; choose the LRU item in this set as the victim */
	ptbl = lruitem;	
       	/* else hit, and ptbl points to matching CT entry */
	  
  	/* Update LRU state: selected item, whether selected because it
	 * matched or because it was LRU and selected as a victim, becomes 
	 * MRU */
       if (ptbl != lruhead)
	{
      	   /* this splices out the matched entry... */
	   if (ptbl->prev)
		ptbl->prev->next = ptbl->next;
           if (ptbl->next)
		ptbl->next->prev = ptbl->prev;
	   /* ...and this puts the matched entry at the head of the list */
	   ptbl->next = lruhead;
	   ptbl->prev = NULL;
	   lruhead->prev = ptbl;
	   dassert(ptbl->prev || ptbl->next);
	   dassert(ptbl->prev != ptbl->next);
	 }
	 /* else ptbl is already MRU item; do nothing */
  }
  else  /* direct mapping */
      ptbl = &pred_cache->table[index];

  return ptbl;
}
 

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
	     enum ss_opcode op)		/* opcode of instruction */
{

  struct dpred_ent *ptbl1 = NULL, *ptbl2 = NULL;
  struct dpred_cache *pred_cache;
  int index, i, j, hist, matched_index;

  /* don't change dpred state for non-predicted instructions */
  if (!is_PRED)
    return;

  /* keep stats about data prediction */
  if (!no_pred)	/* if predict */
    {
  	if (correct)
    		pred->data_hits++;
	else
	    	pred->misses++;
    }
  else if (no_pred == 1) 
		/* no pred by counter(last,2lev) or by state(stride) */
    {
  	if (correct)
    		pred->no_misses++;
	else
	    	pred->no_hits++;
    }
  else if (no_pred == 2) 
		/* no pred by CT miss(last) or by VHT miss(stride,2lev) */
    {	;
    }
  else if (no_pred == 3) 
		/* no pred by VPT miss(last) */
    {	;
    }

   /**** update cache table *****/
  switch (pred->class) {
  case DPredLast:
    /* get CT entry pointer and update LRU */
    pred_cache = pred->datapred.last->config.last.ct;
    ptbl1 = dpred_cache_lru_update(pred_cache, daddr);

    /* update CT */
    if (ptbl1)
     {
      /* update current information */
      if (ptbl1->addr == daddr) /* CT cache hit */
	{
	  if (correct)	/* prediction hit */
	    {
		/* ####### only if ct_size == 2 */
		if (ptbl1->ent.last.ct_counter < 3) 
		   ++ptbl1->ent.last.ct_counter;
	    }
          else	/* prediction miss */
	    {
		/* ####### only if ct_size == 2 */
		if (ptbl1->ent.last.ct_counter > 0) 
		   --ptbl1->ent.last.ct_counter;
	    }
	}
      else	/* cache miss */
	{
	  assert(no_pred == 2);
	  pred->l1_misses++;
	  /* enter a new entry in the table */
	  ptbl1->addr = daddr;
	  ptbl1->op = op;
	  ptbl1->ent.last.ct_counter = 0;
	}
     }

    /* get VPT entry pointer and update LRU */
    pred_cache = pred->datapred.last->config.last.vpt;
    ptbl2 = dpred_cache_lru_update(pred_cache, daddr);

    /* update VPT */
    if (ptbl2)
     {
      /* update current information */
      if (ptbl2->addr == daddr) /* VPT cache hit */
	{
	  if (!correct)	/* prediction miss */
	   ptbl2->ent.last.vpt_data[0] = data;  /* ######## only if hist == 1 */
	}
      else	/* cache miss */
	{
	  assert(no_pred == 3 || no_pred == 2);
	  pred->l2_misses++;
	  /* enter a new entry in the table */
	  ptbl2->addr = daddr;
	  ptbl2->op = op;
	  ptbl2->ent.last.vpt_data[0] = data;  /* ######## only if hist == 1 */
	}
     } 

    break;

  case DPredStride:
    /* get VHT entry pointer and update LRU */
    pred_cache = pred->datapred.stride->config.stride.vht;
    ptbl1 = dpred_cache_lru_update(pred_cache, daddr);

    /* update VHT */
    if (ptbl1)
     {
      /* update current information */
      if (ptbl1->addr == daddr) /* VHT cache hit */
	{

	  if (ptbl1->ent.stride.state == Init)  /* Init state */
            {
 		ptbl1->ent.stride.state = Transient;
	  	ptbl1->ent.stride.value = data;
	  	ptbl1->ent.stride.stride = 0;  /* any stride */
            }
	  else  if (ptbl1->ent.stride.state == Transient)  /* Transient state */
            {
	        int OldStride, CurrentStride;

                OldStride = ptbl1->ent.stride.stride;
		CurrentStride = data - ptbl1->ent.stride.value;
		
		if (OldStride == CurrentStride)  /* same as previous stride */
		  {
 			ptbl1->ent.stride.state = Steady;
	  		ptbl1->ent.stride.value = data;
		  }
	 	else  /* different from previous stride */
		  {
	  		ptbl1->ent.stride.value = data;
	  		ptbl1->ent.stride.stride = CurrentStride;  
		  }
            }  /* Transient state */
	  else  if (ptbl1->ent.stride.state == Steady)  /* Steady state */
            {
	        int OldStride, CurrentStride;

                OldStride = ptbl1->ent.stride.stride;
		CurrentStride = data - ptbl1->ent.stride.value;
		
		if (OldStride == CurrentStride)  /* same as previous stride */
		  {
	  		ptbl1->ent.stride.value = data;
		  }
	 	else  /* different from previous stride */
		  {
 			ptbl1->ent.stride.state = Transient;
	  		ptbl1->ent.stride.value = data;
	  		ptbl1->ent.stride.stride = CurrentStride;  
		  }
            }  /* Steady state */
	}

      else	/* VHT che miss */
	{
	  assert(no_pred == 2);
	  pred->l1_misses++;
	  /* enter a new entry in the table */
	  ptbl1->addr = daddr;
	  ptbl1->op = op;

	  /* set inital state */
	  ptbl1->ent.stride.state = Init;	
	  ptbl1->ent.stride.value = data;
	}
     }

    break;

  case DPred2Level: 
   { 
     int hist_num /* bit # of encoded hist */, hist_temp;

    /* get VHT entry pointer and update LRU */
    pred_cache = pred->datapred.two->config.two.vht;
    ptbl1 = dpred_cache_lru_update(pred_cache, daddr);

    /* get PHT entry pointer */
    if (!pred->datapred.two->config.two.xor)   
        index = ptbl1->ent.two.vht_val.vhp & (pred->datapred.two->config.two.pht_size - 1);
    else { /* xor */
  	unsigned int pctemp;
	pctemp = daddr;
	switch (pred->datapred.two->config.two.xor) {
	  case 6 : pctemp &= 0x3f; break;
	  case 8 : pctemp &= 0xff; break;	  
	  case 12 : pctemp &= 0xfff; break;
	  case 16 : pctemp &= 0xffff; break;
	  default : fatal("illegal xor bit size");
	}
        index = (ptbl1->ent.two.vht_val.vhp ^ pctemp) 
                & (pred->datapred.two->config.two.pht_size - 1);
    }

    ptbl2 = &pred->datapred.two->config.two.pht->table[index];

    hist = pred->datapred.two->config.two.hist;
    switch (hist) { /* determine bit # of hist */
	case 2 : ;
	case 4 : hist_num = 2; hist_temp = 3; break;
	case 6 : ;
	case 8 : hist_num = 3; hist_temp = 7; break;
	case 10 : ;
	case 12 : ;
	case 14 : ;
	case 16 : hist_num = 4; hist_temp = 15; break;
	default : fatal("illegal hist #");
    }


    /* update VHT */
    if (ptbl1)
     {
      /* update current information */
      if (ptbl1->addr == daddr) /* VHT cache hit */
	{
	    /* PHT, counter update */
	    for (i=0; i < hist; i++) /* find the index of data values in VHT */
	    	if (ptbl1->ent.two.vht_val.values[i] == data) 
			break;  
 	    matched_index = i;

	    for (j=0; j < hist; j++) { 
	    	if (matched_index == j) { /* matched counter, increment by 3 */
			if (ptbl2->ent.two.pht_val[j] < 9)
			  ptbl2->ent.two.pht_val[j] += 3;
			else  /* saturated value */
			  ptbl2->ent.two.pht_val[j] = 12;  
		}
		else {  /* unmatched counter, decrement by 1 */
			if (ptbl2->ent.two.pht_val[j] > 0)
			  ptbl2->ent.two.pht_val[j]--;
		}
	    }

   	   if (!no_pred && ptbl2->addr != daddr) { /* prediction & aliasing occur */
 		if (correct)
	   		pred->alias_hits++;
		else
	   		pred->alias_misses++;
    	   }

   	    /* LRU info update */
            if (matched_index < hist) {
		/* if there has the value of resolved data in  VHT */
		for (i=0; i < hist; i++) /* find matched_index in lru_info */
		    if (ptbl1->ent.two.vht_val.lru_info[i] == matched_index)
			break;
		for (j=i; j < hist-1; j++)  /* shift lru info */
		     ptbl1->ent.two.vht_val.lru_info[j]	= ptbl1->ent.two.vht_val.lru_info[j+1];
		ptbl1->ent.two.vht_val.lru_info[hist-1] = matched_index;  /* MRU value */

		/* VHP update */
		ptbl1->ent.two.vht_val.vhp <<= hist_num;
		ptbl1->ent.two.vht_val.vhp |= (matched_index & hist_temp);
		if (correct) /* PC address for aliasing check */
	    		ptbl2->addr = daddr;	
	    }
	    else {
	        /* if there has no value of resolved data in VHT */
	        i = ptbl1->ent.two.vht_val.lru_info[0];  /* LRU victim */
	        ptbl1->ent.two.vht_val.values[i] = data; /* new data, index i */
	        /* LRU info update */
	        for (j=0; j < hist-1; j++)  /* shift lru info */
	      		ptbl1->ent.two.vht_val.lru_info[j] = ptbl1->ent.two.vht_val.lru_info[j+1];
		ptbl1->ent.two.vht_val.lru_info[j] = i; /* MRU */
		/* VHP update */
		ptbl1->ent.two.vht_val.vhp <<= hist_num;
		ptbl1->ent.two.vht_val.vhp |= (i & hist_temp);
		/* new counter in PHT */
		ptbl2->ent.two.pht_val[i] = 0;
	    }
	}  /* VHT hit */
      else	/* VHT cache miss */
	{
	  assert(no_pred == 2);
	  pred->l1_misses++;
	  /* enter a new entry in the table */
	  ptbl1->addr = daddr;
	  ptbl1->op = op;

  	  /* data value initialization */
	  ptbl1->ent.two.vht_val.values[0] = data;  /* new data, index 0 */
	  for (i=1; i < hist; i++) /* find the index of data values in VHT */
	     ptbl1->ent.two.vht_val.values[i] = 0;  /* initialized data */
   	  /* LRU info initialization */
          for (j=0; j < hist-1; j++)  /* shift lru info */
	     ptbl1->ent.two.vht_val.lru_info[j]	= j+1;
	  ptbl1->ent.two.vht_val.lru_info[j] = 0;  /* MRU value */
	  /* VHP initialization */
          ptbl1->ent.two.vht_val.vhp = 0;
	}
     }  /* ptbl1 */

    break;
   }  /* case DPred2Level */

  case DPredHybrid: 
   { 
     int hist_num /* bit # of encoded hist */, hist_temp;

    /* get VHT entry pointer and update LRU */
    pred_cache = pred->datapred.hybrid->config.hybrid.vht;
    ptbl1 = dpred_cache_lru_update(pred_cache, daddr);

    /* get PHT entry pointer */
    if (!pred->datapred.hybrid->config.hybrid.xor)   
        index = ptbl1->ent.hybrid.vht_val.vhp & (pred->datapred.hybrid->config.hybrid.pht_size - 1);
    else { /* xor */
  	unsigned int pctemp;
	pctemp = daddr;
	switch (pred->datapred.hybrid->config.hybrid.xor) {
	  case 6 : pctemp &= 0x3f; break;
	  case 8 : pctemp &= 0xff; break;	  
	  case 12 : pctemp &= 0xfff; break;
	  case 16 : pctemp &= 0xffff; break;
	  default : fatal("illegal xor bit size");
	}
        index = (ptbl1->ent.hybrid.vht_val.vhp ^ pctemp) 
                & (pred->datapred.hybrid->config.hybrid.pht_size - 1);
    }

    ptbl2 = &pred->datapred.hybrid->config.hybrid.pht->table[index];

    hist = pred->datapred.hybrid->config.hybrid.hist;
    switch (hist) { /* determine bit # of hist */
	case 2 : ;
	case 4 : hist_num = 2; hist_temp = 3; break;
	case 6 : ;
	case 8 : hist_num = 3; hist_temp = 7; break;
	case 10 : ;
	case 12 : ;
	case 14 : ;
	case 16 : hist_num = 4; hist_temp = 15; break;
	default : fatal("illegal hist #");
    }

    /* update VHT */
    if (ptbl1)
     {
      /* update current information */
      if (ptbl1->addr == daddr) /* VHT cache hit */
	{
	    /* PHT, counter update */
	    for (i=0; i < hist; i++) /* find the index of data values in VHT */
	    	if (ptbl1->ent.hybrid.vht_val.values[i] == data) 
			break;  
 	    matched_index = i;

	    for (j=0; j < hist; j++) { 
	    	if (matched_index == j) { /* matched counter, increment by 3 */
			if (ptbl2->ent.hybrid.pht_val[j] < 9)
			  ptbl2->ent.hybrid.pht_val[j] += 3;
			else  /* saturated value */
			  ptbl2->ent.hybrid.pht_val[j] = 12;  
		}
		else {  /* unmatched counter, decrement by 1 */
			if (ptbl2->ent.hybrid.pht_val[j] > 0)
			  ptbl2->ent.hybrid.pht_val[j]--;
		}
	    }

	    if (!no_pred && ptbl2->addr != daddr) {  /* aliasing check */ 
    	     if (!pred->datapred.hybrid->config.hybrid.is_stride) {
			/* if two level prediction */ 
  		if (correct)
	   	    pred->alias_hits++;
		else
	   	    pred->alias_misses++;
	     }
	    }

   	    /* LRU info update */
            if (matched_index < hist) {
		/* if there has the value of resolved data in  VHT */
		for (i=0; i < hist; i++) /* find matched_index in lru_info */
		    if (ptbl1->ent.hybrid.vht_val.lru_info[i] == matched_index)
			break;
		for (j=i; j < hist-1; j++)  /* shift lru info */
		     ptbl1->ent.hybrid.vht_val.lru_info[j]	= ptbl1->ent.hybrid.vht_val.lru_info[j+1];
		ptbl1->ent.hybrid.vht_val.lru_info[hist-1] = matched_index;  /* MRU value */

		/* VHP update */
		ptbl1->ent.hybrid.vht_val.vhp <<= hist_num;
		ptbl1->ent.hybrid.vht_val.vhp |= (matched_index & hist_temp);
		if (correct) /* PC address for aliasing check */
	    		ptbl2->addr = daddr;	
	    }
	    else {
	        /* if there has no value of resolved data in VHT */
	        i = ptbl1->ent.hybrid.vht_val.lru_info[0];  /* LRU victim */
	        ptbl1->ent.hybrid.vht_val.values[i] = data; /* new data, index i */
	        /* LRU info update */
	        for (j=0; j < hist-1; j++)  /* shift lru info */
	      		ptbl1->ent.hybrid.vht_val.lru_info[j] = ptbl1->ent.hybrid.vht_val.lru_info[j+1];
		ptbl1->ent.hybrid.vht_val.lru_info[j] = i; /* MRU */
		/* VHP update */
		ptbl1->ent.hybrid.vht_val.vhp <<= hist_num;
		ptbl1->ent.hybrid.vht_val.vhp |= (i & hist_temp);
		/* new counter in PHT */
		ptbl2->ent.hybrid.pht_val[i] = 0;
	    }

          /* stride update */
	  if (ptbl1->ent.hybrid.vht_val.state == Init)  /* Init state */
            {
 		ptbl1->ent.hybrid.vht_val.state = Transient;
	  	ptbl1->ent.hybrid.vht_val.stride = 0;  /* any stride */
            }
	  else  if (ptbl1->ent.hybrid.vht_val.state == Transient)  /* Transient state */
            {
	        int OldStride, CurrentStride, mru2;

		mru2 = ptbl1->ent.hybrid.vht_val.lru_info[hist-2]; /* already data update */
                OldStride = ptbl1->ent.hybrid.vht_val.stride;
		CurrentStride = data - ptbl1->ent.hybrid.vht_val.values[mru2];
	
		if (OldStride == CurrentStride)  /* same as previous stride */
		  {
 			ptbl1->ent.hybrid.vht_val.state = Steady;
		  }
	 	else  /* different from previous stride */
		  {
	  		ptbl1->ent.hybrid.vht_val.stride = CurrentStride;  
		  }
            }  /* Transient state */
	  else  if (ptbl1->ent.hybrid.vht_val.state == Steady)  /* Steady state */
            {
	        int OldStride, CurrentStride, mru2;

		mru2 = ptbl1->ent.hybrid.vht_val.lru_info[hist-2]; /* already data update */
                OldStride = ptbl1->ent.hybrid.vht_val.stride;
		CurrentStride = data - ptbl1->ent.hybrid.vht_val.values[mru2];
		
		if (OldStride != CurrentStride) /* different from previous stride */
		  {
 			ptbl1->ent.hybrid.vht_val.state = Transient;
	  		ptbl1->ent.hybrid.vht_val.stride = CurrentStride;  
		  }
            }  /* Steady state */
	}  /* VHT hit */
      else	/* VHT cache miss */
	{
	  assert(no_pred == 2);
	  pred->l1_misses++;
	  /* enter a new entry in the table */
	  ptbl1->addr = daddr;
	  ptbl1->op = op;

  	  /* data value initialization */
	  ptbl1->ent.hybrid.vht_val.values[0] = data;  /* new data, index 0 */
	  for (i=1; i < hist; i++) /* find the index of data values in VHT */
	     ptbl1->ent.hybrid.vht_val.values[i] = 0;  /* initialized data */
   	  /* LRU info initialization */
          for (j=0; j < hist-1; j++)  /* shift lru info */
	     ptbl1->ent.hybrid.vht_val.lru_info[j]	= j+1;
	  ptbl1->ent.hybrid.vht_val.lru_info[j] = 0;  /* MRU value */
	  /* VHP initialization */
          ptbl1->ent.hybrid.vht_val.vhp = 0;
          /* stride */
	  ptbl1->ent.hybrid.vht_val.state = Init;	
	}
     }  /* ptbl1 */

    break;
   }  /* case DPredHybrid */

  default:
    panic("bogus data predictor class");
  } 
}



/* print the trace of the table content */
void				
dpred_trace(struct dpred *pred,	/* data predictor instance */
	     SS_ADDR_TYPE pcaddr,	/* instruction PC address */
	     SS_WORD_TYPE pdata,	/* predicted data */
	     SS_WORD_TYPE rdata,	/* resolved data */
	     char *no_pred,		/* predictor state pointer */
	     struct dpred_ent *ptbl1,	/* pointer of table 1 entry */
	     struct dpred_ent *ptbl2)	/* pointer of table 2 entry */
{
  int i, index;

  printf("   ");
  if (!*no_pred) /* prediction */
	printf("%s pDATA:%d rDATA:%d - ", 
                pdata == rdata ? "CRT" : "INCRT",
	        pdata, rdata); 
  else
	printf("NOT_PRED rDATA:%d - ", rdata);
	     	
  switch (pred->class) {
  case DPredLast:
    if (*no_pred == 2)
	printf("CT MISS ");
    else
	printf("CT_cnt:%d ", ptbl1->ent.last.ct_counter);

    if (*no_pred == 3) 
	printf("VPT MISS ");
    else
	printf("VPT_data: %d ", ptbl2->ent.last.vpt_data[0]);
    break;
 
  case DPredStride: 
    if (*no_pred == 2)
	printf("VHT MISS ");
    else {
	printf("state: ");
	if (ptbl1->ent.stride.state == Init)
		printf("Init ");
	else if (ptbl1->ent.stride.state == Transient)
		printf("Transient ");
	else if (ptbl1->ent.stride.state == Steady)
		printf("Steady ");
    	printf("stride: %d ", ptbl1->ent.stride.stride);
    }
    break;
 
  case DPred2Level: 
    if (*no_pred == 2)
	printf("VHT MISS ");
    else { 
        /* printf("no_pred: %d\n", *no_pred); */
        assert(ptbl1);
        assert(ptbl2);
	for (i=0; i < pred->datapred.two->config.two.hist; i++)
		printf("%d ", ptbl1->ent.two.vht_val.values[i]);
       	index = ptbl1->ent.two.vht_val.vhp & 
		(pred->datapred.two->config.two.pht_size - 1);
	printf("vhp:%x - ", index);

        if (pred->datapred.two->config.two.xor) { /* xor */
  	  unsigned int pctemp;
	  pctemp = pcaddr;
	  switch (pred->datapred.two->config.two.xor) {
	  	case 6 : pctemp &= 0x3f; break;
	  	case 8 : pctemp &= 0xff; break;	  
	  	case 12 : pctemp &= 0xfff; break;
	  	case 16 : pctemp &= 0xffff; break;
	  	default : fatal("illegal xor bit size");
	  }
          index = (ptbl1->ent.two.vht_val.vhp ^ pctemp) 
                & (pred->datapred.two->config.two.pht_size - 1);
	  printf("vhp^pc: %x - ", index);
        } /* xor */	
	if (!*no_pred && ptbl2->addr != pcaddr)
		printf("ALIAS(%x) ", ptbl2->addr);
	for (i=0; i < pred->datapred.two->config.two.hist; i++)
		printf("%d ", ptbl2->ent.two.pht_val[i]);
    }
    break;

  case DPredHybrid: 
    if (*no_pred == 2)
	printf("VHT MISS ");
    else { 
        /* printf("no_pred: %d\n", *no_pred); */
        assert(ptbl1);
        assert(ptbl2);

        if (pred->datapred.hybrid->config.hybrid.is_stride) {
		if (ptbl1->ent.hybrid.vht_val.state == Init)
			printf("Init ");
		else if (ptbl1->ent.hybrid.vht_val.state == Transient)
			printf("Transient ");
		else if (ptbl1->ent.hybrid.vht_val.state == Steady)
			printf("Steady ");		
	}

	for (i=0; i < pred->datapred.hybrid->config.hybrid.hist; i++)
		printf("%d ", ptbl1->ent.hybrid.vht_val.values[i]);
       	index = ptbl1->ent.hybrid.vht_val.vhp & 
		(pred->datapred.hybrid->config.hybrid.pht_size - 1);
	printf("vhp:%x - ", index);

        if (pred->datapred.hybrid->config.hybrid.xor) { /* xor */
  	  unsigned int pctemp;
	  pctemp = pcaddr;
	  switch (pred->datapred.hybrid->config.hybrid.xor) {
	  	case 6 : pctemp &= 0x3f; break;
	  	case 8 : pctemp &= 0xff; break;	  
	  	case 12 : pctemp &= 0xfff; break;
	  	case 16 : pctemp &= 0xffff; break;
	  	default : fatal("illegal xor bit size");
	  }
          index = (ptbl1->ent.hybrid.vht_val.vhp ^ pctemp) 
                & (pred->datapred.hybrid->config.hybrid.pht_size - 1);
	  printf("vhp^pc: %x - ", index);
        } /* xor */	
	if (!*no_pred && ptbl2->addr != pcaddr)
		printf("ALIAS(%x) ", ptbl2->addr);
	for (i=0; i < pred->datapred.hybrid->config.hybrid.hist; i++)
		printf("%d ", ptbl2->ent.hybrid.pht_val[i]);
    }
    break;

  default:
    panic("bogus data predictor class");
  } 
  printf("\n");
}



