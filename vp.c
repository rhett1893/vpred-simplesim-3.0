/************************************************/
/* vp.c - stride Value Predictor.		*/
/* Technion - Israel Institute of Technology.	*/
/* Department of EE.				*/
/* F. Gabbay, e-mail: fredg@tx.technion.ac.il 	*/
/* Y. Rivlin & O. Pozner, e-mail: s3751427@psl  */
/************************************************/

#include <stdio.h>
//#include "misc.h"
#include "machine.h"
// #include "ss.h"
#include "memory.h"

//#include "define.h"
#include "vp.h"
// #include "ss.h"
int common_ref=0;
int misses_vpred=0;
/* decleration of extern variables (may be moved to the sim_main) */
int                       en_lookup_stride; /* enable X mechanism */
tick_t                    sim_cycle;        /* simulation time */
int                       use_stride;       /* if (use_stride == 1) then
						use stride predictor, else
						use last value predictor */
int                       use_fsm;          /* use Final state machine flag */

int                       vp_replace;       /* Value prediction replace policy*/
int                       start_fsm;        /* FSM intial state */
int                       hash_no_i1;       /* Single precision set number */
int                       hash_no_i2;       /* Double precision set number */
int                       hash_asso_i1;     /* Single precision set 
					       associativity */
int                       hash_asso_i2;     /* Double precision set 
					       associativity */
int predicted_ok=0;     /* total number of correct predictions */
int lookup_accessed=0;
int update_accessed=0;
int allocate_accessed=0;
extern int use_vp;           /* use value prediction */
static Hash_ptr_i1       *hash_table_i1;    /* pointers to VP tables */
static Hash_ptr_i2       *hash_table_i2;


/*allocates the memory needed for the Value Prediction Tables*/ 
void create_stride()
{
	register i;
	hash_no_i1 = MEM_PTAB_SIZE;
  hash_asso_i1 = 4;
	if(use_vp==0) /* check if value prediction is used */
	  return;	
	if((hash_table_i1=(Hash_ptr_i1 *) calloc (hash_no_i1,sizeof(Hash_ptr_i1)))==0){
		fprintf(stderr,"out of memory. terminating program.\n");  /*calloc error*/
		exit(-1);
	}
	for(i=0;i<hash_no_i1;i++)
		if((hash_table_i1[i]=(Hash_i1 *) calloc (hash_asso_i1,sizeof(Hash_i1)))==0){
			fprintf(stderr,"out of memory. terminating program.\n");  /*calloc error*/
			exit(-1);
		}
	if((hash_table_i2=(Hash_ptr_i2 *) calloc (hash_no_i2,sizeof(Hash_ptr_i2)))==0){
		fprintf(stderr,"out of memory. terminating program.\n");  /*calloc error*/
		exit(-1);
	}
	for(i=0;i<hash_no_i2;i++)
		if((hash_table_i2[i]=(Hash_i2 *) calloc (hash_asso_i2,sizeof(Hash_i2)))==0){
			fprintf(stderr,"out of memory. terminating program.\n");  /*calloc error*/
			exit(-1);
		}

} /* end of create_stride */   

/* update an instruction in the Vp table (update output outcomes)*/ 
int update(
 md_addr_t instpc, /* The instruction address */
 md_inst_t inst,   /* The instruction opcode and registers */
 VAL_TAG_TYPE pred_val, /* The given (by lookup) predicted value and flags */
 VAL_TAG_TYPE calc_val, /* The actual (calc) value and flags */
 int my_ref             /* The Time of the update */  )
{
  update_accessed++;
  unsigned int index;  /* set index */
  unsigned int pc;     /* shifted inst address */
  int val[2];          /* instruction output value single/double int */
  register i;
  char flag[2];        /* single/double precision correct prediction flags */
  char found=0;        /* inst "found" status */
  enum md_opcode op;   /* decoded opcode enum */
 
  if(use_vp==0) /* check if using VP */
     return(MISS);	
  pc = (unsigned int)(instpc>>INST_OFFSET); 
  flag[0]=0;
  flag[1]=0; 
  MD_SET_OPCODE(op, inst);  /* recover opcode from inst */  
  //if(SS_OP_FLAGS(op)&F_SINGLE_P){ 
    /* single precision output */
     val[0]=calc_val.value.single_p;
     index = (unsigned int)(pc % hash_no_i1); /* find table entry */
     for(i=0;i<hash_asso_i1;i++)  
       /* find inst in set */
         if((hash_table_i1[index][i].pc==pc)&&(hash_table_i1[index][i].valid)){
	    found=1;
	    break;
	 }			
     if(!found)                   /*miss*/
	return(MISS);
     else{  
        if(val[0]==pred_val.value.single_p){ 
              /*hit*/
	   predicted_ok++;       /* increase total predicated ok counter */
	   /* update the classification fsm state */
	   if(hash_table_i1[index][i].ps<3) 
	      hash_table_i1[index][i].ps++;
	   flag[0]=1;                      
	   hash_table_i1[index][i].pred_correct++; /* statistics */
	 }
         else{ /* hit fault: found, but wrong prediction was given */
	   /* update the classification fsm state */
	    if(hash_table_i1[index][i].ps>0)  
	      hash_table_i1[index][i].ps--;
	    
	 }

       hash_table_i1[index][i].accessed++;    /* statistics */

       if (hash_table_i1[index][i].x<1)
	   panic("invalid X value (lookup stride)");
       
       /* X is decreased only if:
	   A. X was increased in lookup (there was a lookup hit).
	   B. X > 1.   */
       if ((hash_table_i1[index][i].x>1)&&(pred_val.fsm_pred!=MISS)) 
	 hash_table_i1[index][i].x--;     /* lookup X mechanism update */

       /* stride, last_val, last_ref are updated only in order */

       if (hash_table_i1[index][i].last_ref<my_ref) {
         /* update stride value */ 
         hash_table_i1[index][i].stride=val[0]-hash_table_i1[index][i].last_val;
         hash_table_i1[index][i].last_ref=my_ref;
         hash_table_i1[index][i].last_val=val[0]; /* update last value */
       }
       /* else:
	    there was already an update of a command that appears later 
	    in the code (sequentally) - stride, last_val, last_ref are not
	    updated */

       if(flag[0]==1)
	 return(HIT);
       else return(HITFAULT);
     }
  //}
  /*
  Not Using Double Precision Calculations
  */
 //  else if(SS_OP_FLAGS(op)&F_DOUBLE_P){
 //    /* double precision output */
 //     val[0]=calc_val.value.double_p[0];
 //     val[1]=calc_val.value.double_p[1];
 //     index=(unsigned int)(pc % hash_no_i2);
 //     for(i=0;i<hash_asso_i2;i++)
	//    /* find inst in set */
	// if((hash_table_i2[index][i].pc==pc)&&(hash_table_i2[index][i].valid)){
	//    found=1;
	//    break;
	// }
 //     if(!found)                   /*miss*/
	// return(MISS);
 //     else{               
	// for(j=0;j<2;j++) 
 //         if(val[j]==pred_val.value.double_p[j]){
 //            /*hit*/
	//       predicted_ok++; /* statistics */
	//       flag[j]=1;
	//       hash_table_i2[index][i].pred_correct[j]++; /* statistics */
	//    }

 //        hash_table_i2[index][i].accessed++; /* statistics */

	// if (hash_table_i2[index][i].x<1)
 // 	   panic("invalid X value (lookup stride)"); 	

 //        if ((hash_table_i2[index][i].x>1)&&(pred_val.fsm_pred!=MISS)) 
	//    hash_table_i2[index][i].x--;    /* update lookup X mechanism*/

 //        if (hash_table_i2[index][i].last_ref<my_ref){
	//   /* update stride & last_val */
 //          hash_table_i2[index][i].stride[0]=
	//     val[0] - hash_table_i2[index][i].last_val[0];
 //          hash_table_i2[index][i].last_val[0]=val[0];
 //          hash_table_i2[index][i].stride[1]=
	//     val[1] - hash_table_i2[index][i].last_val[1];
 //          hash_table_i2[index][i].last_val[1]=val[1];

	//   /* update last_ref */
	//   hash_table_i2[index][i].last_ref = my_ref;
 //        }  


	//    /* update the classification fsm state */
 //        if((flag[0]==1)&&(flag[1]==1)){ 
	//    if(hash_table_i2[index][i].ps<3)
	//       hash_table_i2[index][i].ps++;
	//    return(HIT);
	// }
	// else{ /* hit fault */
 //           if(hash_table_i2[index][i].ps>0)
	//       hash_table_i2[index][i].ps--;
	//    return(HITFAULT);
 //     	}
 //     }
 //  }
 // else return(MISS);
}

/* recover from a pipe flush - decrease the X value */
int lookup_undo(
 md_addr_t instpc,/* The instruction address */
 md_inst_t inst,  /* The instruction opcode and registers */
 VAL_TAG_TYPE pred_val  /* The instruction predicted value and flags */ ) 
{
  unsigned int index;
  unsigned int pc;
  register i;
  char found=0;
  enum md_opcode op; /* decoded opcode enum */
 
  if(use_vp==0)
     return(MISS);	
  pc = (unsigned int)(instpc>>INST_OFFSET); 
  MD_SET_OPCODE(op, inst); 
  //if(SS_OP_FLAGS(op)&F_SINGLE_P){
     index = (unsigned int)(pc % hash_no_i1);
     for(i=0;i<hash_asso_i1;i++)
         if((hash_table_i1[index][i].pc==pc)&&(hash_table_i1[index][i].valid)){
	    found=1;
	    break;
	 }			
     if(!found)                   /*miss*/
	return(MISS);
     else{  
       if (hash_table_i1[index][i].x<1)
	   panic("invalid X value (lookup stride)");

       if ((hash_table_i1[index][i].x>1)&&(pred_val.fsm_pred!=MISS)) 
	   hash_table_i1[index][i].x--;
       return(HIT);
     }
 // }
 //  else if(SS_OP_FLAGS(op)&F_DOUBLE_P){
 //     index=(unsigned int)(pc % hash_no_i2);
 //     for(i=0;i<hash_asso_i2;i++)
	// if((hash_table_i2[index][i].pc==pc)&&(hash_table_i2[index][i].valid)){
	//    found=1;
	//    break;
	// }
 //     if(!found)                   /*miss*/
	// return(MISS);
 //     else{               
	// if (hash_table_i2[index][i].x<1)
 // 	   panic("invalid X value (lookup stride)"); 
        
 //        if ((hash_table_i2[index][i].x>1)&&(pred_val.fsm_pred!=MISS)) 
	//     hash_table_i2[index][i].x--;
	// return(HIT);
 //     }
 //  }
  //else return(MISS);
}

/* find a replacment candidate. Replacement policy: 
   LRU */ 
unsigned int find_new_place_lru(
int index, /* the table entry in which the record should be place */ 
char data_type /* single or double presicion */)
{
  unsigned int i;
  int min_ref;
  unsigned int new_place;

  switch(data_type){
	  case SINGLE:
	                 /* init min_ref */
			 min_ref=hash_table_i1[index][0].last_ref; 
			 new_place=0;
			 for(i=0;i<hash_asso_i1;i++){
				 if(hash_table_i1[index][i].valid == 0)
				 return(i);
				 else if(hash_table_i1[index][i].last_ref<min_ref){
				 min_ref=hash_table_i1[index][i].last_ref;
				 new_place=i;
				 }
			 }
			 return(new_place);
			 break;
	  // case DOUBLE:
	  //                /* init min_ref */
			//  min_ref=hash_table_i2[index][0].last_ref;
			//  new_place=0;
			//  for(i=0;i<hash_asso_i2;i++){
			// 	 if(hash_table_i2[index][i].valid == 0)
			// 	 return(i);
			// 	 else if(hash_table_i2[index][i].last_ref<min_ref){
			// 	 min_ref=hash_table_i2[index][i].last_ref;
			// 	 new_place=i;
			// 	 }
			//  }
			//  return(new_place);
			//  break;
          default :      panic("unknown data type\n");
      }
}

/* find a replacment candidate. Replacment policy:
   find all the entries with the (same) lowest fsm state (these are less
   predictable insts). from all these replace the LRU entry. */ 
unsigned int find_new_place_fsm(int index, char data_type)
{
  unsigned int i;
  char min_ps;
  int min_ref;
  unsigned int new_place;

  switch(data_type){
	  case SINGLE:
	                 min_ps=hash_table_i1[index][0].ps;
			 min_ref=hash_table_i1[index][0].last_ref;
			 new_place=0;
			 for(i=0;i<hash_asso_i1;i++){
			     if(hash_table_i1[index][i].valid == 0)
				return(i);
			     else if(hash_table_i1[index][i].ps<min_ps){
			         min_ps=hash_table_i1[index][i].ps;
				 min_ref=hash_table_i1[index][i].last_ref;
			         new_place=i;
			     }
			     else if((hash_table_i1[index][i].ps==min_ps)&&
				    (hash_table_i1[index][i].last_ref<min_ref)){
				 min_ref=hash_table_i1[index][i].last_ref;
				 new_place=i;
			     }
			 }
			 return(new_place);
			 break;
	  // case DOUBLE:
			//  min_ref=hash_table_i2[index][0].last_ref;
			//  new_place=0;
			//  for(i=0;i<hash_asso_i2;i++){
			//      if(hash_table_i2[index][i].valid == 0)
			// 	 return(i);
			//      else if(hash_table_i2[index][i].ps<min_ps){
			//          min_ps=hash_table_i2[index][i].ps;
			// 	 min_ref=hash_table_i2[index][i].last_ref;
			//          new_place=i;
			//      }
			//      else if((hash_table_i2[index][i].ps==min_ps)&&
			// 	    (hash_table_i2[index][i].last_ref<min_ref)){
			// 	 min_ref=hash_table_i2[index][i].last_ref;
			// 	 new_place=i;
			//      }
			//  }
			//  return(new_place);
			//  break;
          default :      panic("unknown data type\n");
      }
}


/* allocates a new record */
void allocate(md_addr_t pred_PC, /* allocated inst's address */
md_inst_t inst, /* allocated inst's opcode and registers */
VAL_TAG_TYPE calc_val,  /* allocated inst's last output value and flags */
int my_ref /* allocated inst's fetch/decode time stamp */)
{
allocate_accessed++;
unsigned int index;    /* table index */
unsigned int pc;       /* shifted pc address */
int val[2];            /* the instruction output(s) value */
unsigned int replace;  /* the record in which the pc is allocated */
//register j;
char data_type;        /* instruction precision */
enum md_opcode op;     /* instruction opcode */
  
/* using Value Prediction ? */
if(use_vp==0)     
   return;	
MD_SET_OPCODE(op, inst);
pc = (unsigned int)(pred_PC)>>INST_OFFSET;

// if(SS_OP_FLAGS(op)&F_SINGLE_P){  
  /* single precision */
   val[0]=calc_val.value.single_p;
   data_type=SINGLE;
   index=(unsigned int)(pc % hash_no_i1);

   /* find allocation entry */
   replace=(vp_replace ? find_new_place_lru(index,data_type)
	               : find_new_place_fsm(index,data_type));

   /* update the entry fields */
   hash_table_i1[index][replace].pc = pc;
   hash_table_i1[index][replace].ps = start_fsm; 
   hash_table_i1[index][replace].x = 1;
   hash_table_i1[index][replace].last_val = val[0];
   hash_table_i1[index][replace].stride = 0;
   hash_table_i1[index][replace].accessed = 1;
   hash_table_i1[index][replace].pred_correct = 0;
   hash_table_i1[index][replace].last_alloc = (int) sim_cycle;
   hash_table_i1[index][replace].last_ref = (int) my_ref;
   hash_table_i1[index][replace].valid=1;
//} 
// else if(SS_OP_FLAGS(op)&F_DOUBLE_P){
//    val[0]=calc_val.value.double_p[0]; 
//    val[1]=calc_val.value.double_p[1]; 
//    index=(unsigned int)(pc % hash_no_i2);
//    data_type=DOUBLE;

//    /* find allocation entry */
//    replace=(vp_replace ? find_new_place_lru(index,data_type)
// 	               : find_new_place_fsm(index,data_type));

//    /* update the entry fields */
//    hash_table_i2[index][replace].pc = pc;
//    hash_table_i2[index][replace].ps = start_fsm;
//    hash_table_i2[index][replace].x = 1;
//    hash_table_i2[index][replace].accessed = 1;
//    hash_table_i2[index][replace].last_alloc = (int) sim_cycle;
//    hash_table_i2[index][replace].last_ref = (int) my_ref;
//    hash_table_i2[index][replace].valid = 1;
//    for(j=0;j<2;j++){
//       hash_table_i2[index][replace].last_val[j] = val[j];
//       hash_table_i2[index][replace].stride[j] = 0;
//       hash_table_i2[index][replace].pred_correct[j] = 0;
//    }
// }
}
/* lookup for inst in the VP table and return the predicted output values
    */
void lookup(md_addr_t pred_PC,/* The instruction address */
md_inst_t inst, /* The instruction calculated value and flags */
VAL_TAG_TYPE *pred_val,   /* The instruction predicted value and flags */
struct mem_t *mem  /* The memory block to be accessed for reading the value */ ) 
{
  lookup_accessed++;
  unsigned int index; /* VP table index*/
  unsigned int pc;    /* instruction shifted address */
  int val[2];
  int data;
  register i;
  char X;             /* X (lookup stride) mechanisim value */
  // char flag[2];      
  char found=0;
  enum md_opcode op;  /* decoded opcode enum */

  pred_val->fsm_pred=MISS;  
  VAL_TAG_TYPE calc_val;
  /* using Value Prediction ? */
  if(use_vp==0)
     return;
  pc = (unsigned int)(pred_PC>>INST_OFFSET); /* shifting the instruction
                                                address */ 
   mem_access(mem,Read,pred_PC,&data,sizeof(word_t));
   calc_val.value.single_p = data;
   calc_val.fsm_pred = MISS;
  // flag[0]=0;  clear flags 
  // flag[1]=0; 
  MD_SET_OPCODE(op, inst); /* recover op code from inst */

  //if(SS_OP_FLAGS(op)&F_SINGLE_P){ 
    /* single precision */
     index = (unsigned int)(pc % hash_no_i1);
     for(i=0;i<hash_asso_i1;i++)
       /* find entry in the single precision table */
         if((hash_table_i1[index][i].pc==pc)&&(hash_table_i1[index][i].valid)){
	    found=1;
	    break;
	 }			
     if(!found)                   /*miss*/
	{  pred_val->fsm_pred=MISS;
      allocate(pred_PC,inst,calc_val,common_ref++);
  }
     else{                        
       	X=hash_table_i1[index][i].x; /* set  the current lookup X value */
        hash_table_i1[index][i].x++; /* increase X value */
        /* calculate the return predicted output value and flags */
	val[0]=((en_lookup_stride==0)?hash_table_i1[index][i].stride:X*hash_table_i1[index][i].stride);
	pred_val->value.single_p=(hash_table_i1[index][i].last_val+((use_stride==0)?0:val[0]));	/*hit*/
	/* fsm classification:
	   fsm_pred == 1 => go with prediction
	   fsm_pred == 0 => don't use prediction */
        pred_val->fsm_pred=(use_fsm?hash_table_i1[index][i].ps>>1:HIT); 
        update(pred_PC,inst,*pred_val,calc_val,common_ref++);
     }
  //}
  // else{
  //    if(SS_OP_FLAGS(op)&F_DOUBLE_P){	
  //      /* double presicion */
  //       index = (unsigned int)(pc % hash_no_i2);
  //       /* find entry in the double precision table */
  //       for(i=0;i<hash_asso_i2;i++)
	 //   if((hash_table_i2[index][i].pc==pc)&&(hash_table_i2[index][i].valid)){
	 //      found=1;
	 //      break;
  //  	   }
  //       if(!found){                   /*miss*/
  //          pred_val->fsm_pred=MISS;
  //       }
  //       else{
  //          X=hash_table_i2[index][i].x;  /* set  the current lookuo X value */
  //  	   hash_table_i2[index][i].x++; /* update X for future use */
       
	 //     calculate and return predicted output value and flag 
  // 	   val[0]=((en_lookup_stride==0)?hash_table_i2[index][i].stride[0]:X*hash_table_i2[index][i].stride[0]);
  // 	   val[1]=((en_lookup_stride==0)?hash_table_i2[index][i].stride[1]:X*hash_table_i2[index][i].stride[1]);
  //          pred_val->value.double_p[0]=(hash_table_i2[index][i].last_val[0]+((use_stride==0)?0:val[0]));      /*hit*/
  //          pred_val->value.double_p[1]=(hash_table_i2[index][i].last_val[1]+((use_stride==0)?0:val[1]));     /*hit*/
  //          pred_val->fsm_pred=(use_fsm?hash_table_i2[index][i].ps>>1:HIT); 
  //       }
  //    }
  //    else return; /* value prediction for other instructions is not used*/

  // /* NOTE: the F_SINGLE_P and F_DOUBLE_P flags should be added to the 
  //          instruction set code (ss.def) for all instructions that you
	 //   want to be Value Predicted */ 
  // }
}