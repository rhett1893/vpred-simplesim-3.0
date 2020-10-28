/*
 * sim-dpred.c - sample data value predictor simulator implementation
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

#include "misc.h"
#include "ss.h"
#include "regs.h"
#include "memory.h"
#include "loader.h"
#include "syscall.h"
#include "dlite.h"
#include "options.h"
#include "stats.h"
#include "dpred.h"
#include "sim.h"

/*
 * This file implements a branch predictor analyzer.
 */

/* data value predictor type {last|stride|2lev} */
static char *pred_type;

/* last predictor config (<ct_size> <counter_size> <vpt_size> <hist_size>) */
static int last_nelt = 4;
static int last_config[4] =
  { /* ct_size */1024, /* counter_size */2, /* vpt_size */4096, /* hist */1};

/* stride predictor config (<vht_size>) */
static int stride_nelt = 1;
static int stride_config[1] =
  { /* vht_size */4096};

/* two level predictor config (<vht_size><threshold><pht_size><hist><xor>) */
static int two_nelt = 5;
static int two_config[5] =
  { /* vht_size */4096, /* threhold */3, /* pht_size */4096, /* hist */4, /* xor */0};

/* hybrid predictor config (<vht_size><threshold><pht_size><hist><xor>) */
static int hybrid_nelt = 5;
static int hybrid_config[5] =
  { /* vht_size */4096, /* threhold */6, /* pht_size */4096, /* hist */4, /* xor */0};

static int sh_trace = 0;  /* show trace flag */

/* data predictor */
static struct dpred *pred;

/* track number of insn and refs */
static SS_COUNTER_TYPE sim_num_insn = 0;
static SS_COUNTER_TYPE sim_num_refs = 0;

/* total number of branches executed */
static SS_COUNTER_TYPE sim_num_branches = 0;

/* total number of insn touched to be predicted */
static SS_COUNTER_TYPE sim_num_datapred = 0;


/* register simulator-specific options */
void
sim_reg_options(struct opt_odb_t *odb)
{
  opt_reg_header(odb, 
"sim-dpred: This simulator implements a data value predictor analyzer.\n"
		 );

  /* data predictor options */
  opt_reg_note(odb,
"  Data value predictor configuration examples for last predictor:\n"
"    Configurations:   N, C, M, H\n"
"      N   # entries in CT(Classification Table)\n"
"      C   # bit width of a counter in CT\n"
"      M   # entries in VPT(Value Prediction Table)\n"
"      H   # of history depth\n" 
"      		this determines the # of data values in an entry of VPT\n"
               );

  opt_reg_note(odb,
"  Data value predictor configuration examples for stride predictor:\n"
"    Configurations:   N\n"
"      N   # entries in VHT(Value History Table)\n"
               );

  opt_reg_note(odb,
"  Data value predictor configuration examples for 2-level predictor:\n"
"    Configurations:   N, T, P, H, X\n"
"      N   # entries in first level VHT(Value History Table)\n"
"      T   # counter threshold value(3*,6,9)\n"
"      P   # entries in second level PHT(Pattern History Table)\n"
"      H   # of history depth(2,4*,6,8,10,12,14,16)\n" 
"      		this determines the # of data values in VHT and counters in PHT\n"
"      X   # of bits PC-address for xoring VHP for PHT index\n"
"	 	0*-no, 6, 8, 12, 16\n"
               );

  opt_reg_note(odb,
"  Data value predictor configuration examples for hybrid predictor:\n"
"    Configurations:   N, T, P, H, X\n"
"      N   # entries in first level VHT(Value History Table)\n"
"      T   # counter threshold value(3,6*,9)\n"
"      P   # entries in second level PHT(Pattern History Table)\n"
"      H   # of history depth(2,4*,6,8,10,12,14,16)\n" 
"      		this determines the # of data values in VHT and counters in PHT\n"
"      X   # of bits PC-address for xoring VHP for PHT index\n"
"	 	0*-no, 6, 8, 12, 16\n"
               );

  opt_reg_string(odb, "-dpred",
		 "data value predictor type {last|stride|2lev|hybrid}",
                 &pred_type, /* default */"last",
                 /* print */TRUE, /* format */NULL);

  opt_reg_int_list(odb, "-dpred:last",
		   "last  data value predictor config "
                   "(<ct_size> <counter_size> <vpt_size> <hist_size>)",
		   last_config, last_nelt, &last_nelt,
		   /* default */last_config,
		   /* print */TRUE, /* format */NULL, /* !accrue */FALSE);

  opt_reg_int_list(odb, "-dpred:stride",
		   "stride  data value predictor config "
                   "(<vht_size>)",
		   stride_config, stride_nelt, &stride_nelt,
		   /* default */stride_config,
		   /* print */TRUE, /* format */NULL, /* !accrue */FALSE);

  opt_reg_int_list(odb, "-dpred:2lev",
		   "two level data value predictor config "
                   "(<vht_size> <threshold> <pht_size> <hist><xor>)",
		   two_config, two_nelt, &two_nelt,
		   /* default */two_config,
		   /* print */TRUE, /* format */NULL, /* !accrue */FALSE);

  opt_reg_int_list(odb, "-dpred:hybrid",
		   "hybrid data value predictor config "
                   "(<vht_size> <threshold> <pht_size> <hist><xor>)",
		   hybrid_config, hybrid_nelt, &hybrid_nelt,
		   /* default */hybrid_config,
		   /* print */TRUE, /* format */NULL, /* !accrue */FALSE);

  opt_reg_int(odb, "-trace",
             	"show instruction trace {1 : show trace} ",
		&sh_trace, sh_trace,TRUE,NULL);
}


/* check simulator-specific option values */
void
sim_check_options(struct opt_odb_t *odb, int argc, char **argv)
{
  if (!mystricmp(pred_type, "last"))
    {
      if (last_nelt != 4)
	fatal("bad last predictor config (<ct_size> <counter_size> <vpt_size> <hist_size>)");
      if (last_config[1] != 2)
	fatal("only allowed for counter_size = 2 in this version");
      if (last_config[3] != 1)
	fatal("only allowed for hist_size = 1 in this version");

      pred = dpred_create(DPredLast,
			  /* last ct size */last_config[0],
			  /* last counter size */last_config[1],
			  /* last vpt size */last_config[2],
			  /* last history depth */last_config[3],
			  /* VHP xor flag */0);
    }
  else if (!mystricmp(pred_type, "stride"))
    {
      if (stride_nelt != 1)
	fatal("bad stride predictor config (<vht_size>)");

      pred = dpred_create(DPredStride,
			  /* stride vht size */stride_config[0],
			  0, 0, 0, 0);	
    }
  else if (!mystricmp(pred_type, "2lev"))
    {
      if (two_nelt != 5)
	fatal("bad 2-level predictor config (<vht_size><threshold><pht_size><hist><xor>)");
      pred = dpred_create(DPred2Level,
			  /* two vht size */two_config[0],
			  /* two threshold */two_config[1],
			  /* two pht_size */two_config[2],
			  /* two history depth */two_config[3],
			  /* two xor flag */two_config[4]);
    }
  else if (!mystricmp(pred_type, "hybrid"))
    {
      if (hybrid_nelt != 5)
	fatal("bad hybrid predictor config (<vht_size><threshold><pht_size><hist><xor>)");
      pred = dpred_create(DPredHybrid,
			  /* hybrid vht size */hybrid_config[0],
			  /* hybrid threshold */hybrid_config[1],
			  /* hybrid pht_size */hybrid_config[2],
			  /* hybrid history depth */hybrid_config[3],
			  /* hybrid xor flag */hybrid_config[4]);
    }
  else
    fatal("cannot parse predictor type `%s'", pred_type);

  pred->sh_trace = sh_trace;  /* set the flag of instruction trace */
}


/* register simulator-specific statistics */
void
sim_reg_stats(struct stat_sdb_t *sdb)
{
  stat_reg_counter(sdb, "sim_num_insn",
		   "total number of instructions executed",
		   &sim_num_insn, 0, NULL);
  stat_reg_counter(sdb, "sim_num_refs",
		   "total number of loads and stores executed",
		   &sim_num_refs, 0, NULL);
  stat_reg_int(sdb, "sim_elapsed_time",
	       "total simulation time in seconds",
	       (int *)&sim_elapsed_time, 0, NULL);
  stat_reg_formula(sdb, "sim_inst_rate",
		   "simulation speed (in insts/sec)",
		   "sim_num_insn / sim_elapsed_time", NULL);

  stat_reg_counter(sdb, "sim_num_branches",
                   "total number of branches executed",
                   &sim_num_branches, /* initial value */0, /* format */NULL);
  stat_reg_formula(sdb, "sim_IPB",
                   "instruction per branch",
                   "sim_num_insn / sim_num_branches", /* format */NULL);

  stat_reg_counter(sdb, "sim_num_datapred",
                   "total number of instructions touched to predict",
                   &sim_num_datapred, /* initial value */0, /* format */NULL);

  stat_reg_formula(sdb, "sim_pred_inst_rate",
                   "sim_num_datapred/sim_num_insn",
                   "sim_num_datapred / sim_num_insn", /* format */NULL);

  /* register predictor stats */
  if (pred)
    dpred_reg_stats(pred, sdb);
}

/* initialize the simulator */
void
sim_init(void)
{
  SS_INST_TYPE inst;

  sim_num_insn = 0;
  sim_num_refs = 0;

  regs_PC = ld_prog_entry;

  /* decode all instructions */
  {
    SS_ADDR_TYPE addr;

    if (OP_MAX > 255)
      fatal("cannot perform fast decoding, too many opcodes");

    debug("sim: decoding text segment...");
    for (addr=ld_text_base;
	 addr < (ld_text_base+ld_text_size);
	 addr += SS_INST_SIZE)
      {
	inst = __UNCHK_MEM_ACCESS(SS_INST_TYPE, addr);
	inst.a = SWAP_WORD(inst.a);
	inst.b = SWAP_WORD(inst.b);
	inst.a = (inst.a & ~0xff) | (unsigned int)SS_OP_ENUM(SS_OPCODE(inst));
	__UNCHK_MEM_ACCESS(SS_INST_TYPE, addr) = inst;
      }
  }

  /* initialize the DLite debugger */
  dlite_init(dlite_reg_obj, dlite_mem_obj, dlite_mstate_obj);
}

/* print simulator-specific configuration information */
void
sim_aux_config(FILE *stream)		/* output stream */
{
  /* nothing currently */
}

/* dump simulator-specific auxiliary simulator statistics */
void
sim_aux_stats(FILE *stream)		/* output stream */
{
  /* nada */
}

/* un-initialize simulator-specific state */
void
sim_uninit(void)
{
  /* nada */
}




/*
 * configure the execution engine
 */

/*
 * precise architected register accessors
 */

/* next program counter */
#define SET_NPC(EXPR)		(next_PC = (EXPR))

/* target program counter */
#undef  SET_TPC
#define SET_TPC(EXPR)		(target_PC = (EXPR))

/* current program counter */
#define CPC			(regs_PC)

/* general purpose registers */
#define GPR(N)			(regs_R[N])
#define SET_GPR(N,EXPR)		(regs_R[N] = (EXPR))

/* floating point registers, L->word, F->single-prec, D->double-prec */
#define FPR_L(N)		(regs_F.l[(N)])
#define SET_FPR_L(N,EXPR)	(regs_F.l[(N)] = (EXPR))
#define FPR_F(N)		(regs_F.f[(N)])
#define SET_FPR_F(N,EXPR)	(regs_F.f[(N)] = (EXPR))
#define FPR_D(N)		(regs_F.d[(N) >> 1])
#define SET_FPR_D(N,EXPR)	(regs_F.d[(N) >> 1] = (EXPR))

/* miscellaneous register accessors */
#define SET_HI(EXPR)		(regs_HI = (EXPR))
#define HI			(regs_HI)
#define SET_LO(EXPR)		(regs_LO = (EXPR))
#define LO			(regs_LO)
#define FCC			(regs_FCC)
#define SET_FCC(EXPR)		(regs_FCC = (EXPR))

/* additional definition on sim-bpred */
#define DGPR(N)			(N)
#define DGPR_D(N)		((N)&~1)
#define DCGPR(N)		(SS_COMP_OP != SS_COMP_NOP ? (N) : 0)
#define DFPR_D(N)		(((N)+32)&~1)
#define DFPR_F(N)		(((N)+32)&~1)
#define DFPR_L(N)		(((N)+32)&~1)


#define DNA			(0)
#define DHI			(0+32+32)
#define DLO			(1+32+32)
#define DFCC			(2+32+32)

/* precise architected memory state help functions */
#define __READ_WORD(DST_T, SRC_T, SRC)					\
  ((unsigned int)((DST_T)(SRC_T)MEM_READ_WORD(addr = (SRC))))

#define __READ_HALF(DST_T, SRC_T, SRC)					\
  ((unsigned int)((DST_T)(SRC_T)MEM_READ_HALF(addr = (SRC))))

#define __READ_BYTE(DST_T, SRC_T, SRC)					\
  ((unsigned int)((DST_T)(SRC_T)MEM_READ_BYTE(addr = (SRC))))

/* precise architected memory state accessor macros */
#define READ_WORD(SRC)							\
  __READ_WORD(unsigned int, unsigned int, (SRC))

#define READ_UNSIGNED_HALF(SRC)						\
  __READ_HALF(unsigned int, unsigned short, (SRC))

#define READ_SIGNED_HALF(SRC)						\
  __READ_HALF(signed int, signed short, (SRC))

#define READ_UNSIGNED_BYTE(SRC)						\
  __READ_BYTE(unsigned int, unsigned char, (SRC))

#define READ_SIGNED_BYTE(SRC)						\
  __READ_BYTE(signed int, signed char, (SRC))

#define WRITE_WORD(SRC, DST)						\
  (MEM_WRITE_WORD(addr = (DST), (unsigned int)(SRC)))

#define WRITE_HALF(SRC, DST)						\
  (MEM_WRITE_HALF(addr = (DST), (unsigned short)(unsigned int)(SRC)))

#define WRITE_BYTE(SRC, DST)						\
  (MEM_WRITE_BYTE(addr = (DST), (unsigned char)(unsigned int)(SRC)))

/* system call handler macro */
#define SYSCALL(INST)		(ss_syscall(mem_access, INST))

/* instantiate the helper functions in the '.def' file */
#define DEFINST(OP,MSK,NAME,OPFORM,RES,CLASS,O1,O2,I1,I2,I3,EXPR)
#define DEFLINK(OP,MSK,NAME,MASK,SHIFT)
#define CONNECT(OP)
#define IMPL
#include "ss.def"
#undef DEFINST
#undef DEFLINK
#undef CONNECT
#undef IMPL

/* show instruction trace */ 
void
trace_inst(SS_ADDR_TYPE pc, char *op, char *oprnd, char is_pred)
{
  if (pred->sh_trace)
  	printf("%x %s %s - %s\n", pc, op, oprnd, 
          	is_pred ? "pred_inst" : "no_pred_inst");
}

/* start simulation, program loaded, processor precise state initialized */
void
sim_main(void)
{
  SS_INST_TYPE inst;
  register SS_WORD_TYPE res_DATA;		/* the resolved data executed */
  register SS_ADDR_TYPE addr, next_PC, target_PC;
  enum ss_opcode op;
  register int is_write;

  fprintf(stderr, "sim: ** starting functional simulation **\n");

  /* set up initial default next PC */
  next_PC = regs_PC + SS_INST_SIZE;

  /* check for DLite debugger entry condition */
  if (dlite_check_break(regs_PC, /* no access */0, /* addr */0, 0, 0))
    dlite_main(regs_PC - SS_INST_SIZE, regs_PC, sim_num_insn);


  if (!mystricmp(pred_type, "last"))
  	dpred_data_config(pred->datapred.last, "last", stdout);
  else if (!mystricmp(pred_type, "stride"))
  	dpred_data_config(pred->datapred.stride, "stride", stdout);
  else if (!mystricmp(pred_type, "2lev"))
  	dpred_data_config(pred->datapred.two, "2lev", stdout);
  else if (!mystricmp(pred_type, "hybrid"))
  	dpred_data_config(pred->datapred.hybrid, "hybrid", stdout);
  else
    fatal("cannot parse predictor type `%s' in sim_main()", pred_type);


  while (TRUE)
    {
      /* maintain $r0 semantics */
      regs_R[0] = 0;

      /* keep an instruction count */
      sim_num_insn++;

      /* get the next instruction to execute */
      inst = __UNCHK_MEM_ACCESS(SS_INST_TYPE, regs_PC);

      /* set default reference address and access mode */
      addr = 0; is_write = FALSE;

      /* decode the instruction */
      op = SS_OPCODE(inst);
      switch (op)
	{
#define DEFINST(OP,MSK,NAME,OPFORM,RES,FLAGS,O1,O2,I1,I2,I3,EXPR)	\
	case OP:                                                        \
          EXPR;                                                         \
          if (is_PRED) 							\
           { if (O1 == DHI) res_DATA = HI;				\
             else if (O1 == DLO) res_DATA = LO;				\
             else res_DATA = GPR(O1);                              	\
	     trace_inst(regs_PC,NAME,OPFORM,1);				\
	   }								\
           else								\
	     trace_inst(regs_PC,NAME,OPFORM,0);				\
          break;
#define DEFLINK(OP,MSK,NAME,MASK,SHIFT)                                 \
        case OP:                                                        \
          panic("attempted to execute a linking opcode");
#define CONNECT(OP)
#include "ss.def"
#undef DEFINST
#undef DEFLINK
#undef CONNECT
	default:
	  panic("bogus opcode");
      }

      if (SS_OP_FLAGS(op) & F_MEM)
	{
	  sim_num_refs++;
	  if (SS_OP_FLAGS(op) & F_STORE)
	    is_write = TRUE;
	}

      if (SS_OP_FLAGS(op) & F_CTRL)
	  sim_num_branches++;

      /* if this is a integer register write instruction, make data prediction */
      if (is_PRED)
	{
	  struct dpred_ent *ptbl1, *ptbl2; /* pointers of the table entry */
	  SS_WORD_TYPE pred_DATA;
	  char no_pred;	/* flag for no prediction
			   0 : pred, 1 : no pred by counter(last,2lev) or by state(stride), 
			   2 : no pred by CT miss(last) or by VHT miss(stride,2lev), 
     			   3 : no pred by VPT miss(last),
			   4 : no pred by others */

	  sim_num_datapred++;

	  if (pred)
	    {
	      no_pred = 0;
 
	      /* get the next predicted data value */
	      pred_DATA = dpred_lookup(pred,
				     /* PC addr */regs_PC,
				     /* opcode */op,
				     /* no pred flag */&no_pred,
				     /* table1 pointer */&ptbl1,
				     /* table2 pointer */&ptbl2);

	      /* print the trace of the table entry */
	      if (pred->sh_trace)
	      	dpred_trace(pred,
		         /* PC addr */regs_PC,
			 /* predicted data*/pred_DATA,
		         /* resolved data */res_DATA,
            	         /* no pred flag */&no_pred,
			 /* table1 pointer */ptbl1,
			 /* table2 pointer */ptbl2);
                 
              dpred_update(pred,
			   /* PC addr */regs_PC,
			   /* resolved data */res_DATA,
		           /* no pred flag */no_pred,
			   /* correct pred? */pred_DATA == res_DATA,
			   /* opcode */op);
	    }
	}

      /* check for DLite debugger entry condition */
      if (dlite_check_break(next_PC,
			    is_write ? ACCESS_WRITE : ACCESS_READ,
			    addr, sim_num_insn, sim_num_insn))
	dlite_main(regs_PC, next_PC, sim_num_insn);

      /* go to the next instruction */
      regs_PC = next_PC;
      next_PC += SS_INST_SIZE;
    }
}
