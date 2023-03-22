
/*  A Bison parser, made from ginsh_parser.yy
    by GNU Bison version 1.28  */

#define YYBISON 1  /* Identify Bison output.  */

#define	T_NUMBER	257
#define	T_SYMBOL	258
#define	T_LITERAL	259
#define	T_DIGITS	260
#define	T_QUOTE	261
#define	T_QUOTE2	262
#define	T_QUOTE3	263
#define	T_EQUAL	264
#define	T_NOTEQ	265
#define	T_LESSEQ	266
#define	T_GREATEREQ	267
#define	T_QUIT	268
#define	T_WARRANTY	269
#define	T_PRINT	270
#define	T_IPRINT	271
#define	T_TIME	272
#define	T_XYZZY	273
#define	T_INVENTORY	274
#define	T_LOOK	275
#define	T_SCORE	276
#define	NEG	277

#line 29 "ginsh_parser.yy"

#include "config.h"

#include <sys/resource.h>

#if HAVE_UNISTD_H
#include <sys/types.h>
#include <unistd.h>
#endif

#include <stdexcept>

#include "ginsh.h"

#define YYERROR_VERBOSE 1

// Original readline settings
static int orig_completion_append_character;
#if (GINAC_RL_VERSION_MAJOR < 4) || (GINAC_RL_VERSION_MAJOR == 4 && GINAC_RL_VERSION_MINOR < 2)
static char *orig_basic_word_break_characters;
#else
static const char *orig_basic_word_break_characters;
#endif

// Expression stack for ", "" and """
static void push(const ex &e);
static ex exstack[3];

// Start and end time for the time() function
static struct rusage start_time, end_time;

// Table of functions (a multimap, because one function may appear with different
// numbers of parameters)
typedef ex (*fcnp)(const exprseq &e);
typedef ex (*fcnp2)(const exprseq &e, int serial);

struct fcn_desc {
	fcn_desc() : p(NULL), num_params(0) {}
	fcn_desc(fcnp func, int num) : p(func), num_params(num), is_ginac(false) {}
	fcn_desc(fcnp2 func, int num, int ser) : p((fcnp)func), num_params(num), is_ginac(true), serial(ser) {}

	fcnp p;		// Pointer to function
	int num_params;	// Number of parameters (0 = arbitrary)
	bool is_ginac;	// Flag: function is GiNaC function
	int serial;	// GiNaC function serial number (if is_ginac == true)
};

typedef multimap<string, fcn_desc> fcn_tab;
static fcn_tab fcns;

static fcn_tab::const_iterator find_function(const ex &sym, int req_params);

// Table to map help topics to help strings
typedef multimap<string, string> help_tab;
static help_tab help;

static void insert_fcn_help(const char *name, const char *str);
static void print_help(const string &topic);
static void print_help_topics(void);
#ifndef YYSTYPE
#define YYSTYPE int
#endif
#include <stdio.h>

#ifndef __cplusplus
#ifndef __STDC__
#define const
#endif
#endif



#define	YYFINAL		106
#define	YYFLAG		-32768
#define	YYNTBASE	45

#define YYTRANSLATE(x) ((unsigned)(x) <= 277 ? yytranslate[x] : 54)

static const char yytranslate[] = {     0,
     2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
     2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
     2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
     2,     2,    33,     2,     2,     2,    30,     2,    39,    36,
    37,    28,    26,    44,    27,     2,    29,     2,     2,     2,
     2,     2,     2,     2,     2,     2,     2,    35,    34,    24,
    23,    25,    38,     2,     2,     2,     2,     2,     2,     2,
     2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
     2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
    42,     2,    43,    32,     2,     2,     2,     2,     2,     2,
     2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
     2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
     2,     2,    40,     2,    41,     2,     2,     2,     2,     2,
     2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
     2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
     2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
     2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
     2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
     2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
     2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
     2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
     2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
     2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
     2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
     2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
     2,     2,     2,     2,     2,     1,     3,     4,     5,     6,
     7,     8,     9,    10,    11,    12,    13,    14,    15,    16,
    17,    18,    19,    20,    21,    22,    31
};

#if YYDEBUG != 0
static const short yyprhs[] = {     0,
     0,     1,     4,     6,     9,    12,    18,    24,    27,    30,
    33,    35,    37,    39,    41,    43,    45,    46,    52,    55,
    58,    60,    62,    66,    68,    70,    72,    74,    76,    81,
    85,    89,    93,    97,   101,   105,   109,   113,   117,   121,
   125,   129,   132,   135,   139,   142,   146,   150,   154,   156,
   160,   161,   163,   165,   169,   173,   179,   181
};

static const short yyrhs[] = {    -1,
    45,    46,     0,    34,     0,    48,    34,     0,    48,    35,
     0,    16,    36,    48,    37,    34,     0,    17,    36,    48,
    37,    34,     0,    38,     4,     0,    38,    18,     0,    38,
    38,     0,    14,     0,    15,     0,    19,     0,    20,     0,
    21,     0,    22,     0,     0,    18,    47,    36,    48,    37,
     0,     1,    34,     0,     1,    35,     0,     3,     0,     4,
     0,    39,     4,    39,     0,     5,     0,     6,     0,     7,
     0,     8,     0,     9,     0,     4,    36,    49,    37,     0,
     6,    23,     3,     0,     4,    23,    48,     0,    48,    10,
    48,     0,    48,    11,    48,     0,    48,    24,    48,     0,
    48,    12,    48,     0,    48,    25,    48,     0,    48,    13,
    48,     0,    48,    26,    48,     0,    48,    27,    48,     0,
    48,    28,    48,     0,    48,    29,    48,     0,    27,    48,
     0,    26,    48,     0,    48,    32,    48,     0,    48,    33,
     0,    36,    48,    37,     0,    40,    50,    41,     0,    42,
    52,    43,     0,    48,     0,    49,    44,    48,     0,     0,
    51,     0,    48,     0,    51,    44,    48,     0,    42,    53,
    43,     0,    52,    44,    42,    53,    43,     0,    48,     0,
    53,    44,    48,     0
};

#endif

#if YYDEBUG != 0
static const short yyrline[] = { 0,
   115,   116,   119,   120,   129,   137,   145,   159,   160,   161,
   162,   163,   176,   177,   178,   179,   184,   184,   191,   192,
   195,   196,   197,   198,   199,   200,   201,   202,   203,   211,
   212,   213,   214,   215,   216,   217,   218,   219,   220,   221,
   222,   223,   224,   225,   226,   227,   228,   229,   232,   233,
   236,   237,   240,   241,   244,   245,   248,   249
};
#endif


#if YYDEBUG != 0 || defined (YYERROR_VERBOSE)

static const char * const yytname[] = {   "$","error","$undefined.","T_NUMBER",
"T_SYMBOL","T_LITERAL","T_DIGITS","T_QUOTE","T_QUOTE2","T_QUOTE3","T_EQUAL",
"T_NOTEQ","T_LESSEQ","T_GREATEREQ","T_QUIT","T_WARRANTY","T_PRINT","T_IPRINT",
"T_TIME","T_XYZZY","T_INVENTORY","T_LOOK","T_SCORE","'='","'<'","'>'","'+'",
"'-'","'*'","'/'","'%'","NEG","'^'","'!'","';'","':'","'('","')'","'?'","'\\''",
"'{'","'}'","'['","']'","','","input","line","@1","exp","exprseq","list_or_empty",
"list","matrix","row", NULL
};
#endif

static const short yyr1[] = {     0,
    45,    45,    46,    46,    46,    46,    46,    46,    46,    46,
    46,    46,    46,    46,    46,    46,    47,    46,    46,    46,
    48,    48,    48,    48,    48,    48,    48,    48,    48,    48,
    48,    48,    48,    48,    48,    48,    48,    48,    48,    48,
    48,    48,    48,    48,    48,    48,    48,    48,    49,    49,
    50,    50,    51,    51,    52,    52,    53,    53
};

static const short yyr2[] = {     0,
     0,     2,     1,     2,     2,     5,     5,     2,     2,     2,
     1,     1,     1,     1,     1,     1,     0,     5,     2,     2,
     1,     1,     3,     1,     1,     1,     1,     1,     4,     3,
     3,     3,     3,     3,     3,     3,     3,     3,     3,     3,
     3,     2,     2,     3,     2,     3,     3,     3,     1,     3,
     0,     1,     1,     3,     3,     5,     1,     3
};

static const short yydefact[] = {     1,
     0,     0,    21,    22,    24,    25,    26,    27,    28,    11,
    12,     0,     0,    17,    13,    14,    15,    16,     0,     0,
     3,     0,     0,     0,    51,     0,     2,     0,    19,    20,
     0,     0,     0,     0,     0,     0,    43,    42,     0,     8,
     9,    10,     0,    53,     0,    52,     0,     0,     0,     0,
     0,     0,     0,     0,     0,     0,     0,     0,     0,    45,
     4,     5,    31,    49,     0,    30,     0,     0,     0,    46,
    23,    47,     0,    57,     0,    48,     0,    32,    33,    35,
    37,    34,    36,    38,    39,    40,    41,    44,    29,     0,
     0,     0,     0,    54,    55,     0,     0,    50,     6,     7,
    18,    58,     0,    56,     0,     0
};

static const short yydefgoto[] = {     1,
    27,    36,    74,    65,    45,    46,    48,    75
};

static const short yypact[] = {-32768,
    93,   -23,-32768,   -13,-32768,   -21,-32768,-32768,-32768,-32768,
-32768,   -22,   -14,-32768,-32768,-32768,-32768,-32768,     0,     0,
-32768,     0,    -3,    21,     0,   -10,-32768,   196,-32768,-32768,
     0,     0,    38,     0,     0,     7,   -16,   -16,    49,-32768,
-32768,-32768,     8,   222,     5,    23,     0,   -15,     0,     0,
     0,     0,     0,     0,     0,     0,     0,     0,     0,-32768,
-32768,-32768,   222,   222,   -24,-32768,   112,   140,     0,-32768,
-32768,-32768,     0,   222,    -6,-32768,    29,   130,   130,    37,
    37,    37,    37,    51,    51,   -16,   -16,   -16,-32768,     0,
    53,    54,   168,   222,-32768,     0,     0,   222,-32768,-32768,
-32768,   222,     1,-32768,    85,-32768
};

static const short yypgoto[] = {-32768,
-32768,-32768,    -1,-32768,-32768,-32768,-32768,    -7
};


#define	YYLAST		255


static const short yytable[] = {    28,
    40,    33,     3,     4,     5,     6,     7,     8,     9,    31,
    29,    30,    89,    34,    41,    59,    60,    37,    38,    90,
    39,    35,    32,    44,    43,    19,    20,    76,    77,    63,
    64,    47,    67,    68,    42,    22,    95,    96,    24,    25,
    66,    26,    69,   104,    96,    72,    71,    78,    79,    80,
    81,    82,    83,    84,    85,    86,    87,    88,    49,    50,
    51,    52,    55,    56,    57,    58,    73,    93,    59,    60,
    97,    94,    53,    54,    55,    56,    57,    58,    57,    58,
    59,    60,    59,    60,   106,    70,    99,   100,    98,   103,
     0,     0,   105,     2,   102,     3,     4,     5,     6,     7,
     8,     9,     0,     0,     0,     0,    10,    11,    12,    13,
    14,    15,    16,    17,    18,     0,     0,     0,    19,    20,
     0,    49,    50,    51,    52,     0,    21,     0,    22,     0,
    23,    24,    25,     0,    26,    53,    54,    55,    56,    57,
    58,    51,    52,    59,    60,     0,     0,     0,    91,    49,
    50,    51,    52,    53,    54,    55,    56,    57,    58,     0,
     0,    59,    60,    53,    54,    55,    56,    57,    58,     0,
     0,    59,    60,     0,     0,     0,    92,    49,    50,    51,
    52,     0,     0,     0,     0,     0,     0,     0,     0,     0,
     0,    53,    54,    55,    56,    57,    58,     0,     0,    59,
    60,     0,     0,     0,   101,    49,    50,    51,    52,     0,
     0,     0,     0,     0,     0,     0,     0,     0,     0,    53,
    54,    55,    56,    57,    58,     0,     0,    59,    60,    61,
    62,    49,    50,    51,    52,     0,     0,     0,     0,     0,
     0,     0,     0,     0,     0,    53,    54,    55,    56,    57,
    58,     0,     0,    59,    60
};

static const short yycheck[] = {     1,
     4,    23,     3,     4,     5,     6,     7,     8,     9,    23,
    34,    35,    37,    36,    18,    32,    33,    19,    20,    44,
    22,    36,    36,    25,     4,    26,    27,    43,    44,    31,
    32,    42,    34,    35,    38,    36,    43,    44,    39,    40,
     3,    42,    36,    43,    44,    41,    39,    49,    50,    51,
    52,    53,    54,    55,    56,    57,    58,    59,    10,    11,
    12,    13,    26,    27,    28,    29,    44,    69,    32,    33,
    42,    73,    24,    25,    26,    27,    28,    29,    28,    29,
    32,    33,    32,    33,     0,    37,    34,    34,    90,    97,
    -1,    -1,     0,     1,    96,     3,     4,     5,     6,     7,
     8,     9,    -1,    -1,    -1,    -1,    14,    15,    16,    17,
    18,    19,    20,    21,    22,    -1,    -1,    -1,    26,    27,
    -1,    10,    11,    12,    13,    -1,    34,    -1,    36,    -1,
    38,    39,    40,    -1,    42,    24,    25,    26,    27,    28,
    29,    12,    13,    32,    33,    -1,    -1,    -1,    37,    10,
    11,    12,    13,    24,    25,    26,    27,    28,    29,    -1,
    -1,    32,    33,    24,    25,    26,    27,    28,    29,    -1,
    -1,    32,    33,    -1,    -1,    -1,    37,    10,    11,    12,
    13,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
    -1,    24,    25,    26,    27,    28,    29,    -1,    -1,    32,
    33,    -1,    -1,    -1,    37,    10,    11,    12,    13,    -1,
    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    24,
    25,    26,    27,    28,    29,    -1,    -1,    32,    33,    34,
    35,    10,    11,    12,    13,    -1,    -1,    -1,    -1,    -1,
    -1,    -1,    -1,    -1,    -1,    24,    25,    26,    27,    28,
    29,    -1,    -1,    32,    33
};
/* -*-C-*-  Note some compilers choke on comments on `#line' lines.  */
#line 3 "/usr/lib/bison.simple"
/* This file comes from bison-1.28.  */

/* Skeleton output parser for bison,
   Copyright (C) 1984, 1989, 1990 Free Software Foundation, Inc.

   This program is free software; you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation; either version 2, or (at your option)
   any later version.

   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.

   You should have received a copy of the GNU General Public License
   along with this program; if not, write to the Free Software
   Foundation, Inc., 59 Temple Place - Suite 330,
   Boston, MA 02111-1307, USA.  */

/* As a special exception, when this file is copied by Bison into a
   Bison output file, you may use that output file without restriction.
   This special exception was added by the Free Software Foundation
   in version 1.24 of Bison.  */

/* This is the parser code that is written into each bison parser
  when the %semantic_parser declaration is not specified in the grammar.
  It was written by Richard Stallman by simplifying the hairy parser
  used when %semantic_parser is specified.  */

#ifndef YYSTACK_USE_ALLOCA
#ifdef alloca
#define YYSTACK_USE_ALLOCA
#else /* alloca not defined */
#ifdef __GNUC__
#define YYSTACK_USE_ALLOCA
#define alloca __builtin_alloca
#else /* not GNU C.  */
#if (!defined (__STDC__) && defined (sparc)) || defined (__sparc__) || defined (__sparc) || defined (__sgi) || (defined (__sun) && defined (__i386))
#define YYSTACK_USE_ALLOCA
#include <alloca.h>
#else /* not sparc */
/* We think this test detects Watcom and Microsoft C.  */
/* This used to test MSDOS, but that is a bad idea
   since that symbol is in the user namespace.  */
#if (defined (_MSDOS) || defined (_MSDOS_)) && !defined (__TURBOC__)
#if 0 /* No need for malloc.h, which pollutes the namespace;
	 instead, just don't use alloca.  */
#include <malloc.h>
#endif
#else /* not MSDOS, or __TURBOC__ */
#if defined(_AIX)
/* I don't know what this was needed for, but it pollutes the namespace.
   So I turned it off.   rms, 2 May 1997.  */
/* #include <malloc.h>  */
 #pragma alloca
#define YYSTACK_USE_ALLOCA
#else /* not MSDOS, or __TURBOC__, or _AIX */
#if 0
#ifdef __hpux /* haible@ilog.fr says this works for HPUX 9.05 and up,
		 and on HPUX 10.  Eventually we can turn this on.  */
#define YYSTACK_USE_ALLOCA
#define alloca __builtin_alloca
#endif /* __hpux */
#endif
#endif /* not _AIX */
#endif /* not MSDOS, or __TURBOC__ */
#endif /* not sparc */
#endif /* not GNU C */
#endif /* alloca not defined */
#endif /* YYSTACK_USE_ALLOCA not defined */

#ifdef YYSTACK_USE_ALLOCA
#define YYSTACK_ALLOC alloca
#else
#define YYSTACK_ALLOC malloc
#endif

/* Note: there must be only one dollar sign in this file.
   It is replaced by the list of actions, each action
   as one case of the switch.  */

#define yyerrok		(yyerrstatus = 0)
#define yyclearin	(yychar = YYEMPTY)
#define YYEMPTY		-2
#define YYEOF		0
#define YYACCEPT	goto yyacceptlab
#define YYABORT 	goto yyabortlab
#define YYERROR		goto yyerrlab1
/* Like YYERROR except do call yyerror.
   This remains here temporarily to ease the
   transition to the new meaning of YYERROR, for GCC.
   Once GCC version 2 has supplanted version 1, this can go.  */
#define YYFAIL		goto yyerrlab
#define YYRECOVERING()  (!!yyerrstatus)
#define YYBACKUP(token, value) \
do								\
  if (yychar == YYEMPTY && yylen == 1)				\
    { yychar = (token), yylval = (value);			\
      yychar1 = YYTRANSLATE (yychar);				\
      YYPOPSTACK;						\
      goto yybackup;						\
    }								\
  else								\
    { yyerror ("syntax error: cannot back up"); YYERROR; }	\
while (0)

#define YYTERROR	1
#define YYERRCODE	256

#ifndef YYPURE
#define YYLEX		yylex()
#endif

#ifdef YYPURE
#ifdef YYLSP_NEEDED
#ifdef YYLEX_PARAM
#define YYLEX		yylex(&yylval, &yylloc, YYLEX_PARAM)
#else
#define YYLEX		yylex(&yylval, &yylloc)
#endif
#else /* not YYLSP_NEEDED */
#ifdef YYLEX_PARAM
#define YYLEX		yylex(&yylval, YYLEX_PARAM)
#else
#define YYLEX		yylex(&yylval)
#endif
#endif /* not YYLSP_NEEDED */
#endif

/* If nonreentrant, generate the variables here */

#ifndef YYPURE

int	yychar;			/*  the lookahead symbol		*/
YYSTYPE	yylval;			/*  the semantic value of the		*/
				/*  lookahead symbol			*/

#ifdef YYLSP_NEEDED
YYLTYPE yylloc;			/*  location data for the lookahead	*/
				/*  symbol				*/
#endif

int yynerrs;			/*  number of parse errors so far       */
#endif  /* not YYPURE */

#if YYDEBUG != 0
int yydebug;			/*  nonzero means print parse trace	*/
/* Since this is uninitialized, it does not stop multiple parsers
   from coexisting.  */
#endif

/*  YYINITDEPTH indicates the initial size of the parser's stacks	*/

#ifndef	YYINITDEPTH
#define YYINITDEPTH 200
#endif

/*  YYMAXDEPTH is the maximum size the stacks can grow to
    (effective only if the built-in stack extension method is used).  */

#if YYMAXDEPTH == 0
#undef YYMAXDEPTH
#endif

#ifndef YYMAXDEPTH
#define YYMAXDEPTH 10000
#endif

/* Define __yy_memcpy.  Note that the size argument
   should be passed with type unsigned int, because that is what the non-GCC
   definitions require.  With GCC, __builtin_memcpy takes an arg
   of type size_t, but it can handle unsigned int.  */

#if __GNUC__ > 1		/* GNU C and GNU C++ define this.  */
#define __yy_memcpy(TO,FROM,COUNT)	__builtin_memcpy(TO,FROM,COUNT)
#else				/* not GNU C or C++ */
#ifndef __cplusplus

/* This is the most reliable way to avoid incompatibilities
   in available built-in functions on various systems.  */
static void
__yy_memcpy (to, from, count)
     char *to;
     char *from;
     unsigned int count;
{
  register char *f = from;
  register char *t = to;
  register int i = count;

  while (i-- > 0)
    *t++ = *f++;
}

#else /* __cplusplus */

/* This is the most reliable way to avoid incompatibilities
   in available built-in functions on various systems.  */
static void
__yy_memcpy (char *to, char *from, unsigned int count)
{
  register char *t = to;
  register char *f = from;
  register int i = count;

  while (i-- > 0)
    *t++ = *f++;
}

#endif
#endif

#line 217 "/usr/lib/bison.simple"

/* The user can define YYPARSE_PARAM as the name of an argument to be passed
   into yyparse.  The argument should have type void *.
   It should actually point to an object.
   Grammar actions can access the variable by casting it
   to the proper pointer type.  */

#ifdef YYPARSE_PARAM
#ifdef __cplusplus
#define YYPARSE_PARAM_ARG void *YYPARSE_PARAM
#define YYPARSE_PARAM_DECL
#else /* not __cplusplus */
#define YYPARSE_PARAM_ARG YYPARSE_PARAM
#define YYPARSE_PARAM_DECL void *YYPARSE_PARAM;
#endif /* not __cplusplus */
#else /* not YYPARSE_PARAM */
#define YYPARSE_PARAM_ARG
#define YYPARSE_PARAM_DECL
#endif /* not YYPARSE_PARAM */

/* Prevent warning if -Wstrict-prototypes.  */
#ifdef __GNUC__
#ifdef YYPARSE_PARAM
int yyparse (void *);
#else
int yyparse (void);
#endif
#endif

int
yyparse(YYPARSE_PARAM_ARG)
     YYPARSE_PARAM_DECL
{
  register int yystate;
  register int yyn;
  register short *yyssp;
  register YYSTYPE *yyvsp;
  int yyerrstatus;	/*  number of tokens to shift before error messages enabled */
  int yychar1 = 0;		/*  lookahead token as an internal (translated) token number */

  short	yyssa[YYINITDEPTH];	/*  the state stack			*/
  YYSTYPE yyvsa[YYINITDEPTH];	/*  the semantic value stack		*/

  short *yyss = yyssa;		/*  refer to the stacks thru separate pointers */
  YYSTYPE *yyvs = yyvsa;	/*  to allow yyoverflow to reallocate them elsewhere */

#ifdef YYLSP_NEEDED
  YYLTYPE yylsa[YYINITDEPTH];	/*  the location stack			*/
  YYLTYPE *yyls = yylsa;
  YYLTYPE *yylsp;

#define YYPOPSTACK   (yyvsp--, yyssp--, yylsp--)
#else
#define YYPOPSTACK   (yyvsp--, yyssp--)
#endif

  int yystacksize = YYINITDEPTH;
  int yyfree_stacks = 0;

#ifdef YYPURE
  int yychar;
  YYSTYPE yylval;
  int yynerrs;
#ifdef YYLSP_NEEDED
  YYLTYPE yylloc;
#endif
#endif

  YYSTYPE yyval;		/*  the variable used to return		*/
				/*  semantic values from the action	*/
				/*  routines				*/

  int yylen;

#if YYDEBUG != 0
  if (yydebug)
    fprintf(stderr, "Starting parse\n");
#endif

  yystate = 0;
  yyerrstatus = 0;
  yynerrs = 0;
  yychar = YYEMPTY;		/* Cause a token to be read.  */

  /* Initialize stack pointers.
     Waste one element of value and location stack
     so that they stay on the same level as the state stack.
     The wasted elements are never initialized.  */

  yyssp = yyss - 1;
  yyvsp = yyvs;
#ifdef YYLSP_NEEDED
  yylsp = yyls;
#endif

/* Push a new state, which is found in  yystate  .  */
/* In all cases, when you get here, the value and location stacks
   have just been pushed. so pushing a state here evens the stacks.  */
yynewstate:

  *++yyssp = yystate;

  if (yyssp >= yyss + yystacksize - 1)
    {
      /* Give user a chance to reallocate the stack */
      /* Use copies of these so that the &'s don't force the real ones into memory. */
      YYSTYPE *yyvs1 = yyvs;
      short *yyss1 = yyss;
#ifdef YYLSP_NEEDED
      YYLTYPE *yyls1 = yyls;
#endif

      /* Get the current used size of the three stacks, in elements.  */
      int size = yyssp - yyss + 1;

#ifdef yyoverflow
      /* Each stack pointer address is followed by the size of
	 the data in use in that stack, in bytes.  */
#ifdef YYLSP_NEEDED
      /* This used to be a conditional around just the two extra args,
	 but that might be undefined if yyoverflow is a macro.  */
      yyoverflow("parser stack overflow",
		 &yyss1, size * sizeof (*yyssp),
		 &yyvs1, size * sizeof (*yyvsp),
		 &yyls1, size * sizeof (*yylsp),
		 &yystacksize);
#else
      yyoverflow("parser stack overflow",
		 &yyss1, size * sizeof (*yyssp),
		 &yyvs1, size * sizeof (*yyvsp),
		 &yystacksize);
#endif

      yyss = yyss1; yyvs = yyvs1;
#ifdef YYLSP_NEEDED
      yyls = yyls1;
#endif
#else /* no yyoverflow */
      /* Extend the stack our own way.  */
      if (yystacksize >= YYMAXDEPTH)
	{
	  yyerror("parser stack overflow");
	  if (yyfree_stacks)
	    {
	      free (yyss);
	      free (yyvs);
#ifdef YYLSP_NEEDED
	      free (yyls);
#endif
	    }
	  return 2;
	}
      yystacksize *= 2;
      if (yystacksize > YYMAXDEPTH)
	yystacksize = YYMAXDEPTH;
#ifndef YYSTACK_USE_ALLOCA
      yyfree_stacks = 1;
#endif
      yyss = (short *) YYSTACK_ALLOC (yystacksize * sizeof (*yyssp));
      __yy_memcpy ((char *)yyss, (char *)yyss1,
		   size * (unsigned int) sizeof (*yyssp));
      yyvs = (YYSTYPE *) YYSTACK_ALLOC (yystacksize * sizeof (*yyvsp));
      __yy_memcpy ((char *)yyvs, (char *)yyvs1,
		   size * (unsigned int) sizeof (*yyvsp));
#ifdef YYLSP_NEEDED
      yyls = (YYLTYPE *) YYSTACK_ALLOC (yystacksize * sizeof (*yylsp));
      __yy_memcpy ((char *)yyls, (char *)yyls1,
		   size * (unsigned int) sizeof (*yylsp));
#endif
#endif /* no yyoverflow */

      yyssp = yyss + size - 1;
      yyvsp = yyvs + size - 1;
#ifdef YYLSP_NEEDED
      yylsp = yyls + size - 1;
#endif

#if YYDEBUG != 0
      if (yydebug)
	fprintf(stderr, "Stack size increased to %d\n", yystacksize);
#endif

      if (yyssp >= yyss + yystacksize - 1)
	YYABORT;
    }

#if YYDEBUG != 0
  if (yydebug)
    fprintf(stderr, "Entering state %d\n", yystate);
#endif

  goto yybackup;
 yybackup:

/* Do appropriate processing given the current state.  */
/* Read a lookahead token if we need one and don't already have one.  */
/* yyresume: */

  /* First try to decide what to do without reference to lookahead token.  */

  yyn = yypact[yystate];
  if (yyn == YYFLAG)
    goto yydefault;

  /* Not known => get a lookahead token if don't already have one.  */

  /* yychar is either YYEMPTY or YYEOF
     or a valid token in external form.  */

  if (yychar == YYEMPTY)
    {
#if YYDEBUG != 0
      if (yydebug)
	fprintf(stderr, "Reading a token: ");
#endif
      yychar = YYLEX;
    }

  /* Convert token to internal form (in yychar1) for indexing tables with */

  if (yychar <= 0)		/* This means end of input. */
    {
      yychar1 = 0;
      yychar = YYEOF;		/* Don't call YYLEX any more */

#if YYDEBUG != 0
      if (yydebug)
	fprintf(stderr, "Now at end of input.\n");
#endif
    }
  else
    {
      yychar1 = YYTRANSLATE(yychar);

#if YYDEBUG != 0
      if (yydebug)
	{
	  fprintf (stderr, "Next token is %d (%s", yychar, yytname[yychar1]);
	  /* Give the individual parser a way to print the precise meaning
	     of a token, for further debugging info.  */
#ifdef YYPRINT
	  YYPRINT (stderr, yychar, yylval);
#endif
	  fprintf (stderr, ")\n");
	}
#endif
    }

  yyn += yychar1;
  if (yyn < 0 || yyn > YYLAST || yycheck[yyn] != yychar1)
    goto yydefault;

  yyn = yytable[yyn];

  /* yyn is what to do for this token type in this state.
     Negative => reduce, -yyn is rule number.
     Positive => shift, yyn is new state.
       New state is final state => don't bother to shift,
       just return success.
     0, or most negative number => error.  */

  if (yyn < 0)
    {
      if (yyn == YYFLAG)
	goto yyerrlab;
      yyn = -yyn;
      goto yyreduce;
    }
  else if (yyn == 0)
    goto yyerrlab;

  if (yyn == YYFINAL)
    YYACCEPT;

  /* Shift the lookahead token.  */

#if YYDEBUG != 0
  if (yydebug)
    fprintf(stderr, "Shifting token %d (%s), ", yychar, yytname[yychar1]);
#endif

  /* Discard the token being shifted unless it is eof.  */
  if (yychar != YYEOF)
    yychar = YYEMPTY;

  *++yyvsp = yylval;
#ifdef YYLSP_NEEDED
  *++yylsp = yylloc;
#endif

  /* count tokens shifted since error; after three, turn off error status.  */
  if (yyerrstatus) yyerrstatus--;

  yystate = yyn;
  goto yynewstate;

/* Do the default action for the current state.  */
yydefault:

  yyn = yydefact[yystate];
  if (yyn == 0)
    goto yyerrlab;

/* Do a reduction.  yyn is the number of a rule to reduce with.  */
yyreduce:
  yylen = yyr2[yyn];
  if (yylen > 0)
    yyval = yyvsp[1-yylen]; /* implement default value of the action */

#if YYDEBUG != 0
  if (yydebug)
    {
      int i;

      fprintf (stderr, "Reducing via rule %d (line %d), ",
	       yyn, yyrline[yyn]);

      /* Print the symbols being reduced, and their result.  */
      for (i = yyprhs[yyn]; yyrhs[i] > 0; i++)
	fprintf (stderr, "%s ", yytname[yyrhs[i]]);
      fprintf (stderr, " -> %s\n", yytname[yyr1[yyn]]);
    }
#endif


  switch (yyn) {

case 4:
#line 120 "ginsh_parser.yy"
{
		try {
			cout << yyvsp[-1] << endl;
			push(yyvsp[-1]);
		} catch (exception &e) {
			cerr << e.what() << endl;
			YYERROR;
		}
	;
    break;}
case 5:
#line 129 "ginsh_parser.yy"
{
		try {
			push(yyvsp[-1]);
		} catch (exception &e) {
			std::cerr << e.what() << endl;
			YYERROR;
		}
	;
    break;}
case 6:
#line 137 "ginsh_parser.yy"
{
		try {
			yyvsp[-2].print(print_tree(std::cout));
		} catch (exception &e) {
			std::cerr << e.what() << endl;
			YYERROR;
		}
	;
    break;}
case 7:
#line 145 "ginsh_parser.yy"
{
		try {
			ex e = yyvsp[-2];
			if (!e.info(info_flags::integer))
				throw (std::invalid_argument("argument to iprint() must be an integer"));
			long i = ex_to<numeric>(e).to_long();
			cout << i << endl;
			cout << "#o" << oct << i << endl;
			cout << "#x" << hex << i << dec << endl;
		} catch (exception &e) {
			cerr << e.what() << endl;
			YYERROR;
		}
	;
    break;}
case 8:
#line 159 "ginsh_parser.yy"
{print_help(ex_to<symbol>(yyvsp[0]).get_name());;
    break;}
case 9:
#line 160 "ginsh_parser.yy"
{print_help("time");;
    break;}
case 10:
#line 161 "ginsh_parser.yy"
{print_help_topics();;
    break;}
case 11:
#line 162 "ginsh_parser.yy"
{YYACCEPT;;
    break;}
case 12:
#line 163 "ginsh_parser.yy"
{
		cout << "This program is free software; you can redistribute it and/or modify it under\n";
		cout << "the terms of the GNU General Public License as published by the Free Software\n";
		cout << "Foundation; either version 2 of the License, or (at your option) any later\n";
		cout << "version.\n";
		cout << "This program is distributed in the hope that it will be useful, but WITHOUT\n";
		cout << "ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS\n";
		cout << "FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more\n";
		cout << "details.\n";
		cout << "You should have received a copy of the GNU General Public License along with\n";
		cout << "this program. If not, write to the Free Software Foundation, 675 Mass Ave,\n";
		cout << "Cambridge, MA 02139, USA.\n";
	;
    break;}
case 13:
#line 176 "ginsh_parser.yy"
{cout << "Nothing happens.\n";;
    break;}
case 14:
#line 177 "ginsh_parser.yy"
{cout << "You're not carrying anything.\n";;
    break;}
case 15:
#line 178 "ginsh_parser.yy"
{cout << "You're in a twisty little maze of passages, all alike.\n";;
    break;}
case 16:
#line 179 "ginsh_parser.yy"
{
		cout << "If you were to quit now, you would score ";
		cout << (syms.size() > 350 ? 350 : syms.size());
		cout << " out of a possible 350.\n";
	;
    break;}
case 17:
#line 184 "ginsh_parser.yy"
{getrusage(RUSAGE_SELF, &start_time);;
    break;}
case 18:
#line 184 "ginsh_parser.yy"
{
		getrusage(RUSAGE_SELF, &end_time);
		cout << (end_time.ru_utime.tv_sec - start_time.ru_utime.tv_sec) +
			(end_time.ru_stime.tv_sec - start_time.ru_stime.tv_sec) +
			 double(end_time.ru_utime.tv_usec - start_time.ru_utime.tv_usec) / 1e6 +
			 double(end_time.ru_stime.tv_usec - start_time.ru_stime.tv_usec) / 1e6 << 's' << endl;
	;
    break;}
case 19:
#line 191 "ginsh_parser.yy"
{yyclearin; yyerrok;;
    break;}
case 20:
#line 192 "ginsh_parser.yy"
{yyclearin; yyerrok;;
    break;}
case 21:
#line 195 "ginsh_parser.yy"
{yyval = yyvsp[0];;
    break;}
case 22:
#line 196 "ginsh_parser.yy"
{yyval = yyvsp[0].eval();;
    break;}
case 23:
#line 197 "ginsh_parser.yy"
{yyval = yyvsp[-1];;
    break;}
case 24:
#line 198 "ginsh_parser.yy"
{yyval = yyvsp[0];;
    break;}
case 25:
#line 199 "ginsh_parser.yy"
{yyval = yyvsp[0];;
    break;}
case 26:
#line 200 "ginsh_parser.yy"
{yyval = exstack[0];;
    break;}
case 27:
#line 201 "ginsh_parser.yy"
{yyval = exstack[1];;
    break;}
case 28:
#line 202 "ginsh_parser.yy"
{yyval = exstack[2];;
    break;}
case 29:
#line 203 "ginsh_parser.yy"
{
		fcn_tab::const_iterator i = find_function(yyvsp[-3], yyvsp[-1].nops());
		if (i->second.is_ginac) {
			yyval = ((fcnp2)(i->second.p))(ex_to<exprseq>(yyvsp[-1]), i->second.serial);
		} else {
			yyval = (i->second.p)(ex_to<exprseq>(yyvsp[-1]));
		}
	;
    break;}
case 30:
#line 211 "ginsh_parser.yy"
{yyval = yyvsp[0]; Digits = ex_to<numeric>(yyvsp[0]).to_int();;
    break;}
case 31:
#line 212 "ginsh_parser.yy"
{yyval = yyvsp[0]; const_cast<symbol&>(ex_to<symbol>(yyvsp[-2])).assign(yyvsp[0]);;
    break;}
case 32:
#line 213 "ginsh_parser.yy"
{yyval = yyvsp[-2] == yyvsp[0];;
    break;}
case 33:
#line 214 "ginsh_parser.yy"
{yyval = yyvsp[-2] != yyvsp[0];;
    break;}
case 34:
#line 215 "ginsh_parser.yy"
{yyval = yyvsp[-2] < yyvsp[0];;
    break;}
case 35:
#line 216 "ginsh_parser.yy"
{yyval = yyvsp[-2] <= yyvsp[0];;
    break;}
case 36:
#line 217 "ginsh_parser.yy"
{yyval = yyvsp[-2] > yyvsp[0];;
    break;}
case 37:
#line 218 "ginsh_parser.yy"
{yyval = yyvsp[-2] >= yyvsp[0];;
    break;}
case 38:
#line 219 "ginsh_parser.yy"
{yyval = yyvsp[-2] + yyvsp[0];;
    break;}
case 39:
#line 220 "ginsh_parser.yy"
{yyval = yyvsp[-2] - yyvsp[0];;
    break;}
case 40:
#line 221 "ginsh_parser.yy"
{yyval = yyvsp[-2] * yyvsp[0];;
    break;}
case 41:
#line 222 "ginsh_parser.yy"
{yyval = yyvsp[-2] / yyvsp[0];;
    break;}
case 42:
#line 223 "ginsh_parser.yy"
{yyval = -yyvsp[0];;
    break;}
case 43:
#line 224 "ginsh_parser.yy"
{yyval = yyvsp[0];;
    break;}
case 44:
#line 225 "ginsh_parser.yy"
{yyval = power(yyvsp[-2], yyvsp[0]);;
    break;}
case 45:
#line 226 "ginsh_parser.yy"
{yyval = factorial(yyvsp[-1]);;
    break;}
case 46:
#line 227 "ginsh_parser.yy"
{yyval = yyvsp[-1];;
    break;}
case 47:
#line 228 "ginsh_parser.yy"
{yyval = yyvsp[-1];;
    break;}
case 48:
#line 229 "ginsh_parser.yy"
{yyval = lst_to_matrix(ex_to<lst>(yyvsp[-1]));;
    break;}
case 49:
#line 232 "ginsh_parser.yy"
{yyval = exprseq(yyvsp[0]);;
    break;}
case 50:
#line 233 "ginsh_parser.yy"
{exprseq es(ex_to<exprseq>(yyvsp[-2])); yyval = es.append(yyvsp[0]);;
    break;}
case 51:
#line 236 "ginsh_parser.yy"
{yyval = *new lst;;
    break;}
case 52:
#line 237 "ginsh_parser.yy"
{yyval = yyvsp[0];;
    break;}
case 53:
#line 240 "ginsh_parser.yy"
{yyval = lst(yyvsp[0]);;
    break;}
case 54:
#line 241 "ginsh_parser.yy"
{lst l(ex_to<lst>(yyvsp[-2])); yyval = l.append(yyvsp[0]);;
    break;}
case 55:
#line 244 "ginsh_parser.yy"
{yyval = lst(yyvsp[-1]);;
    break;}
case 56:
#line 245 "ginsh_parser.yy"
{lst l(ex_to<lst>(yyvsp[-4])); yyval = l.append(yyvsp[-1]);;
    break;}
case 57:
#line 248 "ginsh_parser.yy"
{yyval = lst(yyvsp[0]);;
    break;}
case 58:
#line 249 "ginsh_parser.yy"
{lst l(ex_to<lst>(yyvsp[-2])); yyval = l.append(yyvsp[0]);;
    break;}
}
   /* the action file gets copied in in place of this dollarsign */
#line 543 "/usr/lib/bison.simple"

  yyvsp -= yylen;
  yyssp -= yylen;
#ifdef YYLSP_NEEDED
  yylsp -= yylen;
#endif

#if YYDEBUG != 0
  if (yydebug)
    {
      short *ssp1 = yyss - 1;
      fprintf (stderr, "state stack now");
      while (ssp1 != yyssp)
	fprintf (stderr, " %d", *++ssp1);
      fprintf (stderr, "\n");
    }
#endif

  *++yyvsp = yyval;

#ifdef YYLSP_NEEDED
  yylsp++;
  if (yylen == 0)
    {
      yylsp->first_line = yylloc.first_line;
      yylsp->first_column = yylloc.first_column;
      yylsp->last_line = (yylsp-1)->last_line;
      yylsp->last_column = (yylsp-1)->last_column;
      yylsp->text = 0;
    }
  else
    {
      yylsp->last_line = (yylsp+yylen-1)->last_line;
      yylsp->last_column = (yylsp+yylen-1)->last_column;
    }
#endif

  /* Now "shift" the result of the reduction.
     Determine what state that goes to,
     based on the state we popped back to
     and the rule number reduced by.  */

  yyn = yyr1[yyn];

  yystate = yypgoto[yyn - YYNTBASE] + *yyssp;
  if (yystate >= 0 && yystate <= YYLAST && yycheck[yystate] == *yyssp)
    yystate = yytable[yystate];
  else
    yystate = yydefgoto[yyn - YYNTBASE];

  goto yynewstate;

yyerrlab:   /* here on detecting error */

  if (! yyerrstatus)
    /* If not already recovering from an error, report this error.  */
    {
      ++yynerrs;

#ifdef YYERROR_VERBOSE
      yyn = yypact[yystate];

      if (yyn > YYFLAG && yyn < YYLAST)
	{
	  int size = 0;
	  char *msg;
	  int x, count;

	  count = 0;
	  /* Start X at -yyn if nec to avoid negative indexes in yycheck.  */
	  for (x = (yyn < 0 ? -yyn : 0);
	       x < (sizeof(yytname) / sizeof(char *)); x++)
	    if (yycheck[x + yyn] == x)
	      size += strlen(yytname[x]) + 15, count++;
	  msg = (char *) malloc(size + 15);
	  if (msg != 0)
	    {
	      strcpy(msg, "parse error");

	      if (count < 5)
		{
		  count = 0;
		  for (x = (yyn < 0 ? -yyn : 0);
		       x < (sizeof(yytname) / sizeof(char *)); x++)
		    if (yycheck[x + yyn] == x)
		      {
			strcat(msg, count == 0 ? ", expecting `" : " or `");
			strcat(msg, yytname[x]);
			strcat(msg, "'");
			count++;
		      }
		}
	      yyerror(msg);
	      free(msg);
	    }
	  else
	    yyerror ("parse error; also virtual memory exceeded");
	}
      else
#endif /* YYERROR_VERBOSE */
	yyerror("parse error");
    }

  goto yyerrlab1;
yyerrlab1:   /* here on error raised explicitly by an action */

  if (yyerrstatus == 3)
    {
      /* if just tried and failed to reuse lookahead token after an error, discard it.  */

      /* return failure if at end of input */
      if (yychar == YYEOF)
	YYABORT;

#if YYDEBUG != 0
      if (yydebug)
	fprintf(stderr, "Discarding token %d (%s).\n", yychar, yytname[yychar1]);
#endif

      yychar = YYEMPTY;
    }

  /* Else will try to reuse lookahead token
     after shifting the error token.  */

  yyerrstatus = 3;		/* Each real token shifted decrements this */

  goto yyerrhandle;

yyerrdefault:  /* current state does not do anything special for the error token. */

#if 0
  /* This is wrong; only states that explicitly want error tokens
     should shift them.  */
  yyn = yydefact[yystate];  /* If its default is to accept any token, ok.  Otherwise pop it.*/
  if (yyn) goto yydefault;
#endif

yyerrpop:   /* pop the current state because it cannot handle the error token */

  if (yyssp == yyss) YYABORT;
  yyvsp--;
  yystate = *--yyssp;
#ifdef YYLSP_NEEDED
  yylsp--;
#endif

#if YYDEBUG != 0
  if (yydebug)
    {
      short *ssp1 = yyss - 1;
      fprintf (stderr, "Error: state stack now");
      while (ssp1 != yyssp)
	fprintf (stderr, " %d", *++ssp1);
      fprintf (stderr, "\n");
    }
#endif

yyerrhandle:

  yyn = yypact[yystate];
  if (yyn == YYFLAG)
    goto yyerrdefault;

  yyn += YYTERROR;
  if (yyn < 0 || yyn > YYLAST || yycheck[yyn] != YYTERROR)
    goto yyerrdefault;

  yyn = yytable[yyn];
  if (yyn < 0)
    {
      if (yyn == YYFLAG)
	goto yyerrpop;
      yyn = -yyn;
      goto yyreduce;
    }
  else if (yyn == 0)
    goto yyerrpop;

  if (yyn == YYFINAL)
    YYACCEPT;

#if YYDEBUG != 0
  if (yydebug)
    fprintf(stderr, "Shifting error token, ");
#endif

  *++yyvsp = yylval;
#ifdef YYLSP_NEEDED
  *++yylsp = yylloc;
#endif

  yystate = yyn;
  goto yynewstate;

 yyacceptlab:
  /* YYACCEPT comes here.  */
  if (yyfree_stacks)
    {
      free (yyss);
      free (yyvs);
#ifdef YYLSP_NEEDED
      free (yyls);
#endif
    }
  return 0;

 yyabortlab:
  /* YYABORT comes here.  */
  if (yyfree_stacks)
    {
      free (yyss);
      free (yyvs);
#ifdef YYLSP_NEEDED
      free (yyls);
#endif
    }
  return 1;
}
#line 257 "ginsh_parser.yy"

// Error print routine
int yyerror(char *s)
{
	cerr << s << " at " << yytext << endl;
	return 0;
}

// Push expression "e" onto the expression stack (for ", "" and """)
static void push(const ex &e)
{
	exstack[2] = exstack[1];
	exstack[1] = exstack[0];
	exstack[0] = e;
}


/*
 *  Built-in functions
 */

static ex f_collect(const exprseq &e) {return e[0].collect(e[1]);}
static ex f_collect_distributed(const exprseq &e) {return e[0].collect(e[1], true);}
static ex f_degree(const exprseq &e) {return e[0].degree(e[1]);}
static ex f_denom(const exprseq &e) {return e[0].denom();}
static ex f_eval1(const exprseq &e) {return e[0].eval();}
static ex f_evalf1(const exprseq &e) {return e[0].evalf();}
static ex f_evalm(const exprseq &e) {return e[0].evalm();}
static ex f_expand(const exprseq &e) {return e[0].expand();}
static ex f_gcd(const exprseq &e) {return gcd(e[0], e[1]);}
static ex f_has(const exprseq &e) {return e[0].has(e[1]) ? ex(1) : ex(0);}
static ex f_lcm(const exprseq &e) {return lcm(e[0], e[1]);}
static ex f_lcoeff(const exprseq &e) {return e[0].lcoeff(e[1]);}
static ex f_ldegree(const exprseq &e) {return e[0].ldegree(e[1]);}
static ex f_lsolve(const exprseq &e) {return lsolve(e[0], e[1]);}
static ex f_nops(const exprseq &e) {return e[0].nops();}
static ex f_normal1(const exprseq &e) {return e[0].normal();}
static ex f_numer(const exprseq &e) {return e[0].numer();}
static ex f_numer_denom(const exprseq &e) {return e[0].numer_denom();}
static ex f_pow(const exprseq &e) {return pow(e[0], e[1]);}
static ex f_sqrt(const exprseq &e) {return sqrt(e[0]);}
static ex f_sqrfree1(const exprseq &e) {return sqrfree(e[0]);}
static ex f_subs2(const exprseq &e) {return e[0].subs(e[1]);}
static ex f_tcoeff(const exprseq &e) {return e[0].tcoeff(e[1]);}

#define CHECK_ARG(num, type, fcn) if (!is_a<type>(e[num])) throw(std::invalid_argument("argument " #num " to " #fcn "() must be a " #type))

static ex f_charpoly(const exprseq &e)
{
	CHECK_ARG(0, matrix, charpoly);
	CHECK_ARG(1, symbol, charpoly);
	return ex_to<matrix>(e[0]).charpoly(ex_to<symbol>(e[1]));
}

static ex f_coeff(const exprseq &e)
{
	CHECK_ARG(2, numeric, coeff);
	return e[0].coeff(e[1], ex_to<numeric>(e[2]).to_int());
}

static ex f_content(const exprseq &e)
{
	CHECK_ARG(1, symbol, content);
	return e[0].content(ex_to<symbol>(e[1]));
}

static ex f_decomp_rational(const exprseq &e)
{
	CHECK_ARG(1, symbol, decomp_rational);
	return decomp_rational(e[0], ex_to<symbol>(e[1]));
}

static ex f_determinant(const exprseq &e)
{
	CHECK_ARG(0, matrix, determinant);
	return ex_to<matrix>(e[0]).determinant();
}

static ex f_diag(const exprseq &e)
{
	unsigned dim = e.nops();
	matrix &m = *new matrix(dim, dim);
	for (unsigned i=0; i<dim; i++)
		m.set(i, i, e.op(i));
	return m;
}

static ex f_diff2(const exprseq &e)
{
	CHECK_ARG(1, symbol, diff);
	return e[0].diff(ex_to<symbol>(e[1]));
}

static ex f_diff3(const exprseq &e)
{
	CHECK_ARG(1, symbol, diff);
	CHECK_ARG(2, numeric, diff);
	return e[0].diff(ex_to<symbol>(e[1]), ex_to<numeric>(e[2]).to_int());
}

static ex f_divide(const exprseq &e)
{
	ex q;
	if (divide(e[0], e[1], q))
		return q;
	else
		return fail();
}

static ex f_eval2(const exprseq &e)
{
	CHECK_ARG(1, numeric, eval);
	return e[0].eval(ex_to<numeric>(e[1]).to_int());
}

static ex f_evalf2(const exprseq &e)
{
	CHECK_ARG(1, numeric, evalf);
	return e[0].evalf(ex_to<numeric>(e[1]).to_int());
}

static ex f_find(const exprseq &e)
{
	lst found;
	e[0].find(e[1], found);
	return found;
}

static ex f_inverse(const exprseq &e)
{
	CHECK_ARG(0, matrix, inverse);
	return ex_to<matrix>(e[0]).inverse();
}

static ex f_is(const exprseq &e)
{
	CHECK_ARG(0, relational, is);
	return (bool)ex_to<relational>(e[0]) ? ex(1) : ex(0);
}

class apply_map_function : public map_function {
	ex apply;
public:
	apply_map_function(const ex & a) : apply(a) {}
	virtual ~apply_map_function() {}
	ex operator()(const ex & e) { return apply.subs(wild() == e, true); }
};

static ex f_map(const exprseq &e)
{
	apply_map_function fcn(e[1]);
	return e[0].map(fcn);
}

static ex f_match(const exprseq &e)
{
	lst repl_lst;
	if (e[0].match(e[1], repl_lst))
		return repl_lst;
	else
		return fail();
}

static ex f_normal2(const exprseq &e)
{
	CHECK_ARG(1, numeric, normal);
	return e[0].normal(ex_to<numeric>(e[1]).to_int());
}

static ex f_op(const exprseq &e)
{
	CHECK_ARG(1, numeric, op);
	int n = ex_to<numeric>(e[1]).to_int();
	if (n < 0 || n >= (int)e[0].nops())
		throw(std::out_of_range("second argument to op() is out of range"));
	return e[0].op(n);
}

static ex f_prem(const exprseq &e)
{
	CHECK_ARG(2, symbol, prem);
	return prem(e[0], e[1], ex_to<symbol>(e[2]));
}

static ex f_primpart(const exprseq &e)
{
	CHECK_ARG(1, symbol, primpart);
	return e[0].primpart(ex_to<symbol>(e[1]));
}

static ex f_quo(const exprseq &e)
{
	CHECK_ARG(2, symbol, quo);
	return quo(e[0], e[1], ex_to<symbol>(e[2]));
}

static ex f_rem(const exprseq &e)
{
	CHECK_ARG(2, symbol, rem);
	return rem(e[0], e[1], ex_to<symbol>(e[2]));
}

static ex f_series(const exprseq &e)
{
	CHECK_ARG(2, numeric, series);
	return e[0].series(e[1], ex_to<numeric>(e[2]).to_int());
}

static ex f_sqrfree2(const exprseq &e)
{
	CHECK_ARG(1, lst, sqrfree);
	return sqrfree(e[0], ex_to<lst>(e[1]));
}

static ex f_subs3(const exprseq &e)
{
	CHECK_ARG(1, lst, subs);
	CHECK_ARG(2, lst, subs);
	return e[0].subs(ex_to<lst>(e[1]), ex_to<lst>(e[2]));
}

static ex f_trace(const exprseq &e)
{
	CHECK_ARG(0, matrix, trace);
	return ex_to<matrix>(e[0]).trace();
}

static ex f_transpose(const exprseq &e)
{
	CHECK_ARG(0, matrix, transpose);
	return ex_to<matrix>(e[0]).transpose();
}

static ex f_unassign(const exprseq &e)
{
	CHECK_ARG(0, symbol, unassign);
	const_cast<symbol&>(ex_to<symbol>(e[0])).unassign();
	return e[0];
}

static ex f_unit(const exprseq &e)
{
	CHECK_ARG(1, symbol, unit);
	return e[0].unit(ex_to<symbol>(e[1]));
}

static ex f_dummy(const exprseq &e)
{
	throw(std::logic_error("dummy function called (shouldn't happen)"));
}

// Tables for initializing the "fcns" map and the function help topics
struct fcn_init {
	const char *name;
	const fcn_desc desc;
};

static const fcn_init builtin_fcns[] = {
	{"charpoly", fcn_desc(f_charpoly, 2)},
	{"coeff", fcn_desc(f_coeff, 3)},
	{"collect", fcn_desc(f_collect, 2)},
	{"collect_distributed", fcn_desc(f_collect_distributed, 2)},
	{"content", fcn_desc(f_content, 2)},
	{"decomp_rational", fcn_desc(f_decomp_rational, 2)},
	{"degree", fcn_desc(f_degree, 2)},
	{"denom", fcn_desc(f_denom, 1)},
	{"determinant", fcn_desc(f_determinant, 1)},
	{"diag", fcn_desc(f_diag, 0)},
	{"diff", fcn_desc(f_diff2, 2)},
	{"diff", fcn_desc(f_diff3, 3)},
	{"divide", fcn_desc(f_divide, 2)},
	{"eval", fcn_desc(f_eval1, 1)},
	{"eval", fcn_desc(f_eval2, 2)},
	{"evalf", fcn_desc(f_evalf1, 1)},
	{"evalf", fcn_desc(f_evalf2, 2)},
	{"evalm", fcn_desc(f_evalm, 1)},
	{"expand", fcn_desc(f_expand, 1)},
	{"find", fcn_desc(f_find, 2)},
	{"gcd", fcn_desc(f_gcd, 2)},
	{"has", fcn_desc(f_has, 2)},
	{"inverse", fcn_desc(f_inverse, 1)},
	{"is", fcn_desc(f_is, 1)},
	{"lcm", fcn_desc(f_lcm, 2)},
	{"lcoeff", fcn_desc(f_lcoeff, 2)},
	{"ldegree", fcn_desc(f_ldegree, 2)},
	{"lsolve", fcn_desc(f_lsolve, 2)},
	{"map", fcn_desc(f_map, 2)},
	{"match", fcn_desc(f_match, 2)},
	{"nops", fcn_desc(f_nops, 1)},
	{"normal", fcn_desc(f_normal1, 1)},
	{"normal", fcn_desc(f_normal2, 2)},
	{"numer", fcn_desc(f_numer, 1)},
	{"numer_denom", fcn_desc(f_numer_denom, 1)},
	{"op", fcn_desc(f_op, 2)},
	{"pow", fcn_desc(f_pow, 2)},
	{"prem", fcn_desc(f_prem, 3)},
	{"primpart", fcn_desc(f_primpart, 2)},
	{"quo", fcn_desc(f_quo, 3)},
	{"rem", fcn_desc(f_rem, 3)},
	{"series", fcn_desc(f_series, 3)},
	{"sqrfree", fcn_desc(f_sqrfree1, 1)},
	{"sqrfree", fcn_desc(f_sqrfree2, 2)},
	{"sqrt", fcn_desc(f_sqrt, 1)},
	{"subs", fcn_desc(f_subs2, 2)},
	{"subs", fcn_desc(f_subs3, 3)},
	{"tcoeff", fcn_desc(f_tcoeff, 2)},
	{"time", fcn_desc(f_dummy, 0)},
	{"trace", fcn_desc(f_trace, 1)},
	{"transpose", fcn_desc(f_transpose, 1)},
	{"unassign", fcn_desc(f_unassign, 1)},
	{"unit", fcn_desc(f_unit, 2)},
	{NULL, fcn_desc(f_dummy, 0)}	// End marker
};

struct fcn_help_init {
	const char *name;
	const char *help;
};

static const fcn_help_init builtin_help[] = {
	{"acos", "inverse cosine function"},
	{"acosh", "inverse hyperbolic cosine function"},
	{"asin", "inverse sine function"},
	{"asinh", "inverse hyperbolic sine function"},
	{"atan", "inverse tangent function"},
	{"atan2", "inverse tangent function with two arguments"},
	{"atanh", "inverse hyperbolic tangent function"},
	{"beta", "Beta function"},
	{"binomial", "binomial function"},
	{"cos", "cosine function"},
	{"cosh", "hyperbolic cosine function"},
	{"exp", "exponential function"},
	{"factorial", "factorial function"},
	{"lgamma", "natural logarithm of Gamma function"},
	{"tgamma", "Gamma function"},
	{"log", "natural logarithm"},
	{"psi", "psi function\npsi(x) is the digamma function, psi(n,x) the nth polygamma function"},
	{"sin", "sine function"},
	{"sinh", "hyperbolic sine function"},
	{"tan", "tangent function"},
	{"tanh", "hyperbolic tangent function"},
	{"zeta", "zeta function\nzeta(x) is Riemann's zeta function, zeta(n,x) its nth derivative"},
	{"Li2", "dilogarithm"},
	{"Li3", "trilogarithm"},
	{"Order", "order term function (for truncated power series)"},
	{"Derivative", "inert differential operator"},
	{NULL, NULL}	// End marker
};

#include "ginsh_extensions.h"


/*
 *  Add functions to ginsh
 */

// Functions from fcn_init array
static void insert_fcns(const fcn_init *p)
{
	while (p->name) {
		fcns.insert(make_pair(string(p->name), p->desc));
		p++;
	}
}

static ex f_ginac_function(const exprseq &es, int serial)
{
	return function(serial, es).eval(1);
}

// All registered GiNaC functions
void GiNaC::ginsh_get_ginac_functions(void)
{
	vector<function_options>::const_iterator i = function::registered_functions().begin(), end = function::registered_functions().end();
	unsigned serial = 0;
	while (i != end) {
		fcns.insert(make_pair(i->get_name(), fcn_desc(f_ginac_function, i->get_nparams(), serial)));
		++i;
		serial++;
	}
}


/*
 *  Find a function given a name and number of parameters. Throw exceptions on error.
 */

static fcn_tab::const_iterator find_function(const ex &sym, int req_params)
{
	const string &name = ex_to<symbol>(sym).get_name();
	typedef fcn_tab::const_iterator I;
	pair<I, I> b = fcns.equal_range(name);
	if (b.first == b.second)
		throw(std::logic_error("unknown function '" + name + "'"));
	else {
		for (I i=b.first; i!=b.second; i++)
			if ((i->second.num_params == 0) || (i->second.num_params == req_params))
				return i;
	}
	throw(std::logic_error("invalid number of arguments to " + name + "()"));
}


/*
 *  Insert help strings
 */

// Normal help string
static void insert_help(const char *topic, const char *str)
{
	help.insert(make_pair(string(topic), string(str)));
}

// Help string for functions, automatically generates synopsis
static void insert_fcn_help(const char *name, const char *str)
{
	typedef fcn_tab::const_iterator I;
	pair<I, I> b = fcns.equal_range(name);
	if (b.first != b.second) {
		string help_str = string(name) + "(";
		for (int i=0; i<b.first->second.num_params; i++) {
			if (i)
				help_str += ", ";
			help_str += "expression";
		}
		help_str += ") - ";
		help_str += str;
		help.insert(make_pair(string(name), help_str));
	}
}

// Help strings for functions from fcn_help_init array
static void insert_help(const fcn_help_init *p)
{
	while (p->name) {
		insert_fcn_help(p->name, p->help);
		p++;
	}
}


/*
 *  Print help to cout
 */

// Help for a given topic
static void print_help(const string &topic)
{
	typedef help_tab::const_iterator I;
	pair<I, I> b = help.equal_range(topic);
	if (b.first == b.second)
		cout << "no help for '" << topic << "'\n";
	else {
		for (I i=b.first; i!=b.second; i++)
			cout << i->second << endl;
	}
}

// List of help topics
static void print_help_topics(void)
{
	cout << "Available help topics:\n";
	help_tab::const_iterator i;
	string last_name = string("*");
	int num = 0;
	for (i=help.begin(); i!=help.end(); i++) {
		// Don't print duplicates
		if (i->first != last_name) {
			if (num)
				cout << ", ";
			num++;
			cout << i->first;
			last_name = i->first;
		}
	}
	cout << "\nTo get help for a certain topic, type ?topic\n";
}


/*
 *  Function name completion functions for readline
 */

static char *fcn_generator(const char *text, int state)
{
	static int len;				// Length of word to complete
	static fcn_tab::const_iterator index;	// Iterator to function being currently considered

	// If this is a new word to complete, initialize now
	if (state == 0) {
		index = fcns.begin();
		len = strlen(text);
	}

	// Return the next function which partially matches
	while (index != fcns.end()) {
		const char *fcn_name = index->first.c_str();
		++index;
		if (strncmp(fcn_name, text, len) == 0)
			return strdup(fcn_name);
	}
	return NULL;
}

static char **fcn_completion(const char *text, int start, int end)
{
	if (rl_line_buffer[0] == '!') {
		// For shell commands, revert back to filename completion
		rl_completion_append_character = orig_completion_append_character;
		rl_basic_word_break_characters = orig_basic_word_break_characters;
		rl_completer_word_break_characters = rl_basic_word_break_characters;
#if (GINAC_RL_VERSION_MAJOR < 4) || (GINAC_RL_VERSION_MAJOR == 4 && GINAC_RL_VERSION_MINOR < 2)
		return completion_matches(const_cast<char *>(text), (CPFunction *)filename_completion_function);
#else
		return rl_completion_matches(text, rl_filename_completion_function);
#endif
	} else {
		// Otherwise, complete function names
		rl_completion_append_character = '(';
		rl_basic_word_break_characters = " \t\n\"#$%&'()*+,-./:;<=>?@[\\]^`{|}~";
		rl_completer_word_break_characters = rl_basic_word_break_characters;
#if (GINAC_RL_VERSION_MAJOR < 4) || (GINAC_RL_VERSION_MAJOR == 4 && GINAC_RL_VERSION_MINOR < 2)
		return completion_matches(const_cast<char *>(text), (CPFunction *)fcn_generator);
#else
		return rl_completion_matches(text, fcn_generator);
#endif
	}
}

void greeting(void)
{
    cout << "ginsh - GiNaC Interactive Shell (" << PACKAGE << " V" << VERSION << ")" << endl;
    cout << "  __,  _______  Copyright (C) 1999-2001 Johannes Gutenberg University Mainz,\n"
         << " (__) *       | Germany.  This is free software with ABSOLUTELY NO WARRANTY.\n"
         << "  ._) i N a C | You are welcome to redistribute it under certain conditions.\n"
         << "<-------------' For details type `warranty;'.\n" << endl;
    cout << "Type ?? for a list of help topics." << endl;
}

/*
 *  Main program
 */

int main(int argc, char **argv)
{
	// Print banner in interactive mode
	if (isatty(0)) 
		greeting();

	// Init function table
	insert_fcns(builtin_fcns);
	insert_fcns(extended_fcns);
	ginsh_get_ginac_functions();

	// Init help for operators (automatically generated from man page)
	insert_help("operators", "Operators in falling order of precedence:");
#include "ginsh_op_help.h"

	// Init help for built-in functions (automatically generated from man page)
#include "ginsh_fcn_help.h"

	// Help for GiNaC functions is added manually
	insert_help(builtin_help);
	insert_help(extended_help);

	// Init readline completer
	rl_readline_name = argv[0];
#if (GINAC_RL_VERSION_MAJOR < 4) || (GINAC_RL_VERSION_MAJOR == 4 && GINAC_RL_VERSION_MINOR < 2)
	rl_attempted_completion_function = (CPPFunction *)fcn_completion;
#else
	rl_attempted_completion_function = fcn_completion;
#endif
	orig_completion_append_character = rl_completion_append_character;
	orig_basic_word_break_characters = rl_basic_word_break_characters;

	// Init input file list, open first file
	num_files = argc - 1;
	file_list = argv + 1;
	if (num_files) {
		yyin = fopen(*file_list, "r");
		if (yyin == NULL) {
			cerr << "Can't open " << *file_list << endl;
			exit(1);
		}
		num_files--;
		file_list++;
	}

	// Parse input, catch all remaining exceptions
	int result;
again:	try {
		result = yyparse();
	} catch (exception &e) {
		cerr << e.what() << endl;
		goto again;
	}
	return result;
}
