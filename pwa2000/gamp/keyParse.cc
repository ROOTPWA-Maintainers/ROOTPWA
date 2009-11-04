/* A Bison parser, made by GNU Bison 1.875d.  */

/* Skeleton parser for Yacc-like parsing with Bison,
   Copyright (C) 1984, 1989, 1990, 2000, 2001, 2002, 2003, 2004 Free Software Foundation, Inc.

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

/* Written by Richard Stallman by simplifying the original so called
   ``semantic'' parser.  */

/* All symbols defined below should begin with yy or YY, to avoid
   infringing on user name space.  This should be done even for local
   variables, as they might otherwise be expanded by user macros.
   There are some unavoidable exceptions within include files to
   define necessary library symbols; they are noted "INFRINGES ON
   USER NAME SPACE" below.  */

/* Identify Bison output.  */
#define YYBISON 1

/* Skeleton name.  */
#define YYSKELETON_NAME "yacc.c"

/* Pure parsers.  */
#define YYPURE 0

/* Using locations.  */
#define YYLSP_NEEDED 0

/* If NAME_PREFIX is specified substitute the variables and functions
   names.  */
#define yyparse keyparse
#define yylex   keylex
#define yyerror keyerror
#define yylval  keylval
#define yychar  keychar
#define yydebug keydebug
#define yynerrs keynerrs


/* Tokens.  */
#ifndef YYTOKENTYPE
# define YYTOKENTYPE
   /* Put the tokens into the symbol table, so that GDB and other debuggers
      know about them.  */
   enum yytokentype {
     INT = 258,
     FLOAT = 259,
     COMPLEX = 260,
     SIGN = 261,
     STRING = 262,
     PARTICLENAME = 263,
     DEBUG = 264,
     CHANNEL = 265,
     MODE = 266,
     MASSDEP = 267,
     HELICITY = 268
   };
#endif
#define INT 258
#define FLOAT 259
#define COMPLEX 260
#define SIGN 261
#define STRING 262
#define PARTICLENAME 263
#define DEBUG 264
#define CHANNEL 265
#define MODE 266
#define MASSDEP 267
#define HELICITY 268




/* Copy the first part of user declarations.  */
#line 1 "keyParse.yy"

#include <iostream>
#include <particle.h>
#include <wave.h>
#include <keyfile.h>
#include <massDep.h>
#define stoi(x) strcmp((x),"+")?-1:1

using namespace std;

int yylex(void);
void yyerror(char* s);
#define YYDEBUG 1

int nwave;
int debug = 0;
string mode;

wave wv;
extern int lineno;
extern particleDataTable PDGtable;
extern event e;
complex<double> amp;
string t_part_init;
particle* t_part_final;



/* Enabling traces.  */
#ifndef YYDEBUG
# define YYDEBUG 0
#endif

/* Enabling verbose error messages.  */
#ifdef YYERROR_VERBOSE
# undef YYERROR_VERBOSE
# define YYERROR_VERBOSE 1
#else
# define YYERROR_VERBOSE 0
#endif

#if ! defined (YYSTYPE) && ! defined (YYSTYPE_IS_DECLARED)
#line 29 "keyParse.yy"
typedef union YYSTYPE {
        int num;
        std::complex<double>* Cnum;
        double Fnum;
        char string[100];
        decay* Decay;
        particle* Particle;
} YYSTYPE;
/* Line 191 of yacc.c.  */
#line 148 "keyParse.cc"
# define yystype YYSTYPE /* obsolescent; will be withdrawn */
# define YYSTYPE_IS_DECLARED 1
# define YYSTYPE_IS_TRIVIAL 1
#endif



/* Copy the second part of user declarations.  */


/* Line 214 of yacc.c.  */
#line 160 "keyParse.cc"

#if ! defined (yyoverflow) || YYERROR_VERBOSE

# ifndef YYFREE
#  define YYFREE free
# endif
# ifndef YYMALLOC
#  define YYMALLOC malloc
# endif

/* The parser invokes alloca or malloc; define the necessary symbols.  */

# ifdef YYSTACK_USE_ALLOCA
#  if YYSTACK_USE_ALLOCA
#   define YYSTACK_ALLOC alloca
#  endif
# else
#  if defined (alloca) || defined (_ALLOCA_H)
#   define YYSTACK_ALLOC alloca
#  else
#   ifdef __GNUC__
#    define YYSTACK_ALLOC __builtin_alloca
#   endif
#  endif
# endif

# ifdef YYSTACK_ALLOC
   /* Pacify GCC's `empty if-body' warning. */
#  define YYSTACK_FREE(Ptr) do { /* empty */; } while (0)
# else
#  if defined (__STDC__) || defined (__cplusplus)
#   include <stdlib.h> /* INFRINGES ON USER NAME SPACE */
#   define YYSIZE_T size_t
#  endif
#  define YYSTACK_ALLOC YYMALLOC
#  define YYSTACK_FREE YYFREE
# endif
#endif /* ! defined (yyoverflow) || YYERROR_VERBOSE */


#if (! defined (yyoverflow) \
     && (! defined (__cplusplus) \
	 || (defined (YYSTYPE_IS_TRIVIAL) && YYSTYPE_IS_TRIVIAL)))

/* A type that is properly aligned for any stack member.  */
union yyalloc
{
  short int yyss;
  YYSTYPE yyvs;
  };

/* The size of the maximum gap between one aligned stack and the next.  */
# define YYSTACK_GAP_MAXIMUM (sizeof (union yyalloc) - 1)

/* The size of an array large to enough to hold all stacks, each with
   N elements.  */
# define YYSTACK_BYTES(N) \
     ((N) * (sizeof (short int) + sizeof (YYSTYPE))			\
      + YYSTACK_GAP_MAXIMUM)

/* Copy COUNT objects from FROM to TO.  The source and destination do
   not overlap.  */
# ifndef YYCOPY
#  if defined (__GNUC__) && 1 < __GNUC__
#   define YYCOPY(To, From, Count) \
      __builtin_memcpy (To, From, (Count) * sizeof (*(From)))
#  else
#   define YYCOPY(To, From, Count)		\
      do					\
	{					\
	  register YYSIZE_T yyi;		\
	  for (yyi = 0; yyi < (Count); yyi++)	\
	    (To)[yyi] = (From)[yyi];		\
	}					\
      while (0)
#  endif
# endif

/* Relocate STACK from its old location to the new one.  The
   local variables YYSIZE and YYSTACKSIZE give the old and new number of
   elements in the stack, and YYPTR gives the new location of the
   stack.  Advance YYPTR to a properly aligned location for the next
   stack.  */
# define YYSTACK_RELOCATE(Stack)					\
    do									\
      {									\
	YYSIZE_T yynewbytes;						\
	YYCOPY (&yyptr->Stack, Stack, yysize);				\
	Stack = &yyptr->Stack;						\
	yynewbytes = yystacksize * sizeof (*Stack) + YYSTACK_GAP_MAXIMUM; \
	yyptr += yynewbytes / sizeof (*yyptr);				\
      }									\
    while (0)

#endif

#if defined (__STDC__) || defined (__cplusplus)
   typedef signed char yysigned_char;
#else
   typedef short int yysigned_char;
#endif

/* YYFINAL -- State number of the termination state. */
#define YYFINAL  2
/* YYLAST -- Last index in YYTABLE.  */
#define YYLAST   91

/* YYNTOKENS -- Number of terminals. */
#define YYNTOKENS  27
/* YYNNTS -- Number of nonterminals. */
#define YYNNTS  14
/* YYNRULES -- Number of rules. */
#define YYNRULES  37
/* YYNRULES -- Number of states. */
#define YYNSTATES  89

/* YYTRANSLATE(YYLEX) -- Bison symbol number corresponding to YYLEX.  */
#define YYUNDEFTOK  2
#define YYMAXUTOK   268

#define YYTRANSLATE(YYX) 						\
  ((unsigned int) (YYX) <= YYMAXUTOK ? yytranslate[YYX] : YYUNDEFTOK)

/* YYTRANSLATE[YYLEX] -- Bison symbol number corresponding to YYLEX.  */
static const unsigned char yytranslate[] =
{
       0,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
      20,    21,    16,    14,    19,    15,     2,     2,    26,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,    18,
       2,    17,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,    24,     2,    25,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,    22,     2,    23,     2,     2,     2,     2,
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
       2,     2,     2,     2,     2,     2,     1,     2,     3,     4,
       5,     6,     7,     8,     9,    10,    11,    12,    13
};

#if YYDEBUG
/* YYPRHS[YYN] -- Index of the first RHS symbol of rule number YYN in
   YYRHS.  */
static const unsigned char yyprhs[] =
{
       0,     0,     3,     4,     7,    12,    17,    22,    29,    32,
      34,    38,    42,    46,    50,    54,    57,    59,    61,    65,
      69,    75,    83,    91,   102,   109,   116,   118,   121,   127,
     129,   134,   136,   141,   143,   146,   149,   152
};

/* YYRHS -- A `-1'-separated list of the rules' RHS. */
static const yysigned_char yyrhs[] =
{
      28,     0,    -1,    -1,    28,    29,    -1,     9,    17,     3,
      18,    -1,    10,    17,     7,    18,    -1,    11,    17,     7,
      18,    -1,     7,    17,     7,    19,    38,    18,    -1,    30,
      18,    -1,    31,    -1,    20,    30,    21,    -1,    30,    14,
      30,    -1,    30,    15,    30,    -1,     4,    16,    30,    -1,
      40,    16,    30,    -1,    32,    35,    -1,    35,    -1,    33,
      -1,    34,    34,    34,    -1,     7,    17,     3,    -1,    22,
      36,    36,     3,    23,    -1,    22,    36,    36,     7,    17,
       3,    23,    -1,    22,    36,    36,     7,    17,     4,    23,
      -1,    22,    36,    36,     7,    17,     3,     7,    17,     3,
      23,    -1,    22,    36,    36,     3,     3,    23,    -1,    22,
      36,    36,    36,     3,    23,    -1,    37,    -1,    37,    35,
      -1,    37,    35,    12,    17,     7,    -1,    38,    -1,    38,
      13,    17,     3,    -1,    39,    -1,    39,    24,     3,    25,
      -1,     8,    -1,     8,    14,    -1,     8,    15,    -1,     8,
      26,    -1,    20,     4,    19,     4,    21,    -1
};

/* YYRLINE[YYN] -- source line where rule number YYN was defined.  */
static const unsigned short int yyrline[] =
{
       0,    62,    62,    63,    66,    69,    72,    75,    79,    94,
      98,   106,   115,   124,   132,   142,   164,   210,   214,   218,
     225,   235,   253,   270,   294,   304,   318,   321,   328,   354,
     357,   363,   366,   372,   376,   380,   384,   390
};
#endif

#if YYDEBUG || YYERROR_VERBOSE
/* YYTNME[SYMBOL-NUM] -- String name of the symbol SYMBOL-NUM.
   First, the terminals, then, starting at YYNTOKENS, nonterminals. */
static const char *const yytname[] =
{
  "$end", "error", "$undefined", "INT", "FLOAT", "COMPLEX", "SIGN",
  "STRING", "PARTICLENAME", "DEBUG", "CHANNEL", "MODE", "MASSDEP",
  "HELICITY", "'+'", "'-'", "'*'", "'='", "';'", "','", "'('", "')'",
  "'{'", "'}'", "'['", "']'", "'0'", "$accept", "input", "statement",
  "waveexp", "wave", "resonance", "quantum_nums", "quantum_num", "decay",
  "particle", "pstate", "particleID", "particleCharge", "complex", 0
};
#endif

# ifdef YYPRINT
/* YYTOKNUM[YYLEX-NUM] -- Internal token number corresponding to
   token YYLEX-NUM.  */
static const unsigned short int yytoknum[] =
{
       0,   256,   257,   258,   259,   260,   261,   262,   263,   264,
     265,   266,   267,   268,    43,    45,    42,    61,    59,    44,
      40,    41,   123,   125,    91,    93,    48
};
# endif

/* YYR1[YYN] -- Symbol number of symbol that rule YYN derives.  */
static const unsigned char yyr1[] =
{
       0,    27,    28,    28,    29,    29,    29,    29,    29,    30,
      30,    30,    30,    30,    30,    31,    31,    32,    33,    34,
      35,    35,    35,    35,    35,    35,    36,    36,    36,    37,
      37,    38,    38,    39,    39,    39,    39,    40
};

/* YYR2[YYN] -- Number of symbols composing right hand side of rule YYN.  */
static const unsigned char yyr2[] =
{
       0,     2,     0,     2,     4,     4,     4,     6,     2,     1,
       3,     3,     3,     3,     3,     2,     1,     1,     3,     3,
       5,     7,     7,    10,     6,     6,     1,     2,     5,     1,
       4,     1,     4,     1,     2,     2,     2,     5
};

/* YYDEFACT[STATE-NAME] -- Default rule to reduce with in state
   STATE-NUM when YYTABLE doesn't specify something else to do.  Zero
   means the default is an error.  */
static const unsigned char yydefact[] =
{
       2,     0,     1,     0,     0,     0,     0,     0,     0,     0,
       3,     0,     9,     0,    17,     0,    16,     0,     0,     0,
       0,     0,     0,     0,     0,     0,    33,     0,    26,    29,
      31,     0,     0,     8,    15,     0,     0,    13,    19,     0,
       0,     0,     0,     0,     0,    10,    34,    35,    36,     0,
      27,     0,     0,    11,    12,    18,    14,     0,     4,     5,
       6,     0,     0,     0,     0,     0,     0,     0,     0,    37,
       0,    20,     0,     0,     0,    30,    32,     7,    24,     0,
       0,    25,    28,     0,    21,    22,     0,     0,    23
};

/* YYDEFGOTO[NTERM-NUM]. */
static const yysigned_char yydefgoto[] =
{
      -1,     1,    10,    11,    12,    13,    14,    15,    16,    27,
      28,    29,    30,    17
};

/* YYPACT[STATE-NUM] -- Index in YYTABLE of the portion describing
   STATE-NUM.  */
#define YYPACT_NINF -24
static const yysigned_char yypact[] =
{
     -24,     5,   -24,   -10,     1,    24,    38,    39,    13,    46,
     -24,    -7,   -24,    -9,   -24,    23,   -24,    41,    25,    37,
      55,    52,    53,    34,    44,    28,    22,    46,    -9,    49,
      40,    25,    25,   -24,   -24,    23,    25,   -24,   -24,    47,
      45,    50,    51,    61,    64,   -24,   -24,   -24,   -24,    31,
      58,    54,    69,   -24,   -24,   -24,   -24,    46,   -24,   -24,
     -24,    56,    -2,    57,    70,    59,    72,    60,    62,   -24,
      63,   -24,    48,    65,    71,   -24,   -24,   -24,   -24,    -4,
      66,   -24,   -24,    67,   -24,   -24,    76,    68,   -24
};

/* YYPGOTO[NTERM-NUM].  */
static const yysigned_char yypgoto[] =
{
     -24,   -24,   -24,    -8,   -24,   -24,   -24,   -13,    18,   -23,
     -24,    26,   -24,   -24
};

/* YYTABLE[YYPACT[STATE-NUM]].  What to do in state STATE-NUM.  If
   positive, shift that token.  If negative, reduce the rule which
   number is the opposite.  If zero, do what YYDEFACT says.
   If YYTABLE_NINF, syntax error.  */
#define YYTABLE_NINF -1
static const unsigned char yytable[] =
{
      25,    70,    35,    83,    49,     2,    18,    31,    32,     3,
      37,    33,     4,     9,     5,     6,     7,    23,    19,    84,
      24,    71,    55,    53,    54,     8,    64,     9,    56,     3,
      24,    34,    24,     8,    62,     9,    46,    47,    63,    26,
      38,    20,    31,    32,    39,     8,    50,     9,    48,    45,
      18,    79,    80,    43,    26,    21,    22,    36,    40,    41,
      42,    44,    51,    58,    52,    61,    57,    38,    59,    60,
      65,    66,    67,    73,    72,    75,    74,    69,    82,    87,
      77,     0,     0,    68,    86,    76,    78,     0,    81,    85,
       0,    88
};

static const yysigned_char yycheck[] =
{
       8,     3,    15,     7,    27,     0,    16,    14,    15,     4,
      18,    18,     7,    22,     9,    10,    11,     4,    17,    23,
       7,    23,    35,    31,    32,    20,    49,    22,    36,     4,
       7,    13,     7,    20,     3,    22,    14,    15,     7,     8,
       3,    17,    14,    15,     7,    20,    28,    22,    26,    21,
      16,     3,     4,    19,     8,    17,    17,    16,     3,     7,
       7,    17,    13,    18,    24,     4,    19,     3,    18,    18,
      12,    17,     3,     3,    17,     3,    17,    21,     7,     3,
      18,    -1,    -1,    57,    17,    25,    23,    -1,    23,    23,
      -1,    23
};

/* YYSTOS[STATE-NUM] -- The (internal number of the) accessing
   symbol of state STATE-NUM.  */
static const unsigned char yystos[] =
{
       0,    28,     0,     4,     7,     9,    10,    11,    20,    22,
      29,    30,    31,    32,    33,    34,    35,    40,    16,    17,
      17,    17,    17,     4,     7,    30,     8,    36,    37,    38,
      39,    14,    15,    18,    35,    34,    16,    30,     3,     7,
       3,     7,     7,    19,    17,    21,    14,    15,    26,    36,
      35,    13,    24,    30,    30,    34,    30,    19,    18,    18,
      18,     4,     3,     7,    36,    12,    17,     3,    38,    21,
       3,    23,    17,     3,    17,     3,    25,    18,    23,     3,
       4,    23,     7,     7,    23,    23,    17,     3,    23
};

#if ! defined (YYSIZE_T) && defined (__SIZE_TYPE__)
# define YYSIZE_T __SIZE_TYPE__
#endif
#if ! defined (YYSIZE_T) && defined (size_t)
# define YYSIZE_T size_t
#endif
#if ! defined (YYSIZE_T)
# if defined (__STDC__) || defined (__cplusplus)
#  include <stddef.h> /* INFRINGES ON USER NAME SPACE */
#  define YYSIZE_T size_t
# endif
#endif
#if ! defined (YYSIZE_T)
# define YYSIZE_T unsigned int
#endif

#define yyerrok		(yyerrstatus = 0)
#define yyclearin	(yychar = YYEMPTY)
#define YYEMPTY		(-2)
#define YYEOF		0

#define YYACCEPT	goto yyacceptlab
#define YYABORT		goto yyabortlab
#define YYERROR		goto yyerrorlab


/* Like YYERROR except do call yyerror.  This remains here temporarily
   to ease the transition to the new meaning of YYERROR, for GCC.
   Once GCC version 2 has supplanted version 1, this can go.  */

#define YYFAIL		goto yyerrlab

#define YYRECOVERING()  (!!yyerrstatus)

#define YYBACKUP(Token, Value)					\
do								\
  if (yychar == YYEMPTY && yylen == 1)				\
    {								\
      yychar = (Token);						\
      yylval = (Value);						\
      yytoken = YYTRANSLATE (yychar);				\
      YYPOPSTACK;						\
      goto yybackup;						\
    }								\
  else								\
    { 								\
      yyerror ("syntax error: cannot back up");\
      YYERROR;							\
    }								\
while (0)

#define YYTERROR	1
#define YYERRCODE	256

/* YYLLOC_DEFAULT -- Compute the default location (before the actions
   are run).  */

#ifndef YYLLOC_DEFAULT
# define YYLLOC_DEFAULT(Current, Rhs, N)		\
   ((Current).first_line   = (Rhs)[1].first_line,	\
    (Current).first_column = (Rhs)[1].first_column,	\
    (Current).last_line    = (Rhs)[N].last_line,	\
    (Current).last_column  = (Rhs)[N].last_column)
#endif

/* YYLEX -- calling `yylex' with the right arguments.  */

#ifdef YYLEX_PARAM
# define YYLEX yylex (YYLEX_PARAM)
#else
# define YYLEX yylex ()
#endif

/* Enable debugging if requested.  */
#if YYDEBUG

# ifndef YYFPRINTF
#  include <stdio.h> /* INFRINGES ON USER NAME SPACE */
#  define YYFPRINTF fprintf
# endif

# define YYDPRINTF(Args)			\
do {						\
  if (yydebug)					\
    YYFPRINTF Args;				\
} while (0)

# define YYDSYMPRINT(Args)			\
do {						\
  if (yydebug)					\
    yysymprint Args;				\
} while (0)

# define YYDSYMPRINTF(Title, Token, Value, Location)		\
do {								\
  if (yydebug)							\
    {								\
      YYFPRINTF (stderr, "%s ", Title);				\
      yysymprint (stderr, 					\
                  Token, Value);	\
      YYFPRINTF (stderr, "\n");					\
    }								\
} while (0)

/*------------------------------------------------------------------.
| yy_stack_print -- Print the state stack from its BOTTOM up to its |
| TOP (included).                                                   |
`------------------------------------------------------------------*/

#if defined (__STDC__) || defined (__cplusplus)
static void
yy_stack_print (short int *bottom, short int *top)
#else
static void
yy_stack_print (bottom, top)
    short int *bottom;
    short int *top;
#endif
{
  YYFPRINTF (stderr, "Stack now");
  for (/* Nothing. */; bottom <= top; ++bottom)
    YYFPRINTF (stderr, " %d", *bottom);
  YYFPRINTF (stderr, "\n");
}

# define YY_STACK_PRINT(Bottom, Top)				\
do {								\
  if (yydebug)							\
    yy_stack_print ((Bottom), (Top));				\
} while (0)


/*------------------------------------------------.
| Report that the YYRULE is going to be reduced.  |
`------------------------------------------------*/

#if defined (__STDC__) || defined (__cplusplus)
static void
yy_reduce_print (int yyrule)
#else
static void
yy_reduce_print (yyrule)
    int yyrule;
#endif
{
  int yyi;
  unsigned int yylno = yyrline[yyrule];
  YYFPRINTF (stderr, "Reducing stack by rule %d (line %u), ",
             yyrule - 1, yylno);
  /* Print the symbols being reduced, and their result.  */
  for (yyi = yyprhs[yyrule]; 0 <= yyrhs[yyi]; yyi++)
    YYFPRINTF (stderr, "%s ", yytname [yyrhs[yyi]]);
  YYFPRINTF (stderr, "-> %s\n", yytname [yyr1[yyrule]]);
}

# define YY_REDUCE_PRINT(Rule)		\
do {					\
  if (yydebug)				\
    yy_reduce_print (Rule);		\
} while (0)

/* Nonzero means print parse trace.  It is left uninitialized so that
   multiple parsers can coexist.  */
int yydebug;
#else /* !YYDEBUG */
# define YYDPRINTF(Args)
# define YYDSYMPRINT(Args)
# define YYDSYMPRINTF(Title, Token, Value, Location)
# define YY_STACK_PRINT(Bottom, Top)
# define YY_REDUCE_PRINT(Rule)
#endif /* !YYDEBUG */


/* YYINITDEPTH -- initial size of the parser's stacks.  */
#ifndef	YYINITDEPTH
# define YYINITDEPTH 200
#endif

/* YYMAXDEPTH -- maximum size the stacks can grow to (effective only
   if the built-in stack extension method is used).

   Do not make this value too large; the results are undefined if
   SIZE_MAX < YYSTACK_BYTES (YYMAXDEPTH)
   evaluated with infinite-precision integer arithmetic.  */

#if defined (YYMAXDEPTH) && YYMAXDEPTH == 0
# undef YYMAXDEPTH
#endif

#ifndef YYMAXDEPTH
# define YYMAXDEPTH 10000
#endif



#if YYERROR_VERBOSE

# ifndef yystrlen
#  if defined (__GLIBC__) && defined (_STRING_H)
#   define yystrlen strlen
#  else
/* Return the length of YYSTR.  */
static YYSIZE_T
#   if defined (__STDC__) || defined (__cplusplus)
yystrlen (const char *yystr)
#   else
yystrlen (yystr)
     const char *yystr;
#   endif
{
  register const char *yys = yystr;

  while (*yys++ != '\0')
    continue;

  return yys - yystr - 1;
}
#  endif
# endif

# ifndef yystpcpy
#  if defined (__GLIBC__) && defined (_STRING_H) && defined (_GNU_SOURCE)
#   define yystpcpy stpcpy
#  else
/* Copy YYSRC to YYDEST, returning the address of the terminating '\0' in
   YYDEST.  */
static char *
#   if defined (__STDC__) || defined (__cplusplus)
yystpcpy (char *yydest, const char *yysrc)
#   else
yystpcpy (yydest, yysrc)
     char *yydest;
     const char *yysrc;
#   endif
{
  register char *yyd = yydest;
  register const char *yys = yysrc;

  while ((*yyd++ = *yys++) != '\0')
    continue;

  return yyd - 1;
}
#  endif
# endif

#endif /* !YYERROR_VERBOSE */



#if YYDEBUG
/*--------------------------------.
| Print this symbol on YYOUTPUT.  |
`--------------------------------*/

#if defined (__STDC__) || defined (__cplusplus)
static void
yysymprint (FILE *yyoutput, int yytype, YYSTYPE *yyvaluep)
#else
static void
yysymprint (yyoutput, yytype, yyvaluep)
    FILE *yyoutput;
    int yytype;
    YYSTYPE *yyvaluep;
#endif
{
  /* Pacify ``unused variable'' warnings.  */
  (void) yyvaluep;

  if (yytype < YYNTOKENS)
    {
      YYFPRINTF (yyoutput, "token %s (", yytname[yytype]);
# ifdef YYPRINT
      YYPRINT (yyoutput, yytoknum[yytype], *yyvaluep);
# endif
    }
  else
    YYFPRINTF (yyoutput, "nterm %s (", yytname[yytype]);

  switch (yytype)
    {
      default:
        break;
    }
  YYFPRINTF (yyoutput, ")");
}

#endif /* ! YYDEBUG */
/*-----------------------------------------------.
| Release the memory associated to this symbol.  |
`-----------------------------------------------*/

#if defined (__STDC__) || defined (__cplusplus)
static void
yydestruct (int yytype, YYSTYPE *yyvaluep)
#else
static void
yydestruct (yytype, yyvaluep)
    int yytype;
    YYSTYPE *yyvaluep;
#endif
{
  /* Pacify ``unused variable'' warnings.  */
  (void) yyvaluep;

  switch (yytype)
    {

      default:
        break;
    }
}


/* Prevent warnings from -Wmissing-prototypes.  */

#ifdef YYPARSE_PARAM
# if defined (__STDC__) || defined (__cplusplus)
int yyparse (void *YYPARSE_PARAM);
# else
int yyparse ();
# endif
#else /* ! YYPARSE_PARAM */
#if defined (__STDC__) || defined (__cplusplus)
int yyparse (void);
#else
int yyparse ();
#endif
#endif /* ! YYPARSE_PARAM */



/* The lookahead symbol.  */
int yychar;

/* The semantic value of the lookahead symbol.  */
YYSTYPE yylval;

/* Number of syntax errors so far.  */
int yynerrs;



/*----------.
| yyparse.  |
`----------*/

#ifdef YYPARSE_PARAM
# if defined (__STDC__) || defined (__cplusplus)
int yyparse (void *YYPARSE_PARAM)
# else
int yyparse (YYPARSE_PARAM)
  void *YYPARSE_PARAM;
# endif
#else /* ! YYPARSE_PARAM */
#if defined (__STDC__) || defined (__cplusplus)
int
yyparse (std::complex<double>& result)
#else
int
yyparse (std::complex<double>& result)

#endif
#endif
{
  
  register int yystate;
  register int yyn;
  int yyresult;
  /* Number of tokens to shift before error messages enabled.  */
  int yyerrstatus;
  /* Lookahead token as an internal (translated) token number.  */
  int yytoken = 0;

  /* Three stacks and their tools:
     `yyss': related to states,
     `yyvs': related to semantic values,
     `yyls': related to locations.

     Refer to the stacks thru separate pointers, to allow yyoverflow
     to reallocate them elsewhere.  */

  /* The state stack.  */
  short int yyssa[YYINITDEPTH];
  short int *yyss = yyssa;
  register short int *yyssp;

  /* The semantic value stack.  */
  YYSTYPE yyvsa[YYINITDEPTH];
  YYSTYPE *yyvs = yyvsa;
  register YYSTYPE *yyvsp;



#define YYPOPSTACK   (yyvsp--, yyssp--)

  YYSIZE_T yystacksize = YYINITDEPTH;

  /* The variables used to return semantic value and location from the
     action routines.  */
  YYSTYPE yyval;


  /* When reducing, the number of symbols on the RHS of the reduced
     rule.  */
  int yylen;

  YYDPRINTF ((stderr, "Starting parse\n"));

  yystate = 0;
  yyerrstatus = 0;
  yynerrs = 0;
  yychar = YYEMPTY;		/* Cause a token to be read.  */

  /* Initialize stack pointers.
     Waste one element of value and location stack
     so that they stay on the same level as the state stack.
     The wasted elements are never initialized.  */

  yyssp = yyss;
  yyvsp = yyvs;


  goto yysetstate;

/*------------------------------------------------------------.
| yynewstate -- Push a new state, which is found in yystate.  |
`------------------------------------------------------------*/
 yynewstate:
  /* In all cases, when you get here, the value and location stacks
     have just been pushed. so pushing a state here evens the stacks.
     */
  yyssp++;

 yysetstate:
  *yyssp = yystate;

  if (yyss + yystacksize - 1 <= yyssp)
    {
      /* Get the current used size of the three stacks, in elements.  */
      YYSIZE_T yysize = yyssp - yyss + 1;

#ifdef yyoverflow
      {
	/* Give user a chance to reallocate the stack. Use copies of
	   these so that the &'s don't force the real ones into
	   memory.  */
	YYSTYPE *yyvs1 = yyvs;
	short int *yyss1 = yyss;


	/* Each stack pointer address is followed by the size of the
	   data in use in that stack, in bytes.  This used to be a
	   conditional around just the two extra args, but that might
	   be undefined if yyoverflow is a macro.  */
	yyoverflow ("parser stack overflow",
		    &yyss1, yysize * sizeof (*yyssp),
		    &yyvs1, yysize * sizeof (*yyvsp),

		    &yystacksize);

	yyss = yyss1;
	yyvs = yyvs1;
      }
#else /* no yyoverflow */
# ifndef YYSTACK_RELOCATE
      goto yyoverflowlab;
# else
      /* Extend the stack our own way.  */
      if (YYMAXDEPTH <= yystacksize)
	goto yyoverflowlab;
      yystacksize *= 2;
      if (YYMAXDEPTH < yystacksize)
	yystacksize = YYMAXDEPTH;

      {
	short int *yyss1 = yyss;
	union yyalloc *yyptr =
	  (union yyalloc *) YYSTACK_ALLOC (YYSTACK_BYTES (yystacksize));
	if (! yyptr)
	  goto yyoverflowlab;
	YYSTACK_RELOCATE (yyss);
	YYSTACK_RELOCATE (yyvs);

#  undef YYSTACK_RELOCATE
	if (yyss1 != yyssa)
	  YYSTACK_FREE (yyss1);
      }
# endif
#endif /* no yyoverflow */

      yyssp = yyss + yysize - 1;
      yyvsp = yyvs + yysize - 1;


      YYDPRINTF ((stderr, "Stack size increased to %lu\n",
		  (unsigned long int) yystacksize));

      if (yyss + yystacksize - 1 <= yyssp)
	YYABORT;
    }

  YYDPRINTF ((stderr, "Entering state %d\n", yystate));

  goto yybackup;

/*-----------.
| yybackup.  |
`-----------*/
yybackup:

/* Do appropriate processing given the current state.  */
/* Read a lookahead token if we need one and don't already have one.  */
/* yyresume: */

  /* First try to decide what to do without reference to lookahead token.  */

  yyn = yypact[yystate];
  if (yyn == YYPACT_NINF)
    goto yydefault;

  /* Not known => get a lookahead token if don't already have one.  */

  /* YYCHAR is either YYEMPTY or YYEOF or a valid lookahead symbol.  */
  if (yychar == YYEMPTY)
    {
      YYDPRINTF ((stderr, "Reading a token: "));
      yychar = YYLEX;
    }

  if (yychar <= YYEOF)
    {
      yychar = yytoken = YYEOF;
      YYDPRINTF ((stderr, "Now at end of input.\n"));
    }
  else
    {
      yytoken = YYTRANSLATE (yychar);
      YYDSYMPRINTF ("Next token is", yytoken, &yylval, &yylloc);
    }

  /* If the proper action on seeing token YYTOKEN is to reduce or to
     detect an error, take that action.  */
  yyn += yytoken;
  if (yyn < 0 || YYLAST < yyn || yycheck[yyn] != yytoken)
    goto yydefault;
  yyn = yytable[yyn];
  if (yyn <= 0)
    {
      if (yyn == 0 || yyn == YYTABLE_NINF)
	goto yyerrlab;
      yyn = -yyn;
      goto yyreduce;
    }

  if (yyn == YYFINAL)
    YYACCEPT;

  /* Shift the lookahead token.  */
  YYDPRINTF ((stderr, "Shifting token %s, ", yytname[yytoken]));

  /* Discard the token being shifted unless it is eof.  */
  if (yychar != YYEOF)
    yychar = YYEMPTY;

  *++yyvsp = yylval;


  /* Count tokens shifted since error; after three, turn off error
     status.  */
  if (yyerrstatus)
    yyerrstatus--;

  yystate = yyn;
  goto yynewstate;


/*-----------------------------------------------------------.
| yydefault -- do the default action for the current state.  |
`-----------------------------------------------------------*/
yydefault:
  yyn = yydefact[yystate];
  if (yyn == 0)
    goto yyerrlab;
  goto yyreduce;


/*-----------------------------.
| yyreduce -- Do a reduction.  |
`-----------------------------*/
yyreduce:
  /* yyn is the number of a rule to reduce with.  */
  yylen = yyr2[yyn];

  /* If YYLEN is nonzero, implement the default value of the action:
     `$$ = $1'.

     Otherwise, the following line sets YYVAL to garbage.
     This behavior is undocumented and Bison
     users should not rely upon it.  Assigning to YYVAL
     unconditionally makes the parser a bit smaller, and it avoids a
     GCC warning that YYVAL may be used uninitialized.  */
  yyval = yyvsp[1-yylen];


  YY_REDUCE_PRINT (yyn);
  switch (yyn)
    {
        case 4:
#line 66 "keyParse.yy"
    {
                debug = yyvsp[-1].num;
        }
    break;

  case 5:
#line 69 "keyParse.yy"
    {
                wv.channel(yyvsp[-1].string);
        }
    break;

  case 6:
#line 72 "keyParse.yy"
    {
                mode = yyvsp[-1].string;
        }
    break;

  case 7:
#line 75 "keyParse.yy"
    {
                t_part_init = yyvsp[-3].string;
                t_part_final = yyvsp[-1].Particle;
        }
    break;

  case 8:
#line 79 "keyParse.yy"
    {
      result=*yyvsp[-1].Cnum;
	if (mode == "binary") {
                        cout.write((char*) yyvsp[-1].Cnum,sizeof(complex<double>));
		}
                else {
                        cout << "Mass = " << ~(wv.get4P()) << "\t";
                        if ( wv.channel() == "t" ) {
                                cout << "t = " << (e.beam().get4P()-wv.get4P()).lenSq() << "\t";
                        }
                        cout << "Amp = " <<  *yyvsp[-1].Cnum << endl;
                }
                delete yyvsp[-1].Cnum;
        }
    break;

  case 9:
#line 94 "keyParse.yy"
    {
                yyval.Cnum = new complex<double>(*yyvsp[0].Cnum);
                delete(yyvsp[0].Cnum);
        }
    break;

  case 10:
#line 98 "keyParse.yy"
    {
                yyval.Cnum = new complex<double>(*yyvsp[-1].Cnum);
                if (debug) {
                        cout << " ( " << *yyvsp[-1].Cnum << " ) = " << *yyval.Cnum << endl;
                }
                delete(yyvsp[-1].Cnum);
        }
    break;

  case 11:
#line 106 "keyParse.yy"
    {
                yyval.Cnum = new complex<double>(*yyvsp[-2].Cnum + *yyvsp[0].Cnum);
                if (debug) {
                        cout << *yyvsp[-2].Cnum << " + " << *yyvsp[0].Cnum << " = " << *yyval.Cnum << endl;
                }
                delete(yyvsp[-2].Cnum);
                delete(yyvsp[0].Cnum);
        }
    break;

  case 12:
#line 115 "keyParse.yy"
    {
                yyval.Cnum = new complex<double>(*yyvsp[-2].Cnum - *yyvsp[0].Cnum);
                if (debug) {
                        cout << *yyvsp[-2].Cnum << " - " << *yyvsp[0].Cnum << " = " << *yyval.Cnum << endl;
                }
                delete(yyvsp[-2].Cnum);
                delete(yyvsp[0].Cnum);
        }
    break;

  case 13:
#line 124 "keyParse.yy"
    {
                yyval.Cnum = new complex<double>(yyvsp[-2].Fnum * *yyvsp[0].Cnum);
                if (debug) {
                        cout << yyvsp[-2].Fnum << " * " << *yyvsp[0].Cnum << " = " << *yyval.Cnum << endl;
                }
                delete(yyvsp[0].Cnum);
        }
    break;

  case 14:
#line 132 "keyParse.yy"
    {
                yyval.Cnum = new complex<double>(*yyvsp[-2].Cnum * *yyvsp[0].Cnum);
                if (debug) {
                        cout << *yyvsp[-2].Cnum << " * " << *yyvsp[0].Cnum << " = " << *yyval.Cnum << endl;
                }
                delete(yyvsp[-2].Cnum);
                delete(yyvsp[0].Cnum);
        }
    break;

  case 15:
#line 142 "keyParse.yy"
    {
                wv.setDecay(*yyvsp[0].Decay);
                delete yyvsp[0].Decay;
                if (debug) {
                        cout << "@@Found a wave" << endl;
                        wv.print();
                        cout << "@@Filling wave" << endl;
                }
                wv.fill(e,debug);
                if (debug) {
                        cout << "@@Wave before boosts" << endl;
                        wv.print();
                }
                wv.setupFrames(debug);
                if (debug) {
                        cout << "@@Wave after boosts" << endl;
                        wv.print();
                }
                amp  = wv.decayAmp(debug);
                yyval.Cnum = new complex<double>(amp);
                nwave++;
        }
    break;

  case 16:
#line 164 "keyParse.yy"
    {
                wv.setDecay(*yyvsp[0].Decay);
                delete yyvsp[0].Decay;
                if (debug) {
                        cout << "@@Found a wave" << endl;
                        wv.print();
                        cout << "@@Filling wave" << endl;
                }
                wv.fill(e,debug);
                if (debug) {
                        cout << "@@Wave before boosts" << endl;
                        wv.print();
                }
                wv.setupFrames(debug);
                if (debug) {
                        cout << "@@Wave after boosts" << endl;
                        wv.print();
                }
                if (debug) {
                        cout << "This should compute decay amplitude expt wave" << endl;
                }
                double t = 0.0;
                fourVec t_init(0.0,threeVec(0.0,0.0,0.0));
                if (t_part_init == "beam") {
                        t_init = wv.getBeam();
                }
                else if (t_part_init == "target") {
                        t_init = wv.getTarget();
                }
                else {
                        cerr << "unknown initial t specifier: " << t_part_init << endl;
                        abort();
                }
                t = (t_init - *wv.get4P(t_part_final, debug)).lenSq();
                if (debug) {
                        cout << "calulating amplitude with t = " << t << endl;
                }
                delete t_part_final;

                wv.setT(t);
                amp  = wv.decayAmp(debug);
                yyval.Cnum = new complex<double>(amp);
                nwave++;
        }
    break;

  case 17:
#line 210 "keyParse.yy"
    {
            }
    break;

  case 18:
#line 214 "keyParse.yy"
    {
        }
    break;

  case 19:
#line 218 "keyParse.yy"
    {
                if(!strcmp(yyvsp[-2].string,"J")) wv.setJ(yyvsp[0].num);
                if(!strcmp(yyvsp[-2].string,"M")) wv.setM(yyvsp[0].num);
                if(!strcmp(yyvsp[-2].string,"P")) wv.setP(yyvsp[0].num);
        }
    break;

  case 20:
#line 225 "keyParse.yy"
    {
                decay* d = new decay;
                d->addChild(*yyvsp[-3].Particle);
                d->addChild(*yyvsp[-2].Particle);
                delete yyvsp[-3].Particle;
                delete yyvsp[-2].Particle;
                d->setL(yyvsp[-1].num);
                d->calculateS();
                yyval.Decay = d;
        }
    break;

  case 21:
#line 235 "keyParse.yy"
    {
                decay* d = new decay;
                d->addChild(*yyvsp[-5].Particle);
                d->addChild(*yyvsp[-4].Particle);
                delete yyvsp[-5].Particle;
                delete yyvsp[-4].Particle;
                if(!strcmp(yyvsp[-3].string,"l")) {
                        d->setL(yyvsp[-1].num);
                        d->calculateS();
                }
                else {
                        cerr << "unexpected field at line " << lineno << endl;
                        cerr << "found \'" << yyvsp[-3].string << "\'" << endl;
                        cerr << "expected \'l\'" << endl;
                        exit(1);
                }
                yyval.Decay = d;
        }
    break;

  case 22:
#line 253 "keyParse.yy"
    {
                decay* d = new decay;
                d->addChild(*yyvsp[-5].Particle);
                d->addChild(*yyvsp[-4].Particle);
                delete yyvsp[-5].Particle;
                delete yyvsp[-4].Particle;
                if(!strcmp(yyvsp[-3].string,"b")) {
                        wv.setSlope(yyvsp[-1].Fnum);
                }
                else {
                        cerr << "unexpected field at line " << lineno << endl;
                        cerr << "found \'" << yyvsp[-3].string << "\'" << endl;
                        cerr << "expected \'b\'" << endl;
                        exit(1);
                }
                yyval.Decay = d;
        }
    break;

  case 23:
#line 270 "keyParse.yy"
    {
                decay* d = new decay;
                d->addChild(*yyvsp[-8].Particle);
                d->addChild(*yyvsp[-7].Particle);
                delete yyvsp[-8].Particle;
                delete yyvsp[-7].Particle;
                if(!strcmp(yyvsp[-6].string,"l")) {
                        d->setL(yyvsp[-4].num);
                }
                else {
                        cerr << "expecting \'l\' at line " << lineno << endl;
                        cerr << "found \'" << yyvsp[-6].string << "\'" << endl;
                        exit(1);
                }
                if(!strcmp(yyvsp[-3].string,"s")) {
                        d->setS(yyvsp[-1].num);
                }
                else {
                        cerr << "expecting \'l\' at line " << lineno << endl;
                        cerr << "found \'" << yyvsp[-3].string << "\'" << endl;
                        exit(1);
                }
                yyval.Decay = d;
        }
    break;

  case 24:
#line 294 "keyParse.yy"
    {
                decay* d = new decay;
                d->addChild(*yyvsp[-4].Particle);
                d->addChild(*yyvsp[-3].Particle);
                delete yyvsp[-4].Particle;
                delete yyvsp[-3].Particle;
                d->setL(yyvsp[-2].num);
                d->setS(yyvsp[-1].num);
                yyval.Decay = d;
        }
    break;

  case 25:
#line 304 "keyParse.yy"
    {
                decay* d = new decay;
                d->addChild(*yyvsp[-4].Particle);
                d->addChild(*yyvsp[-3].Particle);
                d->addChild(*yyvsp[-2].Particle);
                delete yyvsp[-4].Particle;
                delete yyvsp[-3].Particle;
                delete yyvsp[-2].Particle;
                d->setL(yyvsp[-1].num);
                d->calculateS();
                yyval.Decay = d;
        }
    break;

  case 26:
#line 318 "keyParse.yy"
    {
                yyval.Particle = yyvsp[0].Particle;
        }
    break;

  case 27:
#line 321 "keyParse.yy"
    {
                yyvsp[-1].Particle->setDecay(*yyvsp[0].Decay);
                massDep* bw = new breitWigner();
                yyvsp[-1].Particle->setMassDep(bw);
                delete yyvsp[0].Decay;
                yyval.Particle = yyvsp[-1].Particle;
        }
    break;

  case 28:
#line 328 "keyParse.yy"
    {
                yyvsp[-4].Particle->setDecay(*yyvsp[-3].Decay);
                massDep* md;
                if (!strcmp(yyvsp[0].string,"bw")) {
                        md = new breitWigner();
                }
                else if (!strcmp(yyvsp[0].string,"amp")) {
                        md = new AMP_M();
                }
                else if (!strcmp(yyvsp[0].string,"amp_ves")) {
                        md = new AMP_ves();
                }
		else if (!strcmp(yyvsp[0].string,"amp_kach")) {
                        md = new AMP_kach();
                }
                else if (!strcmp(yyvsp[0].string,"flat")) {
                        md = new flat();
                }
                else {
                        cerr << "unknown mass dependence: " << yyvsp[0].string;
                        cerr << " at line " << lineno << endl;
                        exit(1);
                }
                yyvsp[-4].Particle->setMassDep(md);
                delete yyvsp[-3].Decay;
                yyval.Particle = yyvsp[-4].Particle;
        }
    break;

  case 29:
#line 354 "keyParse.yy"
    {
                yyval.Particle = yyvsp[0].Particle;
        }
    break;

  case 30:
#line 357 "keyParse.yy"
    {
                yyvsp[-3].Particle->addHelicity(yyvsp[0].num);
                yyval.Particle = yyvsp[-3].Particle;
        }
    break;

  case 31:
#line 363 "keyParse.yy"
    {
                yyval.Particle = yyvsp[0].Particle;
        }
    break;

  case 32:
#line 366 "keyParse.yy"
    {
                yyvsp[-3].Particle->Index(yyvsp[-1].num);
                yyval.Particle = yyvsp[-3].Particle;
        }
    break;

  case 33:
#line 372 "keyParse.yy"
    {
                particle* p = new particle(PDGtable.get(yyvsp[0].string),0);
                yyval.Particle = p;
        }
    break;

  case 34:
#line 376 "keyParse.yy"
    {
                particle* p = new particle(PDGtable.get(yyvsp[-1].string),+1);
                yyval.Particle = p;
        }
    break;

  case 35:
#line 380 "keyParse.yy"
    {
                particle* p = new particle(PDGtable.get(yyvsp[-1].string),-1);
                yyval.Particle = p;
        }
    break;

  case 36:
#line 384 "keyParse.yy"
    {
                particle* p = new particle(PDGtable.get(yyvsp[-1].string),0);
                yyval.Particle = p;
        }
    break;

  case 37:
#line 390 "keyParse.yy"
    {
                yyval.Cnum = new complex<double>(yyvsp[-3].Fnum,yyvsp[-1].Fnum);
        }
    break;


    }

/* Line 1010 of yacc.c.  */
#line 1544 "keyParse.cc"

  yyvsp -= yylen;
  yyssp -= yylen;


  YY_STACK_PRINT (yyss, yyssp);

  *++yyvsp = yyval;


  /* Now `shift' the result of the reduction.  Determine what state
     that goes to, based on the state we popped back to and the rule
     number reduced by.  */

  yyn = yyr1[yyn];

  yystate = yypgoto[yyn - YYNTOKENS] + *yyssp;
  if (0 <= yystate && yystate <= YYLAST && yycheck[yystate] == *yyssp)
    yystate = yytable[yystate];
  else
    yystate = yydefgoto[yyn - YYNTOKENS];

  goto yynewstate;


/*------------------------------------.
| yyerrlab -- here on detecting error |
`------------------------------------*/
yyerrlab:
  /* If not already recovering from an error, report this error.  */
  if (!yyerrstatus)
    {
      ++yynerrs;
#if YYERROR_VERBOSE
      yyn = yypact[yystate];

      if (YYPACT_NINF < yyn && yyn < YYLAST)
	{
	  YYSIZE_T yysize = 0;
	  int yytype = YYTRANSLATE (yychar);
	  const char* yyprefix;
	  char *yymsg;
	  int yyx;

	  /* Start YYX at -YYN if negative to avoid negative indexes in
	     YYCHECK.  */
	  int yyxbegin = yyn < 0 ? -yyn : 0;

	  /* Stay within bounds of both yycheck and yytname.  */
	  int yychecklim = YYLAST - yyn;
	  int yyxend = yychecklim < YYNTOKENS ? yychecklim : YYNTOKENS;
	  int yycount = 0;

	  yyprefix = ", expecting ";
	  for (yyx = yyxbegin; yyx < yyxend; ++yyx)
	    if (yycheck[yyx + yyn] == yyx && yyx != YYTERROR)
	      {
		yysize += yystrlen (yyprefix) + yystrlen (yytname [yyx]);
		yycount += 1;
		if (yycount == 5)
		  {
		    yysize = 0;
		    break;
		  }
	      }
	  yysize += (sizeof ("syntax error, unexpected ")
		     + yystrlen (yytname[yytype]));
	  yymsg = (char *) YYSTACK_ALLOC (yysize);
	  if (yymsg != 0)
	    {
	      char *yyp = yystpcpy (yymsg, "syntax error, unexpected ");
	      yyp = yystpcpy (yyp, yytname[yytype]);

	      if (yycount < 5)
		{
		  yyprefix = ", expecting ";
		  for (yyx = yyxbegin; yyx < yyxend; ++yyx)
		    if (yycheck[yyx + yyn] == yyx && yyx != YYTERROR)
		      {
			yyp = yystpcpy (yyp, yyprefix);
			yyp = yystpcpy (yyp, yytname[yyx]);
			yyprefix = " or ";
		      }
		}
	      yyerror (yymsg);
	      YYSTACK_FREE (yymsg);
	    }
	  else
	    yyerror ("syntax error; also virtual memory exhausted");
	}
      else
#endif /* YYERROR_VERBOSE */
	yyerror ("syntax error");
    }



  if (yyerrstatus == 3)
    {
      /* If just tried and failed to reuse lookahead token after an
	 error, discard it.  */

      if (yychar <= YYEOF)
        {
          /* If at end of input, pop the error token,
	     then the rest of the stack, then return failure.  */
	  if (yychar == YYEOF)
	     for (;;)
	       {
		 YYPOPSTACK;
		 if (yyssp == yyss)
		   YYABORT;
		 YYDSYMPRINTF ("Error: popping", yystos[*yyssp], yyvsp, yylsp);
		 yydestruct (yystos[*yyssp], yyvsp);
	       }
        }
      else
	{
	  YYDSYMPRINTF ("Error: discarding", yytoken, &yylval, &yylloc);
	  yydestruct (yytoken, &yylval);
	  yychar = YYEMPTY;

	}
    }

  /* Else will try to reuse lookahead token after shifting the error
     token.  */
  goto yyerrlab1;


/*---------------------------------------------------.
| yyerrorlab -- error raised explicitly by YYERROR.  |
`---------------------------------------------------*/
yyerrorlab:

#ifdef __GNUC__
  /* Pacify GCC when the user code never invokes YYERROR and the label
     yyerrorlab therefore never appears in user code.  */
  if (0)
     goto yyerrorlab;
#endif

  yyvsp -= yylen;
  yyssp -= yylen;
  yystate = *yyssp;
  goto yyerrlab1;


/*-------------------------------------------------------------.
| yyerrlab1 -- common code for both syntax error and YYERROR.  |
`-------------------------------------------------------------*/
yyerrlab1:
  yyerrstatus = 3;	/* Each real token shifted decrements this.  */

  for (;;)
    {
      yyn = yypact[yystate];
      if (yyn != YYPACT_NINF)
	{
	  yyn += YYTERROR;
	  if (0 <= yyn && yyn <= YYLAST && yycheck[yyn] == YYTERROR)
	    {
	      yyn = yytable[yyn];
	      if (0 < yyn)
		break;
	    }
	}

      /* Pop the current state because it cannot handle the error token.  */
      if (yyssp == yyss)
	YYABORT;

      YYDSYMPRINTF ("Error: popping", yystos[*yyssp], yyvsp, yylsp);
      yydestruct (yystos[yystate], yyvsp);
      YYPOPSTACK;
      yystate = *yyssp;
      YY_STACK_PRINT (yyss, yyssp);
    }

  if (yyn == YYFINAL)
    YYACCEPT;

  YYDPRINTF ((stderr, "Shifting error token, "));

  *++yyvsp = yylval;


  yystate = yyn;
  goto yynewstate;


/*-------------------------------------.
| yyacceptlab -- YYACCEPT comes here.  |
`-------------------------------------*/
yyacceptlab:
  yyresult = 0;
  goto yyreturn;

/*-----------------------------------.
| yyabortlab -- YYABORT comes here.  |
`-----------------------------------*/
yyabortlab:
  yyresult = 1;
  goto yyreturn;

#ifndef yyoverflow
/*----------------------------------------------.
| yyoverflowlab -- parser overflow comes here.  |
`----------------------------------------------*/
yyoverflowlab:
  yyerror ("parser stack overflow");
  yyresult = 2;
  /* Fall through.  */
#endif

yyreturn:
#ifndef yyoverflow
  if (yyss != yyssa)
    YYSTACK_FREE (yyss);
#endif
  return yyresult;
}


#line 394 "keyParse.yy"



