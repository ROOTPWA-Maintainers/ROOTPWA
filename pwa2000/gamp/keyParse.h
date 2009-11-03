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
/* Line 1285 of yacc.c.  */
#line 72 "y.tab.h"
# define yystype YYSTYPE /* obsolescent; will be withdrawn */
# define YYSTYPE_IS_DECLARED 1
# define YYSTYPE_IS_TRIVIAL 1
#endif

extern YYSTYPE keylval;



