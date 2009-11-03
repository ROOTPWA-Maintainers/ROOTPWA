%option prefix="key"
%option noyywrap
%option nounput
%{

#include <iostream>
#include <particle.h>
#include <keyParse.h>

using namespace std;

#undef YY_INPUT
#define YY_INPUT(buf,result,max_size)\
{\
        int c = getc(yyin);\
        result = (c==EOF)?YY_NULL:(buf[0]=c,1);\
}

extern particleDataTable PDGtable;
particleData p;
int i;
int lineno=1;
char* fname;

%}

string [a-zA-Z][a-zA-Z0-9()_']*
float   [+-]?([0-9]+\.?)([Ee][+-][0-9]+)?|[+-]?([0-9]*\.[0-9]+)([Ee][+-][0-9]+)?
comment #.*$

%%

[ \t]+  {
                ;
        }
        
{comment}       {
                ;
        }
        
[-+]?[0-9][0-9]*        {
                keylval.num = atoi(keytext);
                return (INT);
        }

{float} {
                keylval.Fnum = atof(keytext);
                return (FLOAT);
        }

debug   {
                return(DEBUG);
        }

channel {
                return(CHANNEL);
        }

mode    {
                return(MODE);
        }

h|(helicity)    {
                return(HELICITY);
        }

(md)|(massdep)  {
                return(MASSDEP);
        }

{string}        {
                sscanf(keytext,"%s",keylval.string);
                p = PDGtable.get(keytext);
                if (p.Name() != "") {
                        return (PARTICLENAME);
                }
                else {
                        return (STRING);
                }
        }

\n      {
                lineno++;
        }

.       {
                return ((int) keytext[0]);
        }

<<EOF>>        {
                yyterminate();
        }
%%


void keyerror(char* s) {
        cerr << fname << ":" << lineno << " " << s << " at " << keytext << endl;
}

