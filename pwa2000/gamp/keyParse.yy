%{
#include <stream.h>
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

%}

%union{
        int num;
        std::complex<double>* Cnum;
        double Fnum;
        char string[100];
        decay* Decay;
        particle* Particle;
}

%token <num> INT
%token <Fnum> FLOAT
%token <Cnum> COMPLEX
%token <string> SIGN
%token <string> STRING
%token <string> PARTICLENAME
%token <string> DEBUG
%token <string> CHANNEL
%token <string> MODE
%token <string> MASSDEP
%token <string> HELICITY

%type <Decay> decay
%type <Particle> pstate
%type <Particle> particle
%type <Particle> particleID
%type <Particle> particleCharge
%type <Cnum> wave
%type <Cnum> waveexp
%type <Cnum> complex
%left '+' '-'
%left '*'

%%
input:    /* empty */
        | input statement 
        ;

statement:   DEBUG '=' INT ';' {
                debug = $3;
        }
        | CHANNEL '=' STRING ';'{
                wv.channel($3);
        }
        | MODE  '=' STRING ';' {
                mode = $3;
        }
        | STRING  '=' STRING ',' particleID ';' {
                t_part_init = $3;
                t_part_final = $5;
        }
        | waveexp ';' {
                if (mode == "binary") {
                        cout.write((char*) $1,sizeof(complex<double>));
                }
                else {
                        cout << "Mass = " << ~(wv.get4P()) << "\t";
                        if ( wv.channel() == "t" ) {
                                cout << "t = " << (e.beam().get4P()-wv.get4P()).lenSq() << "\t";
                        }
                        cout << "Amp = " <<  *$1 << endl;
                }
		$$ = new complex<double>(*$1);
                delete($1);
        }
        ;

waveexp: wave {
                $$ = new complex<double>(*$1);
                delete($1);
        }
        | '(' waveexp ')' {
                $$ = new complex<double>(*$2);
                if (debug) {
                        cout << " ( " << *$2 << " ) = " << *$$ << endl;
                }
                delete($2);
        }
        |
        waveexp '+' waveexp {
                $$ = new complex<double>(*$1 + *$3);
                if (debug) {
                        cout << *$1 << " + " << *$3 << " = " << *$$ << endl;
                }
                delete($1);
                delete($3);
        }
        |
        waveexp '-' waveexp {
                $$ = new complex<double>(*$1 - *$3);
                if (debug) {
                        cout << *$1 << " - " << *$3 << " = " << *$$ << endl;
                }
                delete($1);
                delete($3);
        }
        |
        FLOAT '*' waveexp {
                $$ = new complex<double>($1 * *$3);
                if (debug) {
                        cout << $1 << " * " << *$3 << " = " << *$$ << endl;
                }
                delete($3);
        }
        |
        complex '*' waveexp {
                $$ = new complex<double>(*$1 * *$3);
                if (debug) {
                        cout << *$1 << " * " << *$3 << " = " << *$$ << endl;
                }
                delete($1);
                delete($3);
        }
        ;

wave:   resonance decay {
                wv.setDecay(*$2);
                delete $2;
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
                $$ = new complex<double>(amp);
                nwave++;
        }
        | decay {
                wv.setDecay(*$1);
                delete $1;
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
                $$ = new complex<double>(amp);
                nwave++;
        }
        ;

resonance:  quantum_nums {
            }
            ;

quantum_nums: quantum_num quantum_num quantum_num {
        }
        ;

quantum_num: STRING '=' INT {
                if(!strcmp($1,"J")) wv.setJ($3);
                if(!strcmp($1,"M")) wv.setM($3);
                if(!strcmp($1,"P")) wv.setP($3);
        }
        ;

decay:      '{' particle particle INT '}'  {
                decay* d = new decay;
                d->addChild(*$2);
                d->addChild(*$3);
                delete $2;
                delete $3;
                d->setL($4);
                d->calculateS();
                $$ = d;
        }
        | '{' particle particle STRING '=' INT '}'  {
                decay* d = new decay;
                d->addChild(*$2);
                d->addChild(*$3);
                delete $2;
                delete $3;
                if(!strcmp($4,"l")) {
                        d->setL($6);
                        d->calculateS();
                }
                else {
                        cerr << "unexpected field at line " << lineno << endl;
                        cerr << "found \'" << $4 << "\'" << endl;
                        cerr << "expected \'l\'" << endl;
                        exit(1);
                }
                $$ = d;
        }
        | '{' particle particle STRING '=' FLOAT '}'  {
                decay* d = new decay;
                d->addChild(*$2);
                d->addChild(*$3);
                delete $2;
                delete $3;
                if(!strcmp($4,"b")) {
                        wv.setSlope($6);
                }
                else {
                        cerr << "unexpected field at line " << lineno << endl;
                        cerr << "found \'" << $4 << "\'" << endl;
                        cerr << "expected \'b\'" << endl;
                        exit(1);
                }
                $$ = d;
        }
        | '{' particle particle STRING '=' INT STRING '=' INT '}'  {
                decay* d = new decay;
                d->addChild(*$2);
                d->addChild(*$3);
                delete $2;
                delete $3;
                if(!strcmp($4,"l")) {
                        d->setL($6);
                }
                else {
                        cerr << "expecting \'l\' at line " << lineno << endl;
                        cerr << "found \'" << $4 << "\'" << endl;
                        exit(1);
                }
                if(!strcmp($7,"s")) {
                        d->setS($9);
                }
                else {
                        cerr << "expecting \'l\' at line " << lineno << endl;
                        cerr << "found \'" << $7 << "\'" << endl;
                        exit(1);
                }
                $$ = d;
        }
        | '{' particle particle INT INT '}'  {
                decay* d = new decay;
                d->addChild(*$2);
                d->addChild(*$3);
                delete $2;
                delete $3;
                d->setL($4);
                d->setS($5);
                $$ = d;
        }
        | '{' particle particle particle INT '}'  {
                decay* d = new decay;
                d->addChild(*$2);
                d->addChild(*$3);
                d->addChild(*$4);
                delete $2;
                delete $3;
                delete $4;
                d->setL($5);
                d->calculateS();
                $$ = d;
        }
        ;

particle:   pstate {
                $$ = $1;
        }
        | pstate decay {
                $1->setDecay(*$2);
                massDep* bw = new breitWigner();
                $1->setMassDep(bw);
                delete $2;
                $$ = $1;
        }
        | pstate decay MASSDEP '=' STRING{
                $1->setDecay(*$2);
                massDep* md;
                if (!strcmp($5,"bw")) {
                        md = new breitWigner();
                }
                else if (!strcmp($5,"amp")) {
                        md = new AMP_M();
                }
                else if (!strcmp($5,"amp_ves")) {
                        md = new AMP_ves();
                }
                else if (!strcmp($5,"flat")) {
                        md = new flat();
                }
                else {
                        cerr << "unknown mass dependence: " << $5;
                        cerr << " at line " << lineno << endl;
                        exit(1);
                }
                $1->setMassDep(md);
                delete $2;
                $$ = $1;
        }
        ;

pstate:  particleID {
                $$ = $1;
        }
        | particleID HELICITY '=' INT {
                $1->addHelicity($4);
                $$ = $1;
        }
        ;
        
particleID:  particleCharge {
                $$ = $1;
        }
        | particleCharge '[' INT ']' {
                $1->Index($3);
                $$ = $1;
        }
        ;
        
particleCharge: PARTICLENAME {
                particle* p = new particle(PDGtable.get($1),0);
                $$ = p;
        }
        | PARTICLENAME '+' {
                particle* p = new particle(PDGtable.get($1),+1);
                $$ = p;
        }
        | PARTICLENAME '-' {
                particle* p = new particle(PDGtable.get($1),-1);
                $$ = p;
        }
        | PARTICLENAME '0' {
                particle* p = new particle(PDGtable.get($1),0);
                $$ = p;
        }
        ;

complex:        '(' FLOAT ',' FLOAT ')' {
                $$ = new complex<double>($2,$4);
        }

%%

