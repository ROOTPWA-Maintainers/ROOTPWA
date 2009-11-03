#line 81 "../int.nw"
#include <iostream>
#include <string>
#include <list>
#include <cstdlib>
#include <unistd.h>
#include <integral.h>

#line 158 "../int.nw"
void printUsage(char* prog) {
    cerr << "usage:"  << endl;
    cerr << "  creating     : " << prog \
         << " [-d] [-m max] [-r n] files" << endl;
    cerr << "  adding       : " << prog \
         << " [-d] [-m max] [-r n] -a intfile" << endl;
    cerr << "  display evts : " << prog \
         << " [-q] -i intfile" << endl;
    cerr << "  renormalizing: " << prog \
         << " [-r n] -i intfile" << endl;
    cerr << "  help         : " << prog \
         << " [-h]" << endl;
    cerr << "where:" << endl;
    cerr << "    max    : maximum number of events" << endl;
    cerr << "    n      : number of events to renormalize to" << endl;
    cerr << "    intfile: integral file to read" <<  endl;
    cerr << "             (for adding to or renormalizing)" << endl;
    cerr << endl;
    exit(0);
}



#line 90 "../int.nw"
int main(int argc, char** argv) {

    int debug = 0;
    int maxevents = 0;
    int renorm = 0;
    int add = 1;
    int display_events = 0;
    string oldint;
    integral ni;

    
#line 189 "../int.nw"
    extern int optind;
    extern char* optarg;
    int c;
    if (argc == 1) printUsage(argv[0]);
    while ( (c = getopt(argc,argv, "dm:r:a:i:h:q")) != -1 )
        switch(c) {
            case 'd':
                debug = 1;    
                break;
            case 'q':
                display_events = 1;
                break;
            case 'm':
                maxevents = atoi(optarg);    
                break;
            case 'r':
                renorm = atoi(optarg);    
                break;
            case 'a':
                oldint = optarg;    
                break;
            case 'i':
                add = 0;
                oldint = optarg;    
                break;
            case 'h':
            case '?':
                printUsage(argv[0]);
        }

#line 102 "../int.nw"
    
#line 127 "../int.nw"
    if (oldint.size()) {
        ifstream oldfile(oldint.c_str());
        ni.scan(oldfile);
    }
    else {
        ni.files(argv+optind);
    }


#line 104 "../int.nw"
    if (display_events) {
       ni.print_events();
    } else {

       if (add) {
           
#line 143 "../int.nw"
    ni.max(maxevents);
    try {
        ni.integrate();
    }
    catch (char* m) {
        cerr << m << endl;
        return 0;
    }

#line 110 "../int.nw"
       }
       if (renorm) {
          ni.renormalize(renorm);
       }
       ni.print();
    }
    return 0;
}


