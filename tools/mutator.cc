///////////////////////////////////////////////////////////////////////////
//
//    Copyright 2009 Sebastian Neubert
//
//    This file is part of rootpwa
//
//    rootpwa is free software: you can redistribute it and/or modify
//    it under the terms of the GNU General Public License as published by
//    the Free Software Foundation, either version 3 of the License, or
//    (at your option) any later version.
//
//    rootpwa is distributed in the hope that it will be useful,
//    but WITHOUT ANY WARRANTY; without even the implied warranty of
//    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//    GNU General Public License for more details.
//
//    You should have received a copy of the GNU General Public License
//    along with rootpwa.  If not, see <http://www.gnu.org/licenses/>.
//
///////////////////////////////////////////////////////////////////////////

#include <time.h>
#include "TH1I.h"
#include <algorithm>
#include <assert.h>
#include "wset.h"


void usage(const char* progname){
  cerr << "usage: -P <wave pool> -S <seed>  <ancestorlist>" << endl;
  cerr << "[-F <number of fixed waves (first N waves in list)> default: 1 ]" << endl;
  cerr << "[-E <number of waves to exchange> default: 5 ]" << endl;
  cerr << "[-D <number of waves to drop> ]" << endl;
  cerr << "[-A <number of waves to add> ]" << endl;
  cerr << "[-V <range of number of waves to vay> ]" << endl;
  cerr << "[-p 1 < selective pressure < 2,  default: 1.5 ]" << endl;
  cerr << "[-c 0 < crossover probability < 1 default> 0.75 ]" << endl;
  cerr << progname << endl;
}


///////////////////////////////////////////////////////
//
//         Main starts here
//
///////////////////////////////////////////////////////


int
main(int argc, char** argv){

  char* progname=argv[0];

  vector<TString> inputlistFiles; // ancestor list
  vector<double> fitness;
  TString wavepoolFile;  // pool of waves to choose from
  unsigned int Exch=0;      // number of waves to exchange
  unsigned int Drop=0;      // number of waves to drop
  unsigned int Add=0;
  unsigned int Fix=1;
  unsigned int Vary=0;      // range of waves to add or subtract
  unsigned int seed=0;
  
  double sel=1.5; // selective pressure 1<c<2
  double pc=0.75;  // crossover probability;

extern char *optarg;
  // extern int optind;
  int c;
  unsigned int optcount=1;
  while ((c = getopt(argc, argv, "hE:A:V:D:P:F:S:c:p:")) != -1){
    optcount+=1;
    switch (c) {
    case 'F':
      Fix= atoi(optarg);
      //cerr << "Fix=" << Fix << endl;
      break;
    case 'E':
      Exch= atoi(optarg);
      //cerr << "Exchange=" << Exch << endl;
      break;
    case 'D':
      Drop= atoi(optarg);
      //cerr << "Drop=" << Drop << endl;
      break;
    case 'A':
      Add= atoi(optarg);
      //cerr << "Add=" << Add << endl;
      break;
    case 'V':
      Vary= atoi(optarg);
      //cerr << "Vary=" << Add << endl;
      break;
    case 'p':
      sel = atof(optarg);
      //cerr << "Vary=" << Add << endl;
      break;
    case 'c':
      pc = atof(optarg);
      //cerr << "Vary=" << Add << endl;
      break;
    case 'P':
      wavepoolFile=optarg;
      //cerr << "Pool=" << wavepoolFile << endl;
      break;
    case 'S':
      seed= atoi(optarg);
      //cerr << "Seed=" << seed << endl;
      break;
    case 'h':
      usage(progname);
      return 1;
    }
  }

  //cerr << "Optcount=" << optcount << endl;

  ////////////////////////////////////////////////////////////////
  ///              Read input lists
  ////////////////////////////////////////////////////////////////


  TString ancestorList(argv[optcount]);
  ifstream ancfile(ancestorList.Data());
  while(!ancfile.eof()){
    TString ancestor;
    double fitn=0;
    ancfile >> ancestor >> fitn;
    ancfile.ignore(200,'\n');
    if(ancestor.Length()<1 || fitn==0)continue;
    //cerr << ancestor << " : " << fitn << endl;
    inputlistFiles.push_back(ancestor);
    fitness.push_back(fitn);
  }
  

  if(inputlistFiles.size()<1){
    cerr << "No ancestor wavelists given. Exiting." << endl;
    return 1;
  }
  else {
    cerr << inputlistFiles.size() << " ancestors." << endl;
  }

  TString inputlistFile=inputlistFiles[0];
  

  // initialize Random number generator
  cerr << "Seed =" << seed << endl;
  gRandom->SetSeed(time(NULL)+seed);

  
  // open and read input wavelists (taken from TPWALikelyhood)
  vector<set<wsetentry>*> ancestorlist;
  unsigned int na=inputlistFiles.size();
  unsigned int mingenes=1000;
  unsigned int mini=0;
  for(unsigned int ia=0;ia<na;++ia){
    ancestorlist.push_back(new set<wsetentry>());
    readWavelist(*(ancestorlist[ia]),inputlistFiles[ia]);
    unsigned int n=ancestorlist[ia]->size();
    //cerr << "Ancestor: "<<inputlistFiles[ia]
    //	 << " with "<<n<<" genes." << endl;
    if(n<mingenes){
      mingenes=n;
      mini=ia;
    }
    if(n<Fix){
      cerr << "Ancestor has not enough genes to fix "
	   <<Fix<<". Fixing all" << endl;
      Fix=n;
    }
  }
  
  ////////////////////////////////////////////////////////////////
  ///              Ranked Selection 
  ////////////////////////////////////////////////////////////////


  // build selection histogram (probability distribution of ranked chromosomes)
  assert(na>=2);
  TH1I hDist("hDist","hDist",na,0,na-1);

  // we implement a ranked selection here such that:
  // Pmax= c/na, Pmin=(2-c)/na,  c is selective pressure 1<c<2

  for(unsigned int ia=0;ia<na;++ia){
    // reverse ranking:
    double rank=(double)(na-ia)-1.;
    double npop=(double)na;
    double npopa=npop-1;
    double p=1./npop*(2.-sel+2.*(sel-1)*rank/npopa);
    //cerr << "Ancestor " << ia << "  P="<<p<< endl;
    hDist.Fill(ia,p*npop);
  }

  // Choose two ancestors
  unsigned int x1=(unsigned int)hDist.GetRandom();
  unsigned int x2=(unsigned int)hDist.GetRandom();
  
  assert(x1>=0 && x1<na);
  assert(x2>=0 && x2<na);

  cerr << "Selective Pressure: " << sel << endl;
  cerr << "Choosing ancestors: " << x1 << " and " << x2 << endl;

  set<wsetentry>* ancestor1=ancestorlist[x1];
  set<wsetentry>* ancestor2=ancestorlist[x2];




  ////////////////////////////////////////////////////////////////
  ///              Cross Over
  ////////////////////////////////////////////////////////////////
  
  cerr << "Crossover probability: "<< pc << endl;

  // copy ancestor1
  set<wsetentry> motherlist;
  copy(ancestor1->begin(),ancestor1->end(),inserter(motherlist,motherlist.begin()));

  if(gRandom->Uniform()<pc){

    cerr << " crossing over ... ";
    // get two points in the genome
    set<wsetentry>* smallancestor=NULL;
    set<wsetentry>* largeancestor=NULL;

    if(ancestor1->size()<=ancestor2->size()){
     smallancestor = ancestor1;
     largeancestor = ancestor2;
    }
    else {
      smallancestor = ancestor2;
      largeancestor = ancestor1;
    }



    unsigned int cp1=gRandom->Integer(smallancestor->size()-Fix)+Fix;
    unsigned int cp2=gRandom->Integer(smallancestor->size()-Fix)+Fix;

    unsigned int cpp1 = cp1; // first point which is swapped
    unsigned int cpp2 = cp2; // second point ... 
    if(cp1>cp2){
      cpp1 = cp2;
      cpp2 = cp1;
    }
    
    cerr << " [" << cpp1 << "..." << cpp2 << "]" << endl;
    // move to cpp1
    set<wsetentry>::iterator ita1=motherlist.begin();
    set<wsetentry>::iterator ita2=ancestor2->begin();
    for(unsigned int ip=0; ip<cpp1; ++ip){
      ++ita1;++ita2;
    }
    // now see if we can exchange genes 
    // (we only exchange if the mother would receive a new gene!
    for(unsigned int ip=cpp1; ip<cpp2; ++ip){
      wsetentry gene=*ita2;
      if(findByName(motherlist,gene)==motherlist.end()){
	cerr << "crossing gene " << ip << ": "
	     << gene.name << endl;
	motherlist.erase(*ita1);
	motherlist.insert(gene);
	ita1=find(motherlist.begin(),motherlist.end(),gene);
      }
      ++ita1;++ita2;
    }
  }
  
  ////////////////////////////////////////////////////////////////
  ///              Mutation
  ////////////////////////////////////////////////////////////////

  set<wsetentry> wavepool;
  readWavelist(wavepool,wavepoolFile);

  // calculate dimension:
  unsigned int nwaves=motherlist.size();
  
  
  set<wsetentry> redpool; // reduced_wavepool=wavepool\motherlist
  //for(unsigned int i=0; i<nwaves;++i)cout<<motherlist[i].first<<endl;
  set<wsetentry>::iterator it1=wavepool.begin();
  while(it1!=wavepool.end()){
    if(find(motherlist.begin(),motherlist.end(),*it1)==motherlist.end()){
      redpool.insert(*it1);
    }
    ++it1;
  }

  // we have taken out all waves from the mother list
  unsigned int npool=redpool.size();
  cerr << wavepool.size() << " waves in pool." << endl;
  cerr << nwaves << " waves in ancestor." << endl;
  cerr << npool << " waves in reduced pool." << endl;
  
  set<wsetentry> mutantlist(motherlist);

  cerr << "Add =" << Add << endl;
  cerr << "Drop=" << Drop << endl;
  cerr << "Exch=" << Exch << endl;

  // draw vary value
  if(Vary>0){
    int v= gRandom->Integer(2*Vary)-Vary;
    cerr << "Var ="<<v<<  endl;
    if(v>0)Add+=v;
    if(v<0){
      // frist reduce added
      if((int)Add<-v){
	Drop+=-Add-v;
	Add=0;
      }
      else Add+=v;
    }
  }

  cerr << "Add =" << Add << endl;
  cerr << "Drop=" << Drop << endl;

  unsigned int exchangable=mutantlist.size()-Fix;
  if(Exch+Drop>exchangable){
    Drop=exchangable-Exch;
    cerr << "Adjusting drop to maximal value "<<Drop<<endl;
  }

  // Step1: Remove waves. Prepare for exchange
  for(unsigned int ie=0;ie<Exch+Drop;++ie){
    // choose one wave to remove
    mutantlist.erase(randomEntry(mutantlist,Fix));
  }

  cerr << "removed " << Exch+Drop << " waves" << endl;

  // Setp2: Exchange and add
  for(unsigned int ie=0;ie<Exch+Add;++ie){
    set<wsetentry>::iterator it=randomEntry(redpool);
    wsetentry anentry=*it;
    anentry.index=nwaves+ie;
    mutantlist.insert(anentry);
    redpool.erase(it); // remove from pool to avoid double waves
  }

  cerr << "added " << Exch+Add << " waves" << endl;

  // output
  set<wsetentry>::iterator it=mutantlist.begin();
  int count=0;
  while(it!=mutantlist.end()){
    // make sure that there are no double entries
    set<wsetentry>::iterator it2=it;
    ++it2;
    if(find(it2,mutantlist.end(),*it)!=mutantlist.end()){
      cerr << "Double entry " << it->name << endl;
      ++it;
      continue;
    }
    cout << it->name;
    if((it->threshold)!=0) cout << " " << it->threshold;
    cout<<endl;
    ++count;
    ++it;
  }
  cerr << count << " waves in mutant." << endl;
  return 0;
}
