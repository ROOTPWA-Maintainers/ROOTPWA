#!/bin/bash
# scrip doing things on all mass bins
# Iospin symmetrizing

for i in $1*;
do
  echo " ----- Mass Bin $i";
  #cd $i/PSPAMPS;

  #addamp $ROOTPWA/keyfiles/key5pi/f1list.dat OLD/
  #addamp $ROOTPWA/keyfiles/key5pi/f1plist.dat OLD/
  #addamp $ROOTPWA/keyfiles/key5pi/f21565list.dat OLD/
  #addamp $ROOTPWA/keyfiles/key5pi/f0list.dat OLD/
  #addamp $ROOTPWA/keyfiles/key5pi/b1list.dat OLD/
  #addamp $ROOTPWA/keyfiles/key5pi/etalist.dat OLD/
  #addamp $ROOTPWA/keyfiles/key5pi/eta2list.dat OLD/

  #addamp $ROOTPWA/keyfiles/key5pi/f2list.dat OLD/
  #addamp $ROOTPWA/keyfiles/key5pi/rhoprimelist.dat OLD/


  #int *.amp > norm.int;
  #cp norm.int $i/AMPS;

  cd $i/ACCAMPS;

  addamp $ROOTPWA/keyfiles/key5pi/f1list.dat OLD/
  addamp $ROOTPWA/keyfiles/key5pi/f1plist.dat OLD/
  addamp $ROOTPWA/keyfiles/key5pi/f21565list.dat OLD/
  addamp $ROOTPWA/keyfiles/key5pi/f0list.dat OLD/
  addamp $ROOTPWA/keyfiles/key5pi/b1list.dat OLD/
  addamp $ROOTPWA/keyfiles/key5pi/etalist.dat OLD/
  addamp $ROOTPWA/keyfiles/key5pi/eta2list.dat OLD/

  addamp $ROOTPWA/keyfiles/key5pi/f2list.dat OLD/
  addamp $ROOTPWA/keyfiles/key5pi/rhoprimelist.dat OLD/


  int *.amp > accnorm.int;
  cp accnorm.int $i/AMPS;


 cd $i/AMPS;

  addamp $ROOTPWA/keyfiles/key5pi/f1list.dat OLD/
  addamp $ROOTPWA/keyfiles/key5pi/f1plist.dat OLD/
  addamp $ROOTPWA/keyfiles/key5pi/f21565list.dat OLD/
  addamp $ROOTPWA/keyfiles/key5pi/f0list.dat OLD/
  addamp $ROOTPWA/keyfiles/key5pi/b1list.dat OLD/
  addamp $ROOTPWA/keyfiles/key5pi/etalist.dat OLD/
  addamp $ROOTPWA/keyfiles/key5pi/eta2list.dat OLD/

  addamp $ROOTPWA/keyfiles/key5pi/f2list.dat OLD/
  addamp $ROOTPWA/keyfiles/key5pi/rhoprimelist.dat OLD/

  echo "------------------"

done;
