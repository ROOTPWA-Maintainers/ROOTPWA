TFracNum.o:	TFracNum.cc TFracNum.h	
	g++  -c -g TFracNum.cc -o TFracNum.o
ClebschGordanBox.o:	ClebschGordanBox.cc ClebschGordanBox.h TFracNum.o
	g++  -c ClebschGordanBox.cc -o ClebschGordanBox.o
TJwfTensor.o:	TJwfTensor.cc TJwfTensor.h ClebschGordanBox.o
	g++  -c TJwfTensor.cc -o TJwfTensor.o
TSpinWaveFunction.o: TSpinWaveFunction.cc TSpinWaveFunction.h TJwfTensor.o
	g++  -c TSpinWaveFunction.cc -o TSpinWaveFunction.o
TLSAmpl.o:	TLSAmpl.cc TLSAmpl.h TSpinWaveFunction.o
	g++  -c TLSAmpl.cc -o TLSAmpl.o
TFhh.o:	TFhh.cc TFhh.h TLSAmpl.o
	g++  -c TFhh.cc -o TFhh.o
CalcAmpl++: CalcAmpl.cc TFhh.o 
	g++  CalcAmpl.cc -o CalcAmpl++ TFracNum.o ClebschGordanBox.o TJwfTensor.o TSpinWaveFunction.o TLSAmpl.o TFhh.o 