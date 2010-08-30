################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CPP_SRCS += \
../src/rootpwa_gui.cpp 

CC_SRCS += \
../src/TrpwaMainFrame.cc \
../src/TrpwaSessionManager.cc \
../src/TrpwaWaveSelectFrame.cc \
../src/rpwaDict.cc 

OBJS += \
./src/TrpwaMainFrame.o \
./src/TrpwaSessionManager.o \
./src/TrpwaWaveSelectFrame.o \
./src/rootpwa_gui.o \
./src/rpwaDict.o 

CC_DEPS += \
./src/TrpwaMainFrame.d \
./src/TrpwaSessionManager.d \
./src/TrpwaWaveSelectFrame.d \
./src/rpwaDict.d 

CPP_DEPS += \
./src/rootpwa_gui.d 


# Each subdirectory must supply rules for building sources it contributes
src/%.o: ../src/%.cc
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C++ Compiler'
	g++ -I"/home/Promme/bin/root/include" -I/home/Promme/bin/libconfig/libconfig-1.4.1/lib -O0 -g3 -Wall -c -fmessage-length=0 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@:%.o=%.d)" -o"$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '

src/%.o: ../src/%.cpp
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C++ Compiler'
	g++ -I"/home/Promme/bin/root/include" -I/home/Promme/bin/libconfig/libconfig-1.4.1/lib -O0 -g3 -Wall -c -fmessage-length=0 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@:%.o=%.d)" -o"$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


