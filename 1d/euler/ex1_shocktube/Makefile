# Makefile
COMMON_LIB = ../../../lib
LIB = ../../lib
APPLIB = ../lib

CXX = g++

LINK = time $(CXX) \
-O3 \
-I$(COMMON_LIB) \
-I$(LIB) \
-I$(APPLIB) \
-std=c++14 \
-Wall \
-Wno-unknown-pragmas \
-Wpedantic \
-Wno-variadic-macros \
-Wno-c++11-extra-semi \
-Wno-gnu-zero-variadic-macro-arguments \
-Wno-c++11-long-long \
-DCHECK_BOUNDS

COMPILE = $(LINK) -c $(INCLUDE_LIBS)

OBJECTS = \
main.o \
initial_condition.o \
$(APPLIB)/AppSolver.o \
$(COMMON_LIB)/DogArrays.o \
$(COMMON_LIB)/setup.o \
$(COMMON_LIB)/assert.o \
$(COMMON_LIB)/debug.o \
$(COMMON_LIB)/dog_str.o \
$(COMMON_LIB)/dog_math.o \
$(COMMON_LIB)/IniDocument.o \
$(LIB)/wave_prop_method_1D.o \
$(LIB)/Params1D.o \
$(LIB)/HelperArrays1D.o \
$(LIB)/waveprop1D.o \
$(LIB)/single_time_step.o \

SOURCES = \
main.cpp \
initial_condition.cpp \
$(APPLIB)/AppSolver.cpp \
$(COMMON_LIB)/DogArrays.cpp \
$(COMMON_LIB)/setup.cpp \
$(COMMON_LIB)/assert.cpp \
$(COMMON_LIB)/debug.cpp \
$(COMMON_LIB)/dog_str.cpp \
$(COMMON_LIB)/dog_math.cpp \
$(COMMON_LIB)/IniDocument.cpp \
$(LIB)/wave_prop_method_1D.cpp \
$(LIB)/Params1D.cpp \
$(LIB)/HelperArrays1D.cpp \
$(LIB)/waveprop1D.cpp \
$(LIB)/single_time_step.cpp \

main.exe: $(OBJECTS) $(SOURCES)
	$(LINK) $(OBJECTS) -o main.exe

.cpp.o: ; $(COMPILE) $*.cpp -o $*.o

clean:
	rm -f $(OBJECTS) *.exe *.pdf; rm -f -r output; rm -f -r output_combined; rm -f -r __pycache__;
