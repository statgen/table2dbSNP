# Name of your executable, replace sampleProgram with your program name
EXE=table2dbSNP
########################
# Specify your Files:
# TOOLBASE - files with both .cpp and .h.  Just include the base name.
# Example if you have Class1.cpp, Class1.h, Class2.cpp, and Class2.cpp, specify:
# TOOLBASE = Class1 Class2
TOOLBASE = 
# SRCONLY - .cpp only files.  Specify the entire name, including the .cpp.
# Example, if you have Main.cpp, and File1.cpp, specify:
# SRCONLY = Main.cpp File1.cpp
SRCONLY = Main.cpp
# HDRONLY - .h only files.  Specify the entire name, including the .h.
# Example, if you have Template1.h and Template2.h, specify:
# HDRONLY = Template1.h Template2.h
HDRONLY = 

########################
# Optional settings if you need them.

# Specify the version you want for your program. 
VERSION?=0.1.1

# Specify any additional compile varaibles you want
# Example:
# DATE=$(shell date)
# USER=$(shell whoami)
# USER_COMPILE_VARS = -DDATE="\"${DATE}\"" -DUSER="\"${USER}\"" -lcrypto
USER_COMPILE_VARS = 

# Specify any additional includes you need.
USER_INCLUDES =

# Specify any additional warnings you need to enable.
# Example: 
# USER_WARNINGS ?= -Werror -Wno-unused-variable $(shell if [ X$(CCVERSION) \> X4.2.0 ] ; then echo " -Wno-unused-result" ; fi)
USER_WARNINGS =

# If you don't want to include the current directory (-I.), uncomment this line:
#CURRENT_DIR_INCLUDE=

# Specify any TOOLBASE that needs to be re-compiled when any other file
# is modified.  Do NOT specify this same file in TOOLBASE.
# Example: if you have Executable.cpp and Executable.h, specify:
# COMPILE_ANY_CHANGE = BamExecutable
COMPILE_ANY_CHANGE = 

# If you want to overwrite the default directory for the binaries, from the
# Default setting of 'bin', uncomment the next line and specify a value.
# BIN_DIR = bin

########################
# Library Path:
#
# Default path to the statgen library.  You can either update your copy of
# this Makefile to be where you want to get the Makefile from or you can 
# overwrite LIB_PATH_GENERAL or LIB_PATH_SAMPLE_PROGRAM.
# This design is such that if you have multiple programs, you could just set
# LIB_PATH_GENERAL to the location they should all use to find the library.
# If one program needs a different path, you would set LIB_PATH_SAMPLE_PROGRAM
# (replacing SAMPLE_PROGRAM with the appropriate name) to that different path.
LIB_PATH_GENERAL ?=../libStatGen
LIB_PATH_SAMPLE_PROGRAM ?= $(LIB_PATH_GENERAL)

########################
# Include the base Makefile
include $(LIB_PATH_SAMPLE_PROGRAM)/Makefiles/Makefile.src

########################
# Add any additional targets here.


########################
# This is to handle the case where the statgen library is not already there.
# 
$(LIB_PATH_SAMPLE_PROGRAM)/Makefiles/Makefile.src:
	@echo Unable to locate: $(LIB_PATH_SAMPLE_PROGRAM)/Makefiles/Makefile.src
	@echo To change the location, set LIB_PATH_GENERAL or LIB_PATH_SAMPLE_PROGRAM to the appropriate path to libStatGen.  Or specify \"make LIB_PATH_GENERAL=yourPath\" or \"make LIB_PATH_SAMPLE_PROG=yourPath\"
	@echo Use make download if you have git and want to download the current libStatGen to that location

# Specify to download the statgen library from git (read-only)
download:
	@if test -d $(LIB_PATH_SAMPLE_PROGRAM); \
	then echo $(LIB_PATH_SAMPLE_PROGRAM) already exists; \
	else git clone git://github.com/statgen/libStatGen.git $(LIB_PATH_SAMPLE_PROGRAM); fi
	$(MAKE) -C $(LIB_PATH_SAMPLE_PROGRAM) --no-print-directory all; \
