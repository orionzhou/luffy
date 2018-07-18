BOOST = /soft/sle11/boost/1.47.0
LIBSEQ = /soft/sle11/libsequence/1.7.3

CC      = mpiCC 
CFLAGS  = -Wall -g -wd981 -wd383 -wd1418 -wd177 -wd593 -I$(BOOST)/include -I$(LIBSEQ)/include
LIBS    = -lsequence -lboost_iostreams -lboost_regex -lboost_filesystem -lboost_program_options -lboost_random
LDFLAGS = -L$(BOOST)/lib -L$(LIBSEQ)/include
DEP     =
BIN_DIR = bin
BIN_NAME= matrix taskDiv mf
EXE     = $(addprefix $(BIN_DIR)/,$(BIN_NAME))
vpath %.h src
vpath %.cc src
vpath %.o obj
all: $(BIN_DIR) $(EXE) 

bin/matrix: matrix.cc
	$(CC) $(CFLAGS) src/matrix.cc $(LDFLAGS) $(LIBS) -o $@

bin/mf: mf.cc task_division.h
	$(CC) $(CFLAGS) src/mf.cc $(LDFLAGS) $(LIBS) -o $@

bin/taskDiv: taskDiv.cc task_division.h
	$(CC) $(CFLAGS) src/taskDiv.cc $(LDFLAGS) $(LIBS) -o $@

$(BIN_DIR):
	@echo Bin Dir is $(BIN_DIR)
	@if [ ! -d $(BIN_DIR) ]; then mkdir -p $(BIN_DIR); fi

clean:
	rm -f obj/*.o
	rm -f bin/* 

