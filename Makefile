#
# Makefile for ptaSimulate
#
# Set INCLUDES and LDFLAGS variables to suit your system
#

PREFIX := /usr/local/
PSR_PREFIX := /opt/pulsar/

CC := gcc

CFLAGS := -lm -g -Wall

INCLUDES := -I$(PREFIX)/include -I$(PSR_PREFIX)/include

LFLAGS := -L$(PREFIX)/lib -L$(PSR_PREFIX)/lib -L$(PREFIX)/lib64 -L/usr/lib64

LIBS := -lfftw3 -lfftw3f

SRCS := $(wildcard *.c)

OBJS := ${SRCS:.c=.o}

MAIN := ptaSimulate

.PHONY: all clean install

all: $(MAIN)

$(MAIN): $(OBJS)
	$(CC) $(CFLAGS) $(INCLUDES) -o $(MAIN) $(OBJS) $(LFLAGS) $(LIBS)

.c.o:
	$(CC) $(CFLAGS) $(INCLUDES) -c $< -o $@

install:
	install -m 0755 $(MAIN) $(PREFIX)/bin
 
clean:
	$(RM) $(MAIN)
	$(RM) $(OBJS)

