CC = gcc
LD = gcc
DEFS = /D=HAVE_CONFIG_H /D=ANSI_C
INTDIR = g:\programs\meschachd_co_complex
PROGRAM = meschachd.dll
CFLAGS = 
LDFLAGS = /LD /Fe $(PROGRAM)

OBJS = \
	$(INTDIR)\bdfactor.obj

SRCS = \
	$(INTDIR)\bdfactor.c 


$(PROGRAM):	$(OBJS)
		@echo "Linking $(PROGRAM) ..."
		$(LD) $(ldflags) $(OBJS)
		@echo "done"

bdfactor.obj:	bdfactor.c
		$(CC) $(CFLAGS) /c bdfactor.c