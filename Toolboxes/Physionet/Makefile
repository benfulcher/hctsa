BINDIR=~/bin

all: sampen

# install to user's binary directory
install:
	echo $(BINDIR)
	mv sampen $(BINDIR)

# 'make sampen': compile the program in the current directory
sampen: sampen.c
	$(CC) $(CFLAGS) -o sampen sampen.c -lm

# 'make clean': removes the compiled file from the directory
clean:
	rm sampen
