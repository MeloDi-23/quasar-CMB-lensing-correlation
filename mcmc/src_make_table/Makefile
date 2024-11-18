LIB = -lm

CC = gcc

OBJS01 = create_xitable.o ftio.o 
create_xitable: $(OBJS01)
	$(CC) -o $@ $(OBJS01) $(LIB)

OBJS01a= create_wptable.o ftio.o
create_wptable: $(OBJS01a)
	$(CC) -o $@ $(OBJS01a) $(LIB)

OBJS02 = wf.o 
wf: $(OBJS02)
	$(CC) -o $@ $(OBJS02) $(LIB)

OBJS03 = read2.o ftio.o
read2: $(OBJS03)
	$(CC) -o $@ $(OBJS03) $(LIB)

OBJS04 = read_table.o
read_table: $(OBJS04)
	$(CC) -o $@ $(OBJS04) $(LIB)

OBJS04a = read_wptable.o
read_wptable: $(OBJS04a)
	$(CC) -o $@ $(OBJS04a) $(LIB)

OBJS04b = read_xitable.o
read_xitable: $(OBJS04b)
	$(CC) -o $@ $(OBJS04b) $(LIB)

all: create_xitable create_wptable read2 read_table read_wptable read_xitable

clean:
	rm -f *.o
