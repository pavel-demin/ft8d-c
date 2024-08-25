TARGET = ft8d

OBJECTS = pffft.o ft8d.o

CC = gcc
LD = gcc
RM = rm -f

CFLAGS = -O3 -Wall
LDFLAGS = -lm

all: $(TARGET)

%.o: %.c
	${CC} -c ${CFLAGS} $< -o $@

$(TARGET): $(OBJECTS)
	$(LD) $(OBJECTS) $(LDFLAGS) -o $@

clean:
	$(RM) *.o $(TARGET)
