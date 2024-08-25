TARGET = ft8d

OBJECTS = pffft.o ft8d.o

CC = gcc
LD = gcc
RM = rm -f

CFLAGS = -Wall -O3 -funroll-loops -march=armv7-a -mtune=cortex-a9 -mfpu=neon -mfloat-abi=hard -ffast-math -fsingle-precision-constant -mvectorize-with-neon-quad
LDFLAGS = -lm

all: $(TARGET)

%.o: %.c
	${CC} -c ${CFLAGS} $< -o $@

$(TARGET): $(OBJECTS)
	$(LD) $(OBJECTS) $(LDFLAGS) -o $@

clean:
	$(RM) *.o $(TARGET)
