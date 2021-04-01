EXEC= bendComposition bilayer build brushes brushesAndNanoparticles \
continuumSphereAndBilayer continuumSphereAndLiposome continuumSphereAndLiposomeWrapped \
continuumSphereGrid continuumSphereSurfactant domain gas lipoCyto liposome \
nanoSphereAndBilayer nanoSphereGrid nanoSphereModBilayer rigidNano scl \
setMDConstants setMDConstants_UminAH solventBilayer solvent tetrahedron tube

CC=g++
CPPFLAGS=-O3
TARGETS: $(EXEC)

ifeq ($(PREFIX),)
PREFIX=$(HOME)
endif

.PHONY: all

all: $(TARGETS)

$(TARGET):
	@echo $(CC) -o $(TARGET) $(TARGET).cpp $(CPPFLAGS)

.PHONY: clean

clean:
	rm -f $(EXEC)

install: $(EXEC)
	install -d $(PREFIX)/bin/
	install -m 755 $(EXEC) $(PREFIX)/bin/
