CC= g++ 
CFLAGS= -c
LDFLAGS=
SRCS= main.cpp alignment/Alignment.cpp alignment/Sequence.cpp general/DirichletMixture.cpp  \
        general/UserParameters.cpp general/Utility.cpp BETE/Node.cpp BETE/BETEAlgorithm.cpp \
        BETE/AgglomerationStep.cpp fileInput/FileReader.cpp fileInput/OutputFile.cpp
OBJECTS= $(SRCS:.cpp=.o)
EXECUTABLE= sciphy
DEPEND= makedepend $(CFLAGS)

all: $(SRCS) $(EXECUTABLE)

$(EXECUTABLE): $(OBJECTS)
	$(CC) $(LDFLAGS) $(OBJECTS) -o $@

.cpp.o:
	$(CC) $(CFLAGS) $< -o $@

depend: $(SRCS)
	$(DEPEND) $(SRCS)

clean:
	rm -rf *.o */*.o $(EXECUTABLE) 

