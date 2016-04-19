# Put the directory above 'DNest4' into an environment variable
# called DNEST4_PATH for this to work.

CC = g++
CXXFLAGS = -std=c++11 -O3 -Wall -Wextra -pedantic -DNDEBUG -I/usr/local/include
LIBS = -ldnest4 -lpthread -lgsl -lgslcblas -lm
LINKLIB = -L/usr/local/lib

default:
	$(CC) -I$(DNEST4_PATH) $(CXXFLAGS) -c *.cpp
	$(CC) -L$(DNEST4_PATH)/DNest4/code $(LINKLIB) -o main *.o $(LIBS)
	rm *.o
