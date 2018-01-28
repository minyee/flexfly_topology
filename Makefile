TARGET := libflexfly.so 
SRC := flexfly_topology_simplified.cc flexfly_optical_network.cc flexfly_routing_algo.cc network_manager_node.cc connectivity_matrix.cc my_logp_switch.cc flexfly_topology.cc flexfly_optical_switch.cc flexfly_electrical_switch.cc 

CXX :=    libsst++
CC :=     libsstcc
CXXFLAGS := -fPIC
CPPFLAGS := -I. -I$(HOME)/Programs/install/sst-core/clang/include
#CPPFLAGS := -I. -I$(HOME)/install/sstcore-7.1.0/include
LIBDIR :=  
PREFIX := 
LDFLAGS :=  -Wl,-rpath,$(PREFIX)/lib

OBJ := $(SRC:.cc=.o) 
OBJ := $(OBJ:.cpp=.o)
OBJ := $(OBJ:.c=.o)

.PHONY: clean install 

all: $(TARGET)

$(TARGET): $(OBJ) 
	$(CXX) -o $@ $+ $(LDFLAGS) $(LIBS)  $(CXXFLAGS)

%.o: %.cc 
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -c $< -o $@

%.o: %.c
	$(CC) $(CPPFLAGS) $(CFLAGS) -c $< -o $@

clean: 
	rm -f $(TARGET) $(OBJ) 

install: $(TARGET)
	cp $< $(PREFIX)/bin

