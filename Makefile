CPPFLAGS=-g -Wall 

BOOST_HOME=/home/yichen.lyh/boost_home
TNT_HOME=/home/yichen.lyh/linear_algebra_library/tnt
JAMA_HOME=/home/yichen.lyh/linear_algebra_library/jama

LIBS= -L /usr/include/
INCLUDE = -I ${TNT_HOME} -I ${JAMA_HOME}

CPLUS_INCLUDE_PATH=${BOOST_HOME}/include
export CPLUS_INCLUDE_PATH

.PHONY : clean all

all: $(subst .cpp,.o,$(SOURCES))  pca


%.O: %.cpp
	$(CXX) $(CPPFLAGS) ${INCLUDE} ${LIBS} $^ $@
pca: pca.cpp
	$(CXX) $(CPPFLAGS) ${INCLUDE} $^  ${LIBS} -o $@

clean:
	rm -rf  *.o  pca
