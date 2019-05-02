
APPS = pyMie
DEPENDENT_LIBS = ../libmiescatterforallv01.so

ALL: ${APPS}

%: %.cpp
	swig -python -c++ $@.i
	g++ -c -fPIC -I/usr/include/python2.7 $^ $@_wrap.cxx
	g++ -shared $@.o $@_wrap.o ${DEPENDENT_LIBS} -O3 -lcomplex_bessel -o _$@.so

	rm -r $@.o $@_wrap.o $@_wrap.cxx 

clean:
	rm -r _${APPS}.so
