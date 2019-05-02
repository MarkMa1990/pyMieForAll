APPS = MieScatterForAll_V0_1
LIBAPPS = libmiescatterforallv01.so

all: ${LIBAPPS}

${APPS}.o: ${APPS}.cpp
	g++ -c $^ -fPIC -o $@

${LIBAPPS}: ${APPS}.o
	g++ -O3 -o $@ -shared -lcomplex_bessel $^
	rm -r $^

clean:
	rm -r ${LIBAPPS}
	
