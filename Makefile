CFLAGS=-std=c11 -D_POSIX_C_SOURCE=200809L -lm

tsp: *.c *.h
	gcc $(CFLAGS) *.c -o $@

tsp-debug: *.c *.h
	gcc -g $(CFLAGS) *.c -o $@

tsp-release: *.c *.h
	gcc -O3 $(CFLAGS) *.c -o $@

clean:
	rm -f tsp tsp-debug tsp-release
