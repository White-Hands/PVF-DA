target:
	gcc -o  out 426.c -I ~/.local/include/pbc -L ~/.local/lib -Wl,-rpath ~/.local/lib  -l pbc -lgmp
