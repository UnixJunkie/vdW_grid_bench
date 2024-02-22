
%.cmx : %.ml
	ocamlopt -c $^

grid_init: pqrs.cmx V3.cmx ptable.cmx UFF.cmx mol.cmx grid_init.cmx
	ocamlopt -o $@ $^


test: grid_init
	time ./grid_init -dx 0.1 -l ../data/lig.pqrs -p ../data/prot.pqrs

clean:
	rm -f *.cmo *.cmi *.cmx
