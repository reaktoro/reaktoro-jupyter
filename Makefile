# Note: Tabs are needed in Makefiles; do not let your editor change tabs to spaces in this file!
# Below is an example.
#
# target:
# ****command args
#
# The "****" above denotes a single tab, not spaces.

all: tutorial notebooks

tutorial: clean notebooks
	mkdir -p tutorial
	for f in scripts/notebooks/*.ipynb ; do \
		python utilities/unpair.py $$f --output tutorial/$${f##*/} ; \
	done

notebooks:
	jupytext --set-formats notebooks//ipynb,py:light scripts/*.py

sync:
	jupytext --sync scripts/*.py

clean:
	rm -rf scripts/notebooks
	rm -rf tutorial/*ipynb
