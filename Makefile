# Note: Tabs are needed in Makefiles; do not let your editor change tabs to spaces in this file!
# Below is an example.
#
# target:
# ****command args
#
# The "****" above denotes a single tab, not spaces.

all: tutorial notebooks

tutorial: clean notebooks
	cp -r scripts/notebooks/*ipynb tutorial/
	jupytext --update-metadata '{"jupytext": null}' tutorial/*.ipynb

notebooks:
	jupytext --set-formats notebooks//ipynb,py:light scripts/*.py
	jupytext --set-formats notebooks//ipynb,md scripts/*.md

sync:
	jupytext --sync scripts/*.py

clean:
	rm -rf scripts/notebooks
	rm -rf tutorial/*ipynb
