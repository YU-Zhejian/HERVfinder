help:
	echo "Use 'make dist' to build this program into a wheel"

.PHONY: doc
doc:
	make -C doc

.PHONY: dist
dist:
	python setup.py sdist bdist_wheel
