help:
	echo "Use 'make dist' to build this program into a wheel"

doc:
	make -C doc

dist:
	python setup.py sdist bdist_wheel
