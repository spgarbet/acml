TOP=../..
PKG=${shell cd ${TOP};pwd}
SUITE=run_tests.R
R=R

all: inst test

inst: # Install package
	cd ${TOP}/..;\
	${R} CMD INSTALL ${PKG}

test: # Run unit tests
	export RCMDCHECK=FALSE;\
	${R} --vanilla --slave < ${SUITE}
