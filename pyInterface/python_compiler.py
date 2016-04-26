
import compileall
import os
import sys


if __name__ == "__main__":

	indir = sys.argv[1]
	destdir = sys.argv[2]

	compileall.compile_dir(dir=indir, ddir=destdir, quiet=True)

	os.system('rm -rf ' + destdir + '/*')
	os.system('cd ' + indir + '; for dir in `find * -type d`; do mkdir ' + destdir + '/$dir; done; cd - > /dev/null')
	os.system('cd ' + indir + '; for pycfile in `find . -name "*.pyc"`; do mv $pycfile ' +  destdir + '/$pycfile; done; cd - > /dev/null')
	os.system('cd ' + indir + '; for pyfile in `find . -name "*.py"`; do cp $pyfile ' + destdir + '/$pyfile; done; cd - > /dev/null')
