../configure: aclocal.m4 configure.ac
	/bin/rm -rf autom4te.cache
	autoconf && mv ./configure ..
update: config.sub config.guess

config.guess:
	wget -O config.guess 'http://git.savannah.gnu.org/gitweb/?p=config.git;a=blob_plain;f=config.guess;hb=HEAD'
config.sub:
	wget -O config.sub 'http://git.savannah.gnu.org/gitweb/?p=config.git;a=blob_plain;f=config.sub;hb=HEAD'

.PHONY: config.guess config.sub