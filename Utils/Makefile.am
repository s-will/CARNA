bin_SCRIPTS = carna.pl
CLEANFILES =$(bin_SCRIPTS)
EXTRA_DIST =carna.pl.in

do_subst = sed -e 's,[@]CARNAHOMEDIRECTORY[@],$(bindir),g'

carna.pl: carna.pl.in Makefile
	$(do_subst) <  carna.pl.in > carna.pl
	chmod +x carna.pl
