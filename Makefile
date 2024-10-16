all:
	cd tests; make

publish:
	publish-signalsmith-raw /tmp/wrapped-fft/
	publish-signalsmith-git /tmp/wrapped-fft.git