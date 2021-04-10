HISTOGRAMS = PQS_histogram.png forward_histogram.png
PLOTS =  OU_PQS.png autocorrelations.png OU_forward.png OU_backward.png $(HISTOGRAMS)

.PHONY:default clean

default:$(PLOTS) Makefile

OU_PQS.png:OU_PQS.gpi OU_PQS_data.txt
	gnuplot $<

autocorrelations.png:autocorrelations.gpi autocorrelations.txt
	gnuplot $<

OU_forward.png:OU_forward.gpi OU_forward_data.txt
	gnuplot $<

OU_backward.png:OU_backward.gpi OU_backward_data.txt
	gnuplot $<

#scatterplot.png:scatterplot.gpi scatterplot_data.txt
#	gnuplot $<

$(HISTOGRAMS):histograms.py scatterplot_data.txt gridnumber.txt dt.txt
	python3 histograms.py
clean:
	$(RM) $(PLOTS)
