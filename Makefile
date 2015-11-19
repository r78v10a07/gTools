	
all:	clean
	cd gTools; make 
	cd TPMCalculator; make
	cd TPMOutPutMerger; make

clean:
	cd TPMOutPutMerger; make clean
	cd TPMCalculator; make clean
	cd gTools; make clean
	
	
