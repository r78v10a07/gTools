	
all:	clean
	cd gTools; make 
	cd TPMCalculator; make
	cd TPMOutPutMerger; make

test:	all
	cd gTools; make test 
	cd TPMCalculator; make test
	cd TPMOutPutMerger; make test

clean:
	cd TPMOutPutMerger; make clean
	cd TPMCalculator; make clean
	cd gTools; make clean
	
	
