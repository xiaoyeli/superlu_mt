#!/bin/bash

export ofile=pctest.out			# output file
if [ -f $ofile ]; then
    rm -f $ofile
fi
echo 'Single Complex-precision testing output' > $ofile

export NVAL=(9 19)
export NRHS=(5)
export LWORK=(0 100000000)
export PANELSIZE=(2)
export RELAX=(2)
export NPROCS=(1 4)

export OMP_NUM_THREADS=4

#
# Loop through all matrices ...
#
for m in LAPACK cg20.cua; do
echo $m
    #--------------------------------------------
    # Test matrix types generated in LAPACK-style
    #--------------------------------------------
    if [ $m = LAPACK ]; then
      	echo '' >> $ofile
      	echo '** LAPACK test matrices' >> $ofile
      	for n in $NVAL ; do
            for s in $NRHS ; do
              	for l in $LWORK ; do
		    for p in $NPROCS ; do
	    	      	echo '' >> $ofile
            	      	echo 'n='$n 'nrhs='$s 'lwork='$l 'nprocs='$p >> $ofile
            	      	./pctest -t "LA" -l $l -n $n -s $s -p $p >> $ofile
		    done
              	done
            done
        done
    #--------------------------------------------
    # Test a specified sparse matrix
    #--------------------------------------------
    else
      	echo '' >> $ofile
      	echo '** sparse matrix:' $m >> $ofile
      	for  w in $PANELSIZE ; do
            for r in $RELAX ; do
                for s in $NRHS ; do
                    for l in $LWORK ; do
			for p in $NPROCS ; do
	                    echo '' >> $ofile
                      	    echo 'w='$w 'relax='$r 'nrhs='$s 'lwork='$l \
					'nprocs='$p >> $ofile
            	      	    ./pctest -t "SP" -w $w -r $r -s $s -l $l \
                             		< ../EXAMPLE/$m >> $ofile
	   		done
        	    done
		done
	    done
	done
   fi
done



