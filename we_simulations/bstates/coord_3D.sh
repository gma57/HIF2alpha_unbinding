#!/bin/bash


for BS in 01 02 03 04 05
do

    echo $BS
    python calculateSASA_pcoord.py test.out ${BS}/bstate.ncrst ../common_files/bound_new.prmtop
    cat test.out | tail -n +1 | awk {'print $1'} > test_test.out
    
    COMMAND="         parm ../common_files/bound_new.prmtop\n"
    COMMAND="$COMMAND trajin ${BS}/bstate.ncrst\n"
    COMMAND="$COMMAND reference ../common_files/reference.pdb\n"
    COMMAND="$COMMAND autoimage \n"
    COMMAND="$COMMAND strip :WAT,Na+,Cl- \n"
    COMMAND="$COMMAND rms reference :6-115 \n"
    COMMAND="$COMMAND rms reference @1-8,10,12,13,15,17-19,22,24-26,29 nofit out rmsd.dat \n"
    COMMAND="$COMMAND go\n"
    echo -e $COMMAND | $CPPTRAJ
    cat rmsd.dat | tail -n +2 | awk {'print $2'} > rmsd.txt 
    
    COMMAND="         parm ../common_files/bound_new.prmtop\n"
    COMMAND="$COMMAND trajin ${BS}/bstate.ncrst\n"
    COMMAND="$COMMAND autoimage \n"
    COMMAND="$COMMAND strip :WAT,Na+,Cl- \n"
    COMMAND="$COMMAND nativecontacts mindist :1 :2-118 out dist.dat noimage \n"
    COMMAND="$COMMAND go\n"
    echo -e $COMMAND | $CPPTRAJ
    cat dist.dat | tail -n +2 | awk {'print $4'} > dist.txt 
    
    
    paste test_test.out rmsd.txt dist.txt | awk {'print $1 , $2 , $3'} > pcoord.init
    rm rmsd.txt test.out test_test.out rmsd.dat dist.dat dist.txt
    mv pcoord.init ${BS}/
    
done
