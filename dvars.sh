#!/bin/bash
#
# Script: DVARS.sh
# Purpose: Create standardized version of DVARS
# Author: T. Nichols
# Version: http://github.com/nicholst/DVARS/commit/5814337
#          2017-02-19 08:18:50 +0000
#


###############################################################################
#
# Environment set up
#
###############################################################################

shopt -s nullglob # No-match globbing expands to null
TmpDir=/tmp
Tmp=$TmpDir/`basename $0`-${$}-
trap CleanUp INT

###############################################################################
#
# Functions
#
###############################################################################

Usage() {
cat <<EOF
Usage: `basename $0` [options] fMRI_4Dtimeseries DVARSout
       `basename $0` [options] fMRI_1 fMRI_2 ... fMRI_T DVARSout

Options
    -all Produce 2 additional versions of DVARS.

Creates standardized version of DVARS, normalizing according to the expected 
value of DVARS under an AR1 model. DVARSout, consists of plain text file with
this normalized version of DVARS, scaled so that it is approximately 1 if  
there are no artifacts. 

With the -all option, 2 additional columns are added to DVARSout.  The 2nd column
is the raw DVARS with no scaling (precisely the root mean square (RMS) of the temporal 
difference).  The 3rd is the precision-normalized DVARS:  Before taking the RMS, 
temporal difference images are standardized voxel-wise giving a more precisely
normalized DVARS measure.  A side effect, however, is that high-variance parts of 
the image are down-weighted relative to low-variance areas.
_________________________________________________________________________
\$Id: 22c29d73383968ef25d48d60ad64291a22f27736 $
EOF
exit
}

CleanUp () {
    /bin/rm -f ${Tmp}*
    exit 0
}


###############################################################################
#
# Parse arguments
#
###############################################################################

while (( $# > 1 )) ; do
    case "$1" in
        "-m")
            shift
            brain_mask="$1"
            echo "Using ${brain_mask} as mask"
            break
            ;;
        "-help")
            Usage
            ;;
        "-all")
            shift
            AllVers="$1"
            ;;
        -*)
            echo "ERROR: Unknown option '$1'"
            exit 1
            break
            ;;
        *)
            break
            ;;
    esac
done
Tmp=$TmpDir/DVARS-${$}-

if (( $# < 2 )) ; then
    Usage
elif (( $# == 2 )) ; then
    FUNC="$1"
    OUT="$2"
else
    ((i=0))
    Imgs=()
    while (( $# >= 2 )) ; do
	Imgs[i]="$1"
	((i++))
	shift
    done
    FUNC=$Tmp-1
    fslmerge -t $FUNC "${Imgs[@]}"
    OUT="$1"
fi


###############################################################################
#
# Script Body
#
###############################################################################

echo -n "."

# Find mean over time
fslmaths "$FUNC" -Tmean $Tmp-Mean
# Find the brain

if [ -z $brain_mask ]
then
    brain_mask=$Tmp-MeanBrain
    bet $Tmp-Mean  $brain_mask
fi

echo "Using $brain_mask as mask"

# Compute robust estimate of standard deviation
fslmaths "$FUNC" -Tperc 25 $Tmp-lq -odt float
fslmaths "$FUNC" -Tperc 75 $Tmp-uq -odt float
fslmaths $Tmp-uq -sub $Tmp-lq -div 1.349 $Tmp-SD -odt float

# Compute (non-robust) estimate of lag-1 autocorrelation
fslmaths "$FUNC" -sub $Tmp-Mean -Tar1 $Tmp-AR1 -odt float

# Compute (predicted) standard deviation of temporal difference time series
fslmaths $Tmp-AR1 -mul -1 -add 1 -mul 2 -sqrt -mul $Tmp-SD  $Tmp-DiffSDhat -odt float

# Save mean value
DiffSDmean=$(fslstats $Tmp-DiffSDhat -k $brain_mask -M)

echo -n "."

# Compute temporal difference squared time series
nVol=$(fslnvols "$FUNC")
fslroi "$FUNC" $Tmp-FUNC0 0 $((nVol-1))
fslroi "$FUNC" $Tmp-FUNC1 1 $nVol
fslmaths $Tmp-FUNC0 -sub $Tmp-FUNC1 -sqr $Tmp-DiffSq -odt float

echo -n "."

# Compute DVARS, no standization
fslstats -t $Tmp-DiffSq -k $brain_mask -m > $Tmp-DiffVar.dat

if [ "$AllVers" = "" ] ; then
    # Standardized
    awk '{printf("%g\n",sqrt($1)/'"$DiffSDmean"')}' $Tmp-DiffVar.dat > "$OUT"
else
    # Compute DVARS, based on voxel-wise standardized image
    fslmaths $Tmp-FUNC0 -sub $Tmp-FUNC1 -div $Tmp-DiffSDhat -sqr $Tmp-DiffSqVxStdz
    fslstats -t $Tmp-DiffSqVxStdz -k $brain_mask -M | awk '{print sqrt($1)}' > $Tmp-DiffVxStdzSD.dat

    # Sew it all together
    awk '{printf("%g\t%g\n",sqrt($1)/'"$DiffSDmean"',sqrt($1))}' $Tmp-DiffVar.dat > $Tmp-DVARS
    paste $Tmp-DVARS $Tmp-DiffVxStdzSD.dat > "$OUT"
fi

echo -n "."

echo "."

###############################################################################
#
# Exit & Clean up
#
###############################################################################

#CleanUp

