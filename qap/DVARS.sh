#!/bin/bash
#
# Script: DVARS.sh
# Purpose: Create standardized version of DVARS
# Author: T. Nichols
# Version: $Id: DVARS.sh,v 1.2 2012/10/26 22:17:19 nichols Exp $
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
Usage: `basename $0` [options] mask fMRI_4Dtimeseries DVARSout
       `basename $0` [options] mask fMRI_1 fMRI_2 ... fMRI_T DVARSout

Options
    -all Produce 2 additional versions of DVARS.

Creates standardized version of DVARS, normalizing according to the expected 
standard deviation of DVARS under an AR1 model. DVARSout, consists of plain text 
file with a normalized version of DVARS, scaled so that it is approximately 1 if 
there are no artifacts. 

With the -all option, 2 additional columns are added to DVARSout.  The 2nd column
is the raw DVARS with no scaling (precisely the standard deviation of the temporal 
difference).  The 3rd is the precision-normalized DVARS:  Before taking the SD, 
temporal difference images are standardized voxel-wise giving a more precisely
normalized DVARS measure.  A side effect, however, is that high-variance parts of 
the image are down-weighted relative to low-variance areas.
_________________________________________________________________________
\$Id: DVARS.sh,v 1.2 2012/10/26 22:17:19 nichols Exp $
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

if (( $# < 3 )) ; then
    Usage
elif (( $# == 3 )) ; then
    MASK="$1"
    FUNC="$2"
    OUT="$3"
else
    MASK="$1"
    shift
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

echo $Tmp
echo -n "."

# Find mean over time
fslmaths "$FUNC" -Tmean $Tmp-Mean
# Find the brain
#bet $Tmp-Mean $Tmp-MeanBrain
#MASK="$Tmp-MeanBrain"

# Compute robust estimate of standard deviation
fslmaths "$FUNC" -Tperc 25 $Tmp-lq
fslmaths "$FUNC" -Tperc 75 $Tmp-uq
fslmaths $Tmp-uq -sub $Tmp-lq -div 1.349 $Tmp-SD -odt float

# Compute (non-robust) estimate of lag-1 autocorrelation
fslmaths "$FUNC" -sub $Tmp-Mean -Tar1 $Tmp-AR1 -odt float

# Compute (predicted) standard deviation of temporal difference time series
fslmaths $Tmp-AR1 -mul -1 -add 1 -mul 2 -sqrt -mul $Tmp-SD  $Tmp-DiffSDhat

# Save mean value
DiffSDmean=$(fslstats $Tmp-DiffSDhat -k $MASK -M)

echo -n "."

# Compute temporal difference time series
nVol=$(fslnvols "$FUNC")
fslroi "$FUNC" $Tmp-FUNC0 0 $((nVol-1))
fslroi "$FUNC" $Tmp-FUNC1 1 $nVol

echo -n "."

# Compute DVARS, no standization
fslmaths $Tmp-FUNC0 -sub $Tmp-FUNC1 -mas $MASK $Tmp-Diff -odt float
#fslstats -t $Tmp-Diff -mas $MASK -S > $Tmp-DiffSD.dat
fslstats -t $Tmp-Diff -S > $Tmp-DiffSD.dat

echo -n '.'

if [ "$AllVers" = "" ] ; then
    # Standardized
    awk '{printf("%g\n",$1/'"$DiffSDmean"')}' $Tmp-DiffSD.dat > "$OUT"
else
    # Compute DVARS, based on voxel-wise standardized image
    fslmaths $Tmp-FUNC0 -sub $Tmp-FUNC1 -div $Tmp-DiffSDhat -mas $MASK $Tmp-DiffVxStdz
    #fslstats -t $Tmp-DiffVxStdz -k $MASK -S > $Tmp-DiffVxStdzSD.dat
    fslstats -t $Tmp-DiffVxStdz -S > $Tmp-DiffVxStdzSD.dat

    # Sew it all together
    awk '{printf("%g\t%g\n",$1/'"$DiffSDmean"',$1)}' $Tmp-DiffSD.dat > $Tmp-DVARS
    paste $Tmp-DVARS $Tmp-DiffVxStdzSD.dat > "$OUT"
fi

echo -n "."

echo "."

###############################################################################
#
# Exit & Clean up
#
###############################################################################

CleanUp
