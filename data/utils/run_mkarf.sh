#!/bin/sh
# This script creates the reference arfs.  First it runs marx to create an
# on-axis simulation.  From that, a pcad file is created.
# Then the reference arfs are created using the pcad and marx event file as
# input.

runciao="./runciaotool.sh"
# $runciaox is a script that runs a ciaotool in its own environment

outdir="point"
datadir=".."
bindir="../../src/objs"
pardir="../../par"
marx="$bindir/marx"
marx2fits="$bindir/marx2fits"
marxasp="$bindir/marxasp"
evtfile="marx.evt"
asolfile="marx.asol"
arfprefix="hrma"
aspfile="marx.asp"

cp $pardir/marx.par .
cp $pardir/marxasp.par .

# Run the simulation and create the evt, and asol files.
$marx OutputDir=$outdir DataDirectory=$datadir \
 SourceRA=0 SourceDec=0 RA_Nom=0 Dec_Nom=0 Roll_Nom=0 DitherModel=INTERNAL \
 SpectrumType=FLAT SourceFlux=0.001 MinEnergy=1.5 MaxEnergy=1.5 \
 ExposureTime=10000 GratingType=NONE DetectorType=ACIS-S \
 ACIS_Exposure_Time=3.2 ACIS_Frame_Transfer_Time=0.0

echo $marx2fits "$outdir" "$evtfile"
$marx2fits "$outdir" "$evtfile"
echo $marxasp OutputFile="$asolfile" MarxDir="$outdir"
$marxasp OutputFile="$asolfile" MarxDir="$outdir"

# Now create the aspect histogram and the arfs
detname="ACIS-7;UNIFORM;ideal;bpmask=0"
grating="NONE"
engrid="0.03:12.0:0.003"
obsfile="$evtfile"

echo asphist infile="$asolfile" outfile="$aspfile" evtfile="$evtfile" \
  verbose=0 mode=h clobber=yes
$runciao asphist infile="$asolfile" outfile="$aspfile" evtfile="$evtfile" \
  verbose=0 mode=h clobber=yes

X=4096.5 Y=4096.5

for shell in 1 3 4 6
do
  outfile="$arfprefix$shell.arf"
  echo mkarf mirror="hrma;shell=$shell" detsubsys="$detname" grating="$grating" \
    outfile="$outfile" engrid="$engrid" asphistfile="$aspfile" \
    sourcepixelx=$X sourcepixely=$Y obsfile="$evtfile" \
    maskfile=NONE pbkfile=NONE dafile=NONE verbose=1 mode=h clobber=yes
  $runciao mkarf mirror="hrma;shell=$shell" detsubsys="$detname" grating="$grating" \
    outfile="$outfile" engrid="$engrid" asphistfile="$aspfile" \
    sourcepixelx=$X sourcepixely=$Y obsfile="$evtfile" \
    maskfile=NONE pbkfile=NONE dafile=NONE verbose=1 mode=h clobber=yes
done
