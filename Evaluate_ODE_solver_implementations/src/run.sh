#!/bin/sh

cmake CMakeLists.txt 

make clean
make

./izhikevich -nRS -cstep150 -fsvg/RS_step150.svg -mstdEuler -msympEuler -mSpiNNaker -mEulerNEST -mIzhikevichNEST -mGSL -mHeun -mhighResEuler -madaptiveEuler
# ./izhikevich -nRS -cstep150 -fsvg/RS_step150.svg -mhighResEuler -mstdEuler
#./izhikevich -nRS -cstep150 -fsvg/RS_step150.svg -mhighResEuler -msympEuler
#./izhikevich -nRS -cstep150 -fsvg/RS_step150.svg -mhighResEuler -mSpiNNaker
#./izhikevich -nRS -cstep150 -fsvg/RS_step150.svg -mhighResEuler -mEulerNEST
#./izhikevich -nRS -cstep150 -fsvg/RS_step150.svg -mhighResEuler -mIzhikevichNEST
#./izhikevich -nRS -cstep150 -fsvg/RS_step150.svg -mSpiNNaker -mGSL
#./izhikevich -nRS -cstep150 -fsvg/RS_step150.svg -mhighResEuler -madaptiveEuler
#./izhikevich -nRS -cstep150 -fsvg/RS_step150.svg -mSpiNNaker
#./izhikevich -nCH -cstep150 -fsvg/RS_step150.svg -mhighResEuler -madaptiveEuler -mGSL
#./izhikevich -nRS -cstep150 -fsvg/RS_step150.svg -mGSL -madaptiveEuler -mhighResEuler

#./izhikevich -nRS -cstep150 -fsvg/RS_GSL.svg -mGSL
#./izhikevich -nRS -cstep150 -fsvg/RS_SpiNNaker.svg -mSpiNNaker
#./izhikevich -nCH -cstep150 -fsvg/CH_Type_IntStepSize0.5ms.svg -mstdEuler -mGSL -madaptiveEuler
#./izhikevich -nRS -cstep150 -fsvg/RS_step150.svg -mstdEuler

#./izhikevich -nIB -cstep150 -fsvg/IB_step150.svg
#./izhikevich -nCH -cstep150 -fsvg/CH_step150.svg
#./izhikevich -nFS -cstep150 -fsvg/FS_step150.svg
#./izhikevich -nTC -cstep150 -fsvg/TC_step150.svg
#./izhikevich -nRZ -cstep150 -fsvg/RZ_step150.svg
#./izhikevich -nLTS -cstep150 -fsvg/LTS_step150.svg
