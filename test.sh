#!/bin/bash

echo "Calculating potentials in H2O"

echo "Methanol"
cp meo.energy energy
./csm --c meo --s h2o | tail -3
echo "Methan"
cp met.energy energy
./csm --c met --s h2o | tail -3
echo "Hexane"
cp hex.energy energy
./csm --c hex --s h2o | tail -3
echo "Octane"
cp oct.energy energy
./csm --c oct --s h2o | tail -3
echo "Hydrogen"
cp hyd.energy energy
./csm --c hyd --s h2o | tail -3
