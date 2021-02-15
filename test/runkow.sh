#!/bin/bash

cd kowa
cp $1.cosmo ../out.cosmo
cd ..
cd kowb
cp $1.cosmo ../out2.cosmo
cd ..
./logkow.sh out $2
./logkow.sh out2 $2
