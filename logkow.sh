#!/bin/bash

./csm --c $1 --s wat --model $2 | grep "lnGamma(2)" | gawk '{print $4}' > h2o.out

./csm --c $1 --s noco --model $2 | grep "lnGamma(2)" | gawk '{print $4}' > oct.out

paste h2o.out oct.out | gawk '{print log(((2.71828**$1)*8.37)/((2.71828**$2)*55.5))*0.4342944}'
