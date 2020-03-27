#!/bin/bash
module load R/3.2.3
dir0=${dir0}
name=${name}
k0=${k0}
exact=${exact}
/hpf/largeprojects/agoldenb/aziz/testing_expression/main_disease.r $dir0 $name $k0 $exact
