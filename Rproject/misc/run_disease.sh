name="breast_expr"
dir0="/hpf/largeprojects/agoldenb/aziz/testing_expression/"
k0=100
exact=0

qsub -l mem=40G,vmem=40G,nodes=1:ppn=1,walltime=71:50:00 -N $name"_$k0" -v dir0=$dir0,name="$name",k0="$k0",exact="$exact" /hpf/largeprojects/agoldenb/aziz/testing_expression/run_disease_job.sh

