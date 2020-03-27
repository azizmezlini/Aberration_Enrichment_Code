This folder contains all the code we used to process the data in our paper. 
To simply use our test on your own processed variable you only need the file Rproject/functions.r. Or look for our better documented CRAN package “aziztest”.
The other files and folders contain code that we used for simulations, real data processing and figure generation. We are posting it here for reproducibility purposes.
To use our test, first load the file Rproject/functions.r:  source("Rproject/functions.r")
Then simply call: 

aziz.test(y,x)

Where "x" is the vector containing your numerical variable to be tested and "y"" is your case/control labels 
(0 for controls, 1 for cases).
Additional options include:

rep: number of repetitions, 100000 by default.

w: the weighting scheme. NULL by default means the weight are absolute Z-scores from x, 1 means all weights are equal to one.

unidirectional: indicating whether both direction are tested (0, default), only testing for cases>controls (1) or controls>cases (-1).

doall: perform all permutation (TRUE) or only perform the required number to get an acceptable accuracy on the pvalue estimate. Default=FALSE is faster when testing a large number of variables

Example: res=aziz.test(y,x,w=NULL,rep=1000, doall=TRUE,unidirectional=0)

The test returns the following attributes:
pval: pvalue of association computed from permutations.
es: Max standardized enrichment score (our test statistic)
direction: 1 means cases<controls , 2 means cases>controls
r: proportion of cases in interval of interest
esm: Max standardized enrichment score in both directions
esi: Rank of the sample corresponding to the max standardized enrichment score in both directions
