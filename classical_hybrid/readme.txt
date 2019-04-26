Instructions:
- Download ntruprimeparams.zip from http://scarecryptow.org/publications/ntruprime.html
- Copy the "bkzsim.gp" file to the "over" folder.
- The files used for the overestimetes can be found in the "over" folder, the files used for the undrestimetes in the "under" folder.
- There is one file for each scheme, whose name is starting with the scheme's name.
- Each of those files contains a "find_r" function.
  This function gets as input the scheme's parameters, a list of mitm dimensions r (k_vals), a starting and an ending Hermite delta.
  For each r the function computes the optimal corresponding delta_r and log2 of the runtime of the attack.