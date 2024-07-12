# Data Load For Computational Bio & Pharma - Working Group for Rare Disease

This script loads the data used by the mocomakers meetup  [Computational Bio & Pharma - Working Group for Rare Disease](https://www.meetup.com/mocomakers/events/301844982/).

For access to the input data files, please contact Matt Zamora via `matt@mocomakers.com`. The data files should be placed in the `data` directory.

The output of this python script will upload the 2 data tables into a sqlite database called `rare_disease.db`. The inputs are: 
- `secondary-screen-dose-response-curve-parameters.csv` file that has medication effect on disease and target  
- `OmicsSomaticMutationsMatrixDamaging.csv` file that holds medication and target pivot table indicating level of mutation

It will also create a table that merges the data from both tables and outputs a combined table that holds only the most mutated records for medication effect on disease and target 

LICENSE
-------
<a rel="license" href="http://creativecommons.org/publicdomain/zero/1.0/"> <img src="https://licensebuttons.net/p/zero/1.0/88x31.png" style="border-style: none;" alt="CC0" />  </a>

To the extent possible under law, Jacob Barhak has waived all copyright and 
related or neighboring rights to Data Load For Computational Bio & Pharma - Working Group for Rare Disease
This work is published from: Israel.
