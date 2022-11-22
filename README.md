# Find Weibull

Simple program for finding the Weibull distribution parameters
k shape factor and c scale factor. Usage:<br><br>
cargo run <i>filename</i><br><br>
where <i>filename</i> is a text file from which program loads the first column
with the header. For example from the main project directory:<br>
cargo run data_example/WS125.txt<br><br>
The file <i>WS125.txt</i> with the example wind measurement data set comes
from the measurement mast US Virgin Islands St. Thomas Bovoni and
was downloaded from the site<br>
<https://midcdmz.nrel.gov/apps/sitehome.pl?site=USVILONA>.

## Information about the data set used
### Any publication based in whole or in part on these data sets should cite the data source as:
Roberts, O.; Andreas, A.; (1997). United States Virgin Islands:<br>
St. Thomas & St. Croix (Data); NREL Report No. DA-5500-64451.<br>
<http://dx.doi.org/10.7799/1183464><br>
<https://midcdmz.nrel.gov/><br><br>
Sorting function from<br>
<https://codereview.stackexchange.com/questions/237790/quick-sort-algorithm-in-rust>
Gamma function from<br>
<https://codereview.stackexchange.com/questions/116850/gamma-function-in-rust>




