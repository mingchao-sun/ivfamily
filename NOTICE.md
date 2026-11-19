# Third-party software and data

ivfamily is distributed under the GNU General Public License, version 3 (see `LICENSE`).

## Code

The estimators and tests are implemented from the published formulas in base R. No code of other packages is copied.
Where a statistic is defined by an algorithm rather than a closed form, the implementation follows the published description and
was checked against the reference implementation:

| Statistic | Reference | Checked against |
|---|---|---|
| Kleibergen-Paap rk Wald F | Kleibergen and Paap (2006); the standardisation and singular-value construction follow the Stata command `ranktest` (Kleibergen, Schaffer, Baum) | `ivreg2` |
| Effective F, simplified critical values | Montiel Olea and Pflueger (2013) | `weakivtest` |
| AR, K, CLR tests and confidence sets | Anderson and Rubin (1949), Kleibergen (2002, 2005), Moreira (2003), Finlay and Magnusson (2009) | `weakiv` |
| Union of confidence intervals, LTZ | Conley, Hansen and Rossi (2012) | `plausexog` |

jmvcore, R6 and ggplot2 are used for the jamovi interface and the plots.

## Example data

`cigarettes_sw.csv`, `college_distance.csv` and `fertility.csv` are derived from the data sets `CigarettesSW`, `CollegeDistance`
and `Fertility` of the CRAN package `AER` (Kleiber and Zeileis, GPL-2 | GPL-3), reproduced here under that licence.
Original sources: Stock and Watson (2007) for the cigarette data; Rouse (1995) for the college data (via Stock and Watson);
Angrist and Evans (1998) for the fertility data. `sim_two_endogenous.csv` and `sim_weak.csv` are simulated for this module.
