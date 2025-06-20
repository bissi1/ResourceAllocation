# ResourceAllocation

This project is based on the theory and application of resource allocation optimization models within the context of two-stage stochastic programming. We are updating this repository as our research progresses.

Please refer to the article *"Release Immediately or Sequentially? Strategies for Allocating Scarce Therapeutic Resources during Disease Outbreaks"* by B. Singh and S. Rebennack, and its associated appendix, that will be appearing soon in [IISE Transactions](https://www.tandfonline.com/journals/uiie21). We thank Oliver Rehberg for the organization of the data files. 
For questions, please contact [Bismark Singh](mailto:b.singh@southampton.ac.uk).

In the `Data` folder, there are 2400 spreadsheets denoting the 1200 sampled scenarios considered in the numerical experiments of the article (there are two files for each scenario – the population $S_{c,t}^\omega$ and the benefit $B_{c,t}^\omega$). Within these files, the columns `"c"` correspond to the 254 counties of Texas while the rows `"t"` correspond to the 15 months. The filenames explain whether these are the benefits (e.g., `infectious-001_c4_ic2_benefit_monthly`) or the populations (e.g., `infectious-001_c4_ic2_population_monthly`).

The filenames contain six pandemic instances (`c1`, `c2`, `c3`, `c4`, `c5`, `c6` correspond to 1918, 1928, 1957, 1968, 2009, 2020, respectively – see Table S2), two initial cases (`ic1` and `ic2` refer to Random and Random-weighted, respectively – see Table S3), and a counter index (from `000` to `099`).

An example GAMS code is provided in the file `example_code.gms`.
