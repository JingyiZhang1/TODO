# A Triple-Outcome Double-Criterion Optimal (TODO) Design

R codes to implement triple-outcome double-criterion optimal design in multi-dose randomized trials.

# Description

Detecting the efficacy signal and determining the optimal dose are critical steps to increase the probability of success and expedite the drug development in cancer treatment. After identifying a safe dose range through phase I studies, conducting a multi-dose randomized trial becomes an effective approach to achieve this objective. However, there have been limited formal statistical designs for such multi-dose trials, and dose selection in practice is often ad hoc, relying on descriptive statistics. We propose a Bayesian optimal two-stage design to facilitate rigorous dose monitoring and optimization. Utilizing a flexible Bayesian dynamic linear model for the dose-response relationship, we employ dual criteria to assess dose admissibility and desirability.
%and determine the optimal dose. Additionally, we introduce a triple-outcome trial decision procedure to consider dose selection beyond clinical factors. Under the proposed model and decision rules, we develop a systematic calibration algorithm to determine the sample size and Bayesian posterior probability cutoffs to optimize specific design operating characteristics. Furthermore, we demonstrate how to concurrently assess toxicity and efficacy within the proposed framework using a utility-based risk-benefit trade-off.  To validate the effectiveness of our design, we conduct extensive simulation studies across a variety of scenarios, demonstrating its robust operating characteristics.

# Functions

The repository includes two functions:

- todo2_data_generation.R: The R code that includes the function ```todo2_data_generation()``` to obtain the simulated data and model fitting results of the TODO design with 2 doses.
  
  ```rscript
  todo2_data_generation (rseed, iscena, ntrial, m1, nsample, peff.m, theta0, delta1, sigma1, tau)
  ```

- todo2_futility.R: The R code that includes the function ```todo2_futility()``` to obtain the futility monitoring results of the TODO design with 2 doses.
  
  ```rscript
  todo2_futility (ia12, fdata, ia.post, fa.post1, fa.post2, fa.post3, post, ntrial, peff, theta0, m1, nsample)
  ```

- todo2.R: The R code that includes the function ```todo2()``` to obtain the operating characteristics of the proposed design by simulating trials.
  
  ```rscipt
  todo2 (ic12, fdata, ia.post,fa.post1, fa.post2, fa.post3,  post, ia.mean, fa.mean1, fa.mean2, fa.mean3, ntrial, peff, theta0, delta1, ia12, idr.m, wl, m1, nsample)
  ```

- TODO2-m1-n-a1-a2-optimizaiton.R: The R code for optimizing parameters (m1, n, a1, a2).

- TODO2-c1-c2-optimization.R: The R code for optimizing parameters (c1, c2).

- TODO2-simulation.R: The R code for obtaining the operating characteristics of the proposed design by simulating trials.

- todo3_data_generation.R: The R code that includes the function ```todo3_data_generation()``` to obtain the simulated data and model fitting results of the TODO design with 3 doses.
  
  ```rscript
  todo3_data_generation (rseed, iscena, ntrial, m1, nsample, peff.m, theta0, delta1,sigma1,tau)
  ```

- todo3_futility.R: The R code that includes the function ```todo3_futility()``` to obtain the futility monitoring results of the TODO design with 3 doses.
  
  ```rscript
  todo3_futility (ia12, fdata, ia.post, fa.post1, fa.post2, fa.post3, post1, post2, ntrial, peff, theta0, m1, nsample)
  ```

- todo3.R: The R code that includes the function ```todo3()``` to obtain the operating characteristics of the proposed design by simulating trials.
  
  ```rscipt
  todo3 (ic12, fdata, ia.post, fa.post1, fa.post2, fa.post3, ia.mean, fa.mean1, fa.mean2, fa.mean3, post1, post2, ntrial, peff, theta0, delta1, ia12, idr.m, wl, m1, nsample)
  ```

- TODO3-m1-n-a1-a2-optimizaiton.R: The R code for optimizing parameters (m1, n, a1, a2).

- TODO3-c1-c2-optimization.R: The R code for optimizing parameters (c1, c2).

- TODO3-simulation.R: The R code for obtaining the operating characteristics of the proposed design by simulating trials.  

# Inputs

- `rseed`: The seed for data generation.

- `iscena`: The scenario indicator.

- `ntrial`: The total number of trials to be simulated.

- `n1`: The number of patients enrolled before the interim analysis.

- `nsample`: The maximum sample size per arm.

- `peff.m`: The specified response rate scenarios.

- `theta0`: The null response rate.

- `delta1`: The non-inferiority margin.

- `sigma1`,`tau`: The hyperparameters for the proposed model.

- `ia12`: The indicator of the optimal combination of parameters (a1, a2).

- `ic12`: The indicator of the optimal combination of parameters (c1, c2).

- `idr.m`: The matrix indicating the incorrect decisions for specified scenarios.

- `wl`: The penalty for inconclusive results.

# Outputs

- `todo2_data_generation()` will return the simulated data and model fitting results of the TODO design for two-dose trials, including:
  
  ```
   (1) the simulated data (data);  
   (2) the posterior probability of better than the null response rate at the interim analysis (ia.post);  
   (3) the posterior probability of better than the null response rate at the final analysis if all arms pass the interim futility monitoring (fa.post1);  
   (4) the posterior probability of better than the null response rate at the final analysis if only dose 2 pass the interim futility monitoring (fa.post2);  
   (5) the posterior probability of better than the null response rate at the final analysis if only dose 1 pass the interim futility monitoring (fa.post3);  
   (6) the posterior probability that dose 1 is non-inferior to dose 2 for the two-dose trials, given that both arms pass the final futility monitoring (post);
   (7) the posterior mean response rates at the interim analysis (ia.mean); 
   (8) the posterior mean response rates at the final analysis if all arms pass the interim futility monitoring (fa.mean1);  
   (9) the posterior mean response rates at the final analysis if only dose 2 passes the interim futility monitoring (fa.mean2);  
   (10) the posterior mean response rates at the final analysis if only dose 1 passes the interim futility monitoring (fa.mean3)
  ```

- `todo2_futility()` and `todo3_futility()` will return the simulated data and model-fitting results of the TODO design for two-dose and three-dose trials, including:
  
  ```
   (1) the true response rates (peff);  
   (2) the average number of responders and average sample size (pts);  
   (3) the (m1, n, a1, a2, c1, c2) used in this trial (parameter);  
   (4) the overall power: the probability of regarding both doses as futile at the final analysis (overallpower)
  ```

- `todo2()` and `todo3()` will return the simulated data, model-fitting results and decision-making of the TODO design for two-dose and three-dose trials, including:
  
  ```
   (1) the true response rates (peff);  
   (2) the average number of responders and average sample size (pts);  
   (3) the (m1, n, a1, a2, c1, c2) used in this trial (parameter);  
   (4) the selection rate of each dose, the probability of having inconclusive results, the probability of making incorrect decisions, the weighted loss, the average sample size for each dose, the probability of making go decision at the final analysis, and the probability of not having an effective dose at the final analysis (present);
   (5) the average bias across doses (bias);
   (6) the average mean square error across doses (MSE)
  ```

# Authors and Reference

- Jingyi Zhang, Heng Zhou, Nolan A. Wages, Zifang Guo, Fang Liu, Thomas Jemielita, Fangrong Yan, and Ruitao Lin
- Zhang, J., Zhou, H., Wages, N. A., Guo, Z., Liu, F., Jemielita, T., Yan, F., and Lin, R. (2024) “TODO: A Triple-Outcome Double-Criterion Optimal Design for Dose Testing-and-Optimization in Multi-Dose Randomized Trials”, Statistics in Medicine, revision submitted to journal.
