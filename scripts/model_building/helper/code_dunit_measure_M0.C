double tol = 1e-300; // tolerance value to avoid log(0)
double vtol = 1e-5; // avoid variance to be 0
double m,v;
lik = 0;

// daily cases, N=1
m = C_reported_new > 0 ? C_reported_new : 0;
v = m*tau_cases[u*tau_cases_unit] > 1 ? m*tau_cases[u*tau_cases_unit] : 1;
if (daily_cases > tol) {lik += log(pnorm(daily_cases+0.5, m, sqrt(v), 1, 0) - pnorm(daily_cases-0.5, m, sqrt(v), 1, 0) + tol);}

// daily deaths, N=1
m = D_new > 0 ? D_new : 0;
v = m*tau_deaths[u*tau_deaths_unit] > 1 ? m*tau_deaths[u*tau_deaths_unit] : 1;
if (daily_deaths > tol) {lik += log(pnorm(daily_deaths+0.5, m, sqrt(v), 1, 0) - pnorm(daily_deaths-0.5, m, sqrt(v), 1, 0) + tol);}

// daily sequenced cases all, N=4
// strain 1
m = C_1_all_sequenced_new > 0 ? C_1_all_sequenced_new : 0;
v = m*tau_seq_all[u*tau_seq_all_unit] > 1 ? m*tau_seq_all[u*tau_seq_all_unit] : 1;
if(C_1_all > tol) {lik += log(pnorm(C_1_all+0.5, m, sqrt(v), 1, 0) - pnorm(C_1_all-0.5, m, sqrt(v), 1, 0) + tol);}
// strain 2
m = C_2_all_sequenced_new > 0 ? C_2_all_sequenced_new : 0;
v = m*tau_seq_all[u*tau_seq_all_unit] > 1 ? m*tau_seq_all[u*tau_seq_all_unit] : 1;
if(C_2_all > tol) {lik += log(pnorm(C_2_all+0.5, m, sqrt(v), 1, 0) - pnorm(C_2_all-0.5, m, sqrt(v), 1, 0) + tol)*1.5;}
// strain 3
m = C_3_all_sequenced_new > 0 ? C_3_all_sequenced_new : 0;
v = m*tau_seq_all[u*tau_seq_all_unit] > 1 ? m*tau_seq_all[u*tau_seq_all_unit] : 1;
if(C_3_all > tol) {lik += log(pnorm(C_3_all+0.5, m, sqrt(v), 1, 0) - pnorm(C_3_all-0.5, m, sqrt(v), 1, 0) + tol)*2;}
// strain 4
m = C_4_all_sequenced_new > 0 ? C_4_all_sequenced_new : 0;
v = m*tau_seq_all[u*tau_seq_all_unit] > 1 ? m*tau_seq_all[u*tau_seq_all_unit] : 1;
if(C_4_all > tol) {lik += log(pnorm(C_4_all+0.5, m, sqrt(v), 1, 0) - pnorm(C_4_all-0.5, m, sqrt(v), 1, 0) + tol);} 

// daily sequenced cases community, N=4
// strain 1
m = C_1_c_sequenced_new > 0 ? C_1_c_sequenced_new : 0;
v = m*tau_seq_unit[u*tau_seq_unit_unit] > 1 ? m*tau_seq_unit[u*tau_seq_unit_unit] : 1;
if(C_1_c > tol) {lik += log(pnorm(C_1_c+0.5, m, sqrt(v), 1, 0) - pnorm(C_1_c - 0.5, m, sqrt(v), 1, 0) + tol)/N_units_TO_BE_REPLACED;}
// strain 2
m = C_2_c_sequenced_new > 0 ? C_2_c_sequenced_new : 0;
v = m*tau_seq_unit[u*tau_seq_unit_unit] > 1 ? m*tau_seq_unit[u*tau_seq_unit_unit] : 1;
if(C_2_c > tol) {lik += log(pnorm(C_2_c+0.5, m, sqrt(v), 1, 0) - pnorm(C_2_c - 0.5, m, sqrt(v), 1, 0) + tol)/N_units_TO_BE_REPLACED;}
// strain 3
m = C_3_c_sequenced_new > 0 ? C_3_c_sequenced_new : 0;
v = m*tau_seq_unit[u*tau_seq_unit_unit] > 1 ? m*tau_seq_unit[u*tau_seq_unit_unit] : 1;
if(C_3_c > tol) {lik += log(pnorm(C_3_c+0.5, m, sqrt(v), 1, 0) - pnorm(C_3_c - 0.5, m, sqrt(v), 1, 0) + tol)/N_units_TO_BE_REPLACED;}
// strain 4
m = C_4_c_sequenced_new > 0 ? C_4_c_sequenced_new : 0;
v = m*tau_seq_unit[u*tau_seq_unit_unit] > 1 ? m*tau_seq_unit[u*tau_seq_unit_unit] : 1;
if(C_4_c > tol) {lik += log(pnorm(C_4_c+0.5, m, sqrt(v), 1, 0) - pnorm(C_4_c - 0.5, m, sqrt(v), 1, 0) + tol)/N_units_TO_BE_REPLACED;}

// daily sequenced cases imported, N=116
code_dunit_measure_C_k_i_o_TO_BE_REPLACED;

if(!give_log) lik = (lik > log(tol)) ? exp(lik) : tol;
