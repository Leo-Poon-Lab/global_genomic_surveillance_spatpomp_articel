double tol = 1e-300; // tolerance value to avoid log(0)
double vtol = 1e-5; // tolerance value to avoid numerical instability
double m,v;
int u;
lik = 0;

for (u = 0; u < U; u++) {
  // // daily cases, N=1
  // m = C_reported_new[u] > 0 ? C_reported_new[u] : 0;
  // v = m*tau_cases[u*tau_cases_unit] > 1 ? m*tau_cases[u*tau_cases_unit] : 1;
  // if (daily_cases[u] > tol) {
  //   lik += log(pnorm(daily_cases[u]+0.5,m,sqrt(v),1,0)-pnorm(daily_cases[u]-0.5,m,sqrt(v),1,0) + tol);
  // }

  // // daily deaths, N=1
  // m = D_new[u] > 0 ? D_new[u] : 0;
  // v = m*tau_deaths[u*tau_deaths_unit] > 1 ? m*tau_deaths[u*tau_deaths_unit] : 1;
  // if (daily_deaths[u] > tol) {
  //   lik += log(pnorm(daily_deaths[u]+0.5,m,sqrt(v),1,0)-pnorm(daily_deaths[u]-0.5,m,sqrt(v),1,0) + tol);
  // }

  // daily sequenced cases all, N=4
  // // strain 1
  // m = C_1_all_sequenced_new[u] > 0 ? C_1_all_sequenced_new[u] : 0;
  // v = m*tau_seq_all[u*tau_seq_all_unit] > 1 ? m*tau_seq_all[u*tau_seq_all_unit] : 1;
  // if (C_1_all[u] > tol) {
  //   lik += (log(pnorm(C_1_all[u]+0.5,m,sqrt(v),1,0)-pnorm(C_1_all[u]-0.5,m,sqrt(v),1,0) + tol));
  // }
  // // strain 2
  // m = C_2_all_sequenced_new[u] > 0 ? C_2_all_sequenced_new[u] : 0;
  // v = m*tau_seq_all[u*tau_seq_all_unit] > 1 ? m*tau_seq_all[u*tau_seq_all_unit] : 1;
  // if (C_2_all[u] > tol) {
  //   lik += (log(pnorm(C_2_all[u]+0.5,m,sqrt(v),1,0)-pnorm(C_2_all[u]-0.5,m,sqrt(v),1,0) + tol));
  // }
  // strain 3
  m = C_3_all_sequenced_new[u] > 0 ? C_3_all_sequenced_new[u] : 0;
  v = m*tau_seq_all[u*tau_seq_all_unit] > 1 ? m*tau_seq_all[u*tau_seq_all_unit] : 1;
  if (C_3_all[u] > tol) {
    lik += (log(pnorm(C_3_all[u]+0.5,m,sqrt(v),1,0)-pnorm(C_3_all[u]-0.5,m,sqrt(v),1,0) + tol));
  }
  // // strain 4
  // m = C_4_all_sequenced_new[u] > 0 ? C_4_all_sequenced_new[u] : 0;
  // v = m*tau_seq_all[u*tau_seq_all_unit] > 1 ? m*tau_seq_all[u*tau_seq_all_unit] : 1;
  // if (C_4_all[u] > tol) {
  //   lik += (log(pnorm(C_4_all[u]+0.5,m,sqrt(v),1,0)-pnorm(C_4_all[u]-0.5,m,sqrt(v),1,0) + tol));
  // }

  // // daily sequenced cases community, N=4
  // // strain 1
  // m = C_1_c_sequenced_new[u] > 0 ? C_1_c_sequenced_new[u] : 0;
  // v = m*tau_seq_unit[u*tau_seq_unit_unit] > 1 ? m*tau_seq_unit[u*tau_seq_unit_unit] : 1;
  // if (C_1_c[u] > tol) {
  //   lik += (log(pnorm(C_1_c[u]+0.5,m,sqrt(v),1,0)-pnorm(C_1_c[u]-0.5,m,sqrt(v),1,0) + tol))/N_units_TO_BE_REPLACED;
  // }
  // // strain 2
  // m = C_2_c_sequenced_new[u] > 0 ? C_2_c_sequenced_new[u] : 0;
  // v = m*tau_seq_unit[u*tau_seq_unit_unit] > 1 ? m*tau_seq_unit[u*tau_seq_unit_unit] : 1;
  // if (C_2_c[u] > tol) {
  //   lik += (log(pnorm(C_2_c[u]+0.5,m,sqrt(v),1,0)-pnorm(C_2_c[u]-0.5,m,sqrt(v),1,0) + tol))/N_units_TO_BE_REPLACED;
  // }
  // // strain 3
  // m = C_3_c_sequenced_new[u] > 0 ? C_3_c_sequenced_new[u] : 0;
  // v = m*tau_seq_unit[u*tau_seq_unit_unit] > 1 ? m*tau_seq_unit[u*tau_seq_unit_unit] : 1;
  // if (C_3_c[u] > tol) {
  //   lik += (log(pnorm(C_3_c[u]+0.5,m,sqrt(v),1,0)-pnorm(C_3_c[u]-0.5,m,sqrt(v),1,0) + tol))/N_units_TO_BE_REPLACED;
  // }
  // // strain 4
  // m = C_4_c_sequenced_new[u] > 0 ? C_4_c_sequenced_new[u] : 0;
  // v = m*tau_seq_unit[u*tau_seq_unit_unit] > 1 ? m*tau_seq_unit[u*tau_seq_unit_unit] : 1;
  // if (C_4_c[u] > tol) {
  //   lik += (log(pnorm(C_4_c[u]+0.5,m,sqrt(v),1,0)-pnorm(C_4_c[u]-0.5,m,sqrt(v),1,0) + tol))/N_units_TO_BE_REPLACED;
  // }

  // daily sequenced cases imported, N=116
  

}

if(!give_log) lik = (lik > log(tol)) ? exp(lik) : tol;
