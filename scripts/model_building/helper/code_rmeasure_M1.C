double tol = 1e-300; // tolerance value to avoid log(0)
double vtol = 1e-5; // avoid variance to be 0
int u;
double m,v;

for (u = 0; u < U; u++) {
  // daily cases, N=1
  m = C_reported_new[u];
  if (m > 0.0) {
    v = m*tau_cases[u*tau_cases_unit] > 1 ? m*tau_cases[u*tau_cases_unit] : 1;
    daily_cases[u] = nearbyint(rnorm(m, sqrt(v)));
  } else {
    daily_cases[u] = 0.0;
  }

  // daily deaths, N=1
  m = D_new[u];
  if (m > 0.0) {
    v = m*tau_deaths[u*tau_deaths_unit] > 1 ? m*tau_deaths[u*tau_deaths_unit] : 1;
    daily_deaths[u] = nearbyint(rnorm(m, sqrt(v)));
  } else {
    daily_deaths[u] = 0.0;
  }

  // daily sequenced cases all, N=4
  // strain 1
  m = C_1_all_sequenced_new[u];
  if (m > 0.0) {
    v = m*tau_seq_all[u*tau_seq_all_unit] > 1 ? m*tau_seq_all[u*tau_seq_all_unit] : 1;
    C_1_all[u] = nearbyint(rnorm(m, sqrt(v)));
  } else {
    C_1_all[u] = 0.0;
  }
  // strain 2
  m = C_2_all_sequenced_new[u];
  if (m > 0.0) {
    v = m*tau_seq_all[u*tau_seq_all_unit] > 1 ? m*tau_seq_all[u*tau_seq_all_unit] : 1;
    C_2_all[u] = nearbyint(rnorm(m, sqrt(v)));
  } else {
    C_2_all[u] = 0.0;
  }
  // strain 3
  m = C_3_all_sequenced_new[u];
  if (m > 0.0) {
    v = m*tau_seq_all[u*tau_seq_all_unit] > 1 ? m*tau_seq_all[u*tau_seq_all_unit] : 1;
    C_3_all[u] = nearbyint(rnorm(m, sqrt(v)));
  } else {
    C_3_all[u] = 0.0;
  }
  // strain 4
  m = C_4_all_sequenced_new[u];
  if (m > 0.0) {
    v = m*tau_seq_all[u*tau_seq_all_unit] > 1 ? m*tau_seq_all[u*tau_seq_all_unit] : 1;
    C_4_all[u] = nearbyint(rnorm(m, sqrt(v)));
  } else {
    C_4_all[u] = 0.0;
  }

} 
