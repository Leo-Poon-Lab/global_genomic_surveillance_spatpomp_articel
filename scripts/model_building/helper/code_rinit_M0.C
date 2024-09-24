int u;
for(u = 0; u < U; u++) {
  S[u] = nearbyint(S_init[u]);
  V[u] = nearbyint(V_init[u]);
  R[u] = nearbyint(R_init[u]);
  D[u] = 0;

  E_1_c[u] = nearbyint(E_1_c_init[u]);
  E_2_c[u] = nearbyint(E_2_c_init[u]);
  E_3_c[u] = nearbyint(E_3_c_init[u]);
  E_4_c[u] = nearbyint(E_4_c_init[u]);
  I_1_c[u] = nearbyint(I_1_c_init[u]);
  I_2_c[u] = nearbyint(I_2_c_init[u]);
  I_3_c[u] = nearbyint(I_3_c_init[u]);
  I_4_c[u] = nearbyint(I_4_c_init[u]);

  E_1_i[u] = nearbyint(E_1_i_init[u]);
  E_2_i[u] = nearbyint(E_2_i_init[u]);
  E_3_i[u] = nearbyint(E_3_i_init[u]);
  E_4_i[u] = nearbyint(E_4_i_init[u]);
  I_1_i[u] = nearbyint(I_1_i_init[u]);
  I_2_i[u] = nearbyint(I_2_i_init[u]);
  I_3_i[u] = nearbyint(I_3_i_init[u]);
  I_4_i[u] = nearbyint(I_4_i_init[u]);

  E_by_direct_in_TO_BE_REPLACED;
  I_by_direct_in_TO_BE_REPLACED;
  E_by_transit_first_TO_BE_REPLACED;
  I_by_transit_first_TO_BE_REPLACED;

  C_1_c_infected_new[u] = 0;
  C_2_c_infected_new[u] = 0;
  C_3_c_infected_new[u] = 0;
  C_4_c_infected_new[u] = 0;
  C_1_c_detected_new[u] = 0;
  C_2_c_detected_new[u] = 0;
  C_3_c_detected_new[u] = 0;
  C_4_c_detected_new[u] = 0;
  C_1_c_detected_new[u] = 0;
  C_2_c_detected_new[u] = 0;
  C_3_c_detected_new[u] = 0;
  C_4_c_detected_new[u] = 0;

  New_cases_infected_TO_BE_REPLACED;
  New_cases_detected_TO_BE_REPLACED;
  New_cases_sequenced_TO_BE_REPLACED;

  C_1_all_sequenced_new[u] = 0;
  C_2_all_sequenced_new[u] = 0;
  C_3_all_sequenced_new[u] = 0;
  C_4_all_sequenced_new[u] = 0;

  C_reported_new[u] = 0;
  D_new[u] = 0;
}
