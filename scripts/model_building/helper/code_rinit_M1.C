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

  E_1_i_direct_in[u] = nearbyint(E_1_i_init[u]*0.86);
  E_2_i_direct_in[u] = nearbyint(E_2_i_init[u]*0.86);
  E_3_i_direct_in[u] = nearbyint(E_3_i_init[u]*0.86);
  E_4_i_direct_in[u] = nearbyint(E_4_i_init[u]*0.86);
  I_1_i_direct_in[u] = nearbyint(I_1_i_init[u]*0.86);
  I_2_i_direct_in[u] = nearbyint(I_2_i_init[u]*0.86);
  I_3_i_direct_in[u] = nearbyint(I_3_i_init[u]*0.86);
  I_4_i_direct_in[u] = nearbyint(I_4_i_init[u]*0.86);

  E_1_i_transit_first[u] = nearbyint(E_1_i_init[u]*0.14);
  E_2_i_transit_first[u] = nearbyint(E_2_i_init[u]*0.14);
  E_3_i_transit_first[u] = nearbyint(E_3_i_init[u]*0.14);
  E_4_i_transit_first[u] = nearbyint(E_4_i_init[u]*0.14);
  I_1_i_transit_first[u] = nearbyint(I_1_i_init[u]*0.14);
  I_2_i_transit_first[u] = nearbyint(I_2_i_init[u]*0.14);
  I_3_i_transit_first[u] = nearbyint(I_3_i_init[u]*0.14);
  I_4_i_transit_first[u] = nearbyint(I_4_i_init[u]*0.14);

  C_1_c_infected_new[u] = 0;
  C_2_c_infected_new[u] = 0;
  C_3_c_infected_new[u] = 0;
  C_4_c_infected_new[u] = 0;
  C_1_c_detected_new[u] = 0;
  C_2_c_detected_new[u] = 0;
  C_3_c_detected_new[u] = 0;
  C_4_c_detected_new[u] = 0;
  C_1_c_sequenced_new[u] = 0;
  C_2_c_sequenced_new[u] = 0;
  C_3_c_sequenced_new[u] = 0;
  C_4_c_sequenced_new[u] = 0;

  C_1_i_infected_new[u] = 0;
  C_2_i_infected_new[u] = 0;
  C_3_i_infected_new[u] = 0;
  C_4_i_infected_new[u] = 0;
  C_1_i_detected_new[u] = 0;
  C_2_i_detected_new[u] = 0;
  C_3_i_detected_new[u] = 0;
  C_4_i_detected_new[u] = 0;
  C_1_i_sequenced_new[u] = 0;
  C_2_i_sequenced_new[u] = 0;
  C_3_i_sequenced_new[u] = 0;
  C_4_i_sequenced_new[u] = 0;

  C_1_all_infected_new[u] = 0;
  C_2_all_infected_new[u] = 0;
  C_3_all_infected_new[u] = 0;
  C_4_all_infected_new[u] = 0;
  C_1_all_detected_new[u] = 0;
  C_2_all_detected_new[u] = 0;
  C_3_all_detected_new[u] = 0;
  C_4_all_detected_new[u] = 0;
  C_1_all_sequenced_new[u] = 0;
  C_2_all_sequenced_new[u] = 0;
  C_3_all_sequenced_new[u] = 0;
  C_4_all_sequenced_new[u] = 0;

  C_reported_new[u] = 0;
  D_new[u] = 0;
}
