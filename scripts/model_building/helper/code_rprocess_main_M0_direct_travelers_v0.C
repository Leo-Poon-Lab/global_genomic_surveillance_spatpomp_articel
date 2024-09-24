if(v==0){ // v==0
  // routes_E_direct_in
  // // strain 1
  rate[T_E1i0_origin_0] = v_by_g[day][v][u] * (E_1_c[v]) / population[v];
  dN[T_E1i0_origin_0] = rpois(rate[T_E1i0_origin_0]*dt);
  // // strain 2
  rate[T_E2i0_origin_0] = v_by_g[day][v][u] * (E_2_c[v]) / population[v];
  dN[T_E2i0_origin_0] = rpois(rate[T_E2i0_origin_0]*dt);
  // // strain 3
  rate[T_E3i0_origin_0] = v_by_g[day][v][u] * (E_3_c[v]) / population[v];
  dN[T_E3i0_origin_0] = rpois(rate[T_E3i0_origin_0]*dt);
  // // strain 4
  rate[T_E4i0_origin_0] = v_by_g[day][v][u] * (E_4_c[v]) / population[v];
  dN[T_E4i0_origin_0] = rpois(rate[T_E4i0_origin_0]*dt);

  // routes_E_direct_out
  // // strain 1
  rate[E1i0_P_origin_0] = rate[T_E1i0_origin_0] * q_u;
  dN[E1i0_P_origin_0] = rpois(rate[E1i0_P_origin_0]*dt);
  if(dN[E1i0_P_origin_0] > E_1_i_direct_in_origin_0_[u]){
    dN[E1i0_P_origin_0] = E_1_i_direct_in_origin_0_[u];
  }
  rate[E1i0_I1i0_origin_0] = epsilon_one_u;
  reulermultinom(1,E_1_i_direct_in_origin_0_[u]-dN[E1i0_P_origin_0],&rate[E1i0_I1i0_origin_0],dt,&dN[E1i0_I1i0_origin_0]);
  // // strain 2
  rate[E2i0_P_origin_0] = rate[T_E2i0_origin_0] * q_u;
  dN[E2i0_P_origin_0] = rpois(rate[E2i0_P_origin_0]*dt);
  if(dN[E2i0_P_origin_0] > E_2_i_direct_in_origin_0_[u]){
    dN[E2i0_P_origin_0] = E_2_i_direct_in_origin_0_[u];
  }
  rate[E2i0_I2i0_origin_0] = epsilon_two_u;
  reulermultinom(1,E_2_i_direct_in_origin_0_[u]-dN[E2i0_P_origin_0],&rate[E2i0_I2i0_origin_0],dt,&dN[E2i0_I2i0_origin_0]);
  // // strain 3
  rate[E3i0_P_origin_0] = rate[T_E3i0_origin_0] * q_u;
  dN[E3i0_P_origin_0] = rpois(rate[E3i0_P_origin_0]*dt);
  if(dN[E3i0_P_origin_0] > E_3_i_direct_in_origin_0_[u]){
    dN[E3i0_P_origin_0] = E_3_i_direct_in_origin_0_[u];
  }
  rate[E3i0_I3i0_origin_0] = epsilon_three_u;
  reulermultinom(1,E_3_i_direct_in_origin_0_[u]-dN[E3i0_P_origin_0],&rate[E3i0_I3i0_origin_0],dt,&dN[E3i0_I3i0_origin_0]);
  // // strain 4
  rate[E4i0_P_origin_0] = rate[T_E4i0_origin_0] * q_u;
  dN[E4i0_P_origin_0] = rpois(rate[E4i0_P_origin_0]*dt);
  if(dN[E4i0_P_origin_0] > E_4_i_direct_in_origin_0_[u]){
    dN[E4i0_P_origin_0] = E_4_i_direct_in_origin_0_[u];
  }
  rate[E4i0_I4i0_origin_0] = epsilon_four_u;
  reulermultinom(1,E_4_i_direct_in_origin_0_[u]-dN[E4i0_P_origin_0],&rate[E4i0_I4i0_origin_0],dt,&dN[E4i0_I4i0_origin_0]);

  // routes_I_direct_in
  // // strain 1
  rate[T_I1i0_origin_0] = v_by_g[day][v][u] * (I_1_c[v]) / population[v];
  dN[T_I1i0_origin_0] = rpois(rate[T_I1i0_origin_0]*dt);
  // // strain 2
  rate[T_I2i0_origin_0] = v_by_g[day][v][u] * (I_2_c[v]) / population[v];
  dN[T_I2i0_origin_0] = rpois(rate[T_I2i0_origin_0]*dt);
  // // strain 3
  rate[T_I3i0_origin_0] = v_by_g[day][v][u] * (I_3_c[v]) / population[v];
  dN[T_I3i0_origin_0] = rpois(rate[T_I3i0_origin_0]*dt);
  // // strain 4
  rate[T_I4i0_origin_0] = v_by_g[day][v][u] * (I_4_c[v]) / population[v];
  dN[T_I4i0_origin_0] = rpois(rate[T_I4i0_origin_0]*dt);

  // routes_I_direct_out
  // // strain 1
  rate[I1i0_P_origin_0] = rate[T_I1i0_origin_0] * q_u;
  dN[I1i0_P_origin_0] = rpois(rate[I1i0_P_origin_0]*dt);
  if(dN[I1i0_P_origin_0] > I_1_i_direct_in_origin_0_[u]){
    dN[I1i0_P_origin_0] = I_1_i_direct_in_origin_0_[u];
  }
  rate[I1i0_D_origin_0] = infection_fatality[u] * gamma_u;
  rate[I1i0_R_origin_0] = (1-infection_fatality[u]) * gamma_u;
  reulermultinom(2,I_1_i_direct_in_origin_0_[u]-dN[I1i0_P_origin_0],&rate[I1i0_D_origin_0],dt,&dN[I1i0_D_origin_0]);
  // // strain 2
  rate[I2i0_P_origin_0] = rate[T_I2i0_origin_0] * q_u;
  dN[I2i0_P_origin_0] = rpois(rate[I2i0_P_origin_0]*dt);
  if(dN[I2i0_P_origin_0] > I_2_i_direct_in_origin_0_[u]){
    dN[I2i0_P_origin_0] = I_2_i_direct_in_origin_0_[u];
  }
  rate[I2i0_D_origin_0] = infection_fatality[u] * gamma_u;
  rate[I2i0_R_origin_0] = (1-infection_fatality[u]) * gamma_u;
  reulermultinom(2,I_2_i_direct_in_origin_0_[u]-dN[I2i0_P_origin_0],&rate[I2i0_D_origin_0],dt,&dN[I2i0_D_origin_0]);
  // // strain 3
  rate[I3i0_P_origin_0] = rate[T_I3i0_origin_0] * q_u;
  dN[I3i0_P_origin_0] = rpois(rate[I3i0_P_origin_0]*dt);
  if(dN[I3i0_P_origin_0] > I_3_i_direct_in_origin_0_[u]){
    dN[I3i0_P_origin_0] = I_3_i_direct_in_origin_0_[u];
  }
  rate[I3i0_D_origin_0] = infection_fatality[u] * gamma_u;
  rate[I3i0_R_origin_0] = (1-infection_fatality[u]) * gamma_u;
  reulermultinom(2,I_3_i_direct_in_origin_0_[u]-dN[I3i0_P_origin_0],&rate[I3i0_D_origin_0],dt,&dN[I3i0_D_origin_0]);
  // // strain 4
  rate[I4i0_P_origin_0] = rate[T_I4i0_origin_0] * q_u;
  dN[I4i0_P_origin_0] = rpois(rate[I4i0_P_origin_0]*dt);
  if(dN[I4i0_P_origin_0] > I_4_i_direct_in_origin_0_[u]){
    dN[I4i0_P_origin_0] = I_4_i_direct_in_origin_0_[u];
  }
  rate[I4i0_D_origin_0] = infection_fatality[u] * gamma_u;
  rate[I4i0_R_origin_0] = (1-infection_fatality[u]) * gamma_u;
  reulermultinom(2,I_4_i_direct_in_origin_0_[u]-dN[I4i0_P_origin_0],&rate[I4i0_D_origin_0],dt,&dN[I4i0_D_origin_0]);

  // assigning dNs
  // // strain 1
  E_1_i_direct_in_origin_0_[u] += nearbyint(dN[T_E1i0_origin_0] - dN[E1i0_P_origin_0] - dN[E1i0_I1i0_origin_0]);
  E_1_i_direct_in_origin_0_[u] = E_1_i_direct_in_origin_0_[u] < 0 ? 0 : E_1_i_direct_in_origin_0_[u];
  I_1_i_direct_in_origin_0_[u] += nearbyint(dN[E1i0_I1i0_origin_0] + dN[T_I1i0_origin_0] - dN[I1i0_P_origin_0] - dN[I1i0_D_origin_0] - dN[I1i0_R_origin_0]);
  I_1_i_direct_in_origin_0_[u] = I_1_i_direct_in_origin_0_[u] < 0 ? 0 : I_1_i_direct_in_origin_0_[u];
  for (int i = 0; i < U; i++) { // this part is for P_E1i1 and P_I1i1
    // Assigning E transit travelers to the next level
    E_1_i_transit_first_origin_0_[i] += nearbyint(dN[E1i0_P_origin_0] * v_by_g[day][u][i] /total_out_flow_t[u]);
    // Assigning I transit travelers to the next level
    I_1_i_transit_first_origin_0_[i] += nearbyint(dN[I1i0_P_origin_0] * v_by_g[day][u][i] /total_out_flow_t[u]);
  }
  // // strain 2
  E_2_i_direct_in_origin_0_[u] += nearbyint(dN[T_E2i0_origin_0] - dN[E2i0_P_origin_0] - dN[E2i0_I2i0_origin_0]);
  E_2_i_direct_in_origin_0_[u] = E_2_i_direct_in_origin_0_[u] < 0 ? 0 : E_2_i_direct_in_origin_0_[u];
  I_2_i_direct_in_origin_0_[u] += nearbyint(dN[E2i0_I2i0_origin_0] + dN[T_I2i0_origin_0] - dN[I2i0_P_origin_0] - dN[I2i0_D_origin_0] - dN[I2i0_R_origin_0]);
  I_2_i_direct_in_origin_0_[u] = I_2_i_direct_in_origin_0_[u] < 0 ? 0 : I_2_i_direct_in_origin_0_[u];
  for (int i = 0; i < U; i++) { // this part is for P_E2i1 and P_I2i1
    // Assigning E transit travelers to the next level
    E_2_i_transit_first_origin_0_[i] += nearbyint(dN[E2i0_P_origin_0] * v_by_g[day][u][i] /total_out_flow_t[u]);
    // Assigning I transit travelers to the next level
    I_2_i_transit_first_origin_0_[i] += nearbyint(dN[I2i0_P_origin_0] * v_by_g[day][u][i] /total_out_flow_t[u]);
  }
  // // strain 3
  E_3_i_direct_in_origin_0_[u] += nearbyint(dN[T_E3i0_origin_0] - dN[E3i0_P_origin_0] - dN[E3i0_I3i0_origin_0]);
  E_3_i_direct_in_origin_0_[u] = E_3_i_direct_in_origin_0_[u] < 0 ? 0 : E_3_i_direct_in_origin_0_[u];
  I_3_i_direct_in_origin_0_[u] += nearbyint(dN[E3i0_I3i0_origin_0] + dN[T_I3i0_origin_0] - dN[I3i0_P_origin_0] - dN[I3i0_D_origin_0] - dN[I3i0_R_origin_0]);
  I_3_i_direct_in_origin_0_[u] = I_3_i_direct_in_origin_0_[u] < 0 ? 0 : I_3_i_direct_in_origin_0_[u];
  for (int i = 0; i < U; i++) { // this part is for P_E3i1 and P_I3i1
    // Assigning E transit travelers to the next level
    E_3_i_transit_first_origin_0_[i] += nearbyint(dN[E3i0_P_origin_0] * v_by_g[day][u][i] /total_out_flow_t[u]);
    // Assigning I transit travelers to the next level
    I_3_i_transit_first_origin_0_[i] += nearbyint(dN[I3i0_P_origin_0] * v_by_g[day][u][i] /total_out_flow_t[u]);
  }
  // // strain 4
  E_4_i_direct_in_origin_0_[u] += nearbyint(dN[T_E4i0_origin_0] - dN[E4i0_P_origin_0] - dN[E4i0_I4i0_origin_0]);
  E_4_i_direct_in_origin_0_[u] = E_4_i_direct_in_origin_0_[u] < 0 ? 0 : E_4_i_direct_in_origin_0_[u];
  I_4_i_direct_in_origin_0_[u] += nearbyint(dN[E4i0_I4i0_origin_0] + dN[T_I4i0_origin_0] - dN[I4i0_P_origin_0] - dN[I4i0_D_origin_0] - dN[I4i0_R_origin_0]);
  I_4_i_direct_in_origin_0_[u] = I_4_i_direct_in_origin_0_[u] < 0 ? 0 : I_4_i_direct_in_origin_0_[u];
  for (int i = 0; i < U; i++) { // this part is for P_E4i1 and P_I4i1
    // Assigning E transit travelers to the next level
    E_4_i_transit_first_origin_0_[i] += nearbyint(dN[E4i0_P_origin_0] * v_by_g[day][u][i] /total_out_flow_t[u]);
    // Assigning I transit travelers to the next level
    I_4_i_transit_first_origin_0_[i] += nearbyint(dN[I4i0_P_origin_0] * v_by_g[day][u][i] /total_out_flow_t[u]);
  }
}
