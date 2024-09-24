// routes_E_direct_in
for (v = 0; v < U ; v++) {
  rate[T_E1i0] += v_by_g[day][v][u] * (E_1_c[v]) / population[v];
  rate[T_E2i0] += v_by_g[day][v][u] * (E_2_c[v]) / population[v];
  rate[T_E3i0] += v_by_g[day][v][u] * (E_3_c[v]) / population[v];
  rate[T_E4i0] += v_by_g[day][v][u] * (E_4_c[v]) / population[v];
}
dN[T_E1i0] += rpois(rate[T_E1i0]*dt);
dN[T_E2i0] += rpois(rate[T_E2i0]*dt);
dN[T_E3i0] += rpois(rate[T_E3i0]*dt);
dN[T_E4i0] += rpois(rate[T_E4i0]*dt);

// routes_E_direct_out
// // strain 1
rate[E1i0_P] = rate[T_E1i0] * q_u;
dN[E1i0_P] = rpois(rate[E1i0_P]*dt);
if(dN[E1i0_P] > E_1_i_direct_in[u]){
  dN[E1i0_P] = E_1_i_direct_in[u];
}
rate[E1i0_I1i0] = epsilon_one_u;
reulermultinom(1,E_1_i_direct_in[u]-dN[E1i0_P],&rate[E1i0_I1i0],dt,&dN[E1i0_I1i0]);
// // strain 2
rate[E2i0_P] = rate[T_E2i0] * q_u;
dN[E2i0_P] = rpois(rate[E2i0_P]*dt);
if(dN[E2i0_P] > E_2_i_direct_in[u]){
  dN[E2i0_P] = E_2_i_direct_in[u];
}
rate[E2i0_I2i0] = epsilon_two_u;
reulermultinom(1,E_2_i_direct_in[u]-dN[E2i0_P],&rate[E2i0_I2i0],dt,&dN[E2i0_I2i0]);
// // strain 3
rate[E3i0_P] = rate[T_E3i0] * q_u;
dN[E3i0_P] = rpois(rate[E3i0_P]*dt);
if(dN[E3i0_P] > E_3_i_direct_in[u]){
  dN[E3i0_P] = E_3_i_direct_in[u];
}
rate[E3i0_I3i0] = epsilon_three_u;
reulermultinom(1,E_3_i_direct_in[u]-dN[E3i0_P],&rate[E3i0_I3i0],dt,&dN[E3i0_I3i0]);
// // strain 4
rate[E4i0_P] = rate[T_E4i0] * q_u;
dN[E4i0_P] = rpois(rate[E4i0_P]*dt);
if(dN[E4i0_P] > E_4_i_direct_in[u]){
  dN[E4i0_P] = E_4_i_direct_in[u];
}
rate[E4i0_I4i0] = epsilon_four_u;
reulermultinom(1,E_4_i_direct_in[u]-dN[E4i0_P],&rate[E4i0_I4i0],dt,&dN[E4i0_I4i0]);

// routes_I_direct_in
for (v = 0; v < U ; v++) {
  rate[T_I1i0] += v_by_g[day][v][u] * (I_1_c[v]) / population[v];
  rate[T_I2i0] += v_by_g[day][v][u] * (I_2_c[v]) / population[v];
  rate[T_I3i0] += v_by_g[day][v][u] * (I_3_c[v]) / population[v];
  rate[T_I4i0] += v_by_g[day][v][u] * (I_4_c[v]) / population[v];
}
dN[T_I1i0] += rpois(rate[T_I1i0]*dt);
dN[T_I2i0] += rpois(rate[T_I2i0]*dt);
dN[T_I3i0] += rpois(rate[T_I3i0]*dt);
dN[T_I4i0] += rpois(rate[T_I4i0]*dt);

// routes_I_direct_out
// // strain 1
rate[I1i0_P] = rate[T_I1i0] * q_u;
dN[I1i0_P] = rpois(rate[I1i0_P]*dt);
if(dN[I1i0_P] > I_1_i_direct_in[u]){
  dN[I1i0_P] = I_1_i_direct_in[u];
}
rate[I1i0_D] = infection_fatality[u] * gamma_u;
rate[I1i0_R] = (1-infection_fatality[u]) * gamma_u;
reulermultinom(2,I_1_i_direct_in[u]-dN[I1i0_P],&rate[I1i0_D],dt,&dN[I1i0_D]);
// // strain 2
rate[I2i0_P] = rate[T_I2i0] * q_u;
dN[I2i0_P] = rpois(rate[I2i0_P]*dt);
if(dN[I2i0_P] > I_2_i_direct_in[u]){
  dN[I2i0_P] = I_2_i_direct_in[u];
}
rate[I2i0_D] = infection_fatality[u] * gamma_u;
rate[I2i0_R] = (1-infection_fatality[u]) * gamma_u;
reulermultinom(2,I_2_i_direct_in[u]-dN[I2i0_P],&rate[I2i0_D],dt,&dN[I2i0_D]);
// // strain 3
rate[I3i0_P] = rate[T_I3i0] * q_u;
dN[I3i0_P] = rpois(rate[I3i0_P]*dt);
if(dN[I3i0_P] > I_3_i_direct_in[u]){
  dN[I3i0_P] = I_3_i_direct_in[u];
}
rate[I3i0_D] = infection_fatality[u] * gamma_u;
rate[I3i0_R] = (1-infection_fatality[u]) * gamma_u;
reulermultinom(2,I_3_i_direct_in[u]-dN[I3i0_P],&rate[I3i0_D],dt,&dN[I3i0_D]);
// // strain 4
rate[I4i0_P] = rate[T_I4i0] * q_u;
dN[I4i0_P] = rpois(rate[I4i0_P]*dt);
if(dN[I4i0_P] > I_4_i_direct_in[u]){
  dN[I4i0_P] = I_4_i_direct_in[u];
}
rate[I4i0_D] = infection_fatality[u] * gamma_u;
rate[I4i0_R] = (1-infection_fatality[u]) * gamma_u;
reulermultinom(2,I_4_i_direct_in[u]-dN[I4i0_P],&rate[I4i0_D],dt,&dN[I4i0_D]);

// assigning dNs
// // strain 1
E_1_i_direct_in[u] += nearbyint(dN[T_E1i0] - dN[E1i0_P] - dN[E1i0_I1i0]);
E_1_i_direct_in[u] = E_1_i_direct_in[u] < 0 ? 0 : E_1_i_direct_in[u];
I_1_i_direct_in[u] += nearbyint(dN[E1i0_I1i0] + dN[T_I1i0] - dN[I1i0_P] - dN[I1i0_D] - dN[I1i0_R]);
I_1_i_direct_in[u] = I_1_i_direct_in[u] < 0 ? 0 : I_1_i_direct_in[u];
for (int i = 0; i < U; i++) { // this part is for P_E1i1 and P_I1i1
  // Assigning E transit travelers to the next level
  E_1_i_transit_first[i] += nearbyint(dN[E1i0_P] * v_by_g[day][u][i] /total_out_flow_t[u]);
  // Assigning I transit travelers to the next level
  I_1_i_transit_first[i] += nearbyint(dN[I1i0_P] * v_by_g[day][u][i] /total_out_flow_t[u]);
}
// // strain 2
E_2_i_direct_in[u] += nearbyint(dN[T_E2i0] - dN[E2i0_P] - dN[E2i0_I2i0]);
E_2_i_direct_in[u] = E_2_i_direct_in[u] < 0 ? 0 : E_2_i_direct_in[u];
I_2_i_direct_in[u] += nearbyint(dN[E2i0_I2i0] + dN[T_I2i0] - dN[I2i0_P] - dN[I2i0_D] - dN[I2i0_R]);
I_2_i_direct_in[u] = I_2_i_direct_in[u] < 0 ? 0 : I_2_i_direct_in[u];
for (int i = 0; i < U; i++) { // this part is for P_E2i1 and P_I2i1
  // Assigning E transit travelers to the next level
  E_2_i_transit_first[i] += nearbyint(dN[E2i0_P] * v_by_g[day][u][i] /total_out_flow_t[u]);
  // Assigning I transit travelers to the next level
  I_2_i_transit_first[i] += nearbyint(dN[I2i0_P] * v_by_g[day][u][i] /total_out_flow_t[u]);
}
// // strain 3
E_3_i_direct_in[u] += nearbyint(dN[T_E3i0] - dN[E3i0_P] - dN[E3i0_I3i0]);
E_3_i_direct_in[u] = E_3_i_direct_in[u] < 0 ? 0 : E_3_i_direct_in[u];
I_3_i_direct_in[u] += nearbyint(dN[E3i0_I3i0] + dN[T_I3i0] - dN[I3i0_P] - dN[I3i0_D] - dN[I3i0_R]);
I_3_i_direct_in[u] = I_3_i_direct_in[u] < 0 ? 0 : I_3_i_direct_in[u];
for (int i = 0; i < U; i++) { // this part is for P_E3i1 and P_I3i1
  // Assigning E transit travelers to the next level
  E_3_i_transit_first[i] += nearbyint(dN[E3i0_P] * v_by_g[day][u][i] /total_out_flow_t[u]);
  // Assigning I transit travelers to the next level
  I_3_i_transit_first[i] += nearbyint(dN[I3i0_P] * v_by_g[day][u][i] /total_out_flow_t[u]);
}
// // strain 4
E_4_i_direct_in[u] += nearbyint(dN[T_E4i0] - dN[E4i0_P] - dN[E4i0_I4i0]);
E_4_i_direct_in[u] = E_4_i_direct_in[u] < 0 ? 0 : E_4_i_direct_in[u];
I_4_i_direct_in[u] += nearbyint(dN[E4i0_I4i0] + dN[T_I4i0] - dN[I4i0_P] - dN[I4i0_D] - dN[I4i0_R]);
I_4_i_direct_in[u] = I_4_i_direct_in[u] < 0 ? 0 : I_4_i_direct_in[u];
for (int i = 0; i < U; i++) { // this part is for P_E4i1 and P_I4i1
  // Assigning E transit travelers to the next level
  E_4_i_transit_first[i] += nearbyint(dN[E4i0_P] * v_by_g[day][u][i] /total_out_flow_t[u]);
  // Assigning I transit travelers to the next level
  I_4_i_transit_first[i] += nearbyint(dN[I4i0_P] * v_by_g[day][u][i] /total_out_flow_t[u]);
}
