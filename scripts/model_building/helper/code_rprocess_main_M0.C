double rate[num_of_transition_TO_BE_REPLACED], dN[num_of_transition_TO_BE_REPLACED];
double dw=0;
double iota_u=0, I_tilde1=0, I_tilde2=0, I_tilde3=0, I_tilde4=0;
double eta_one_u=0, epsilon_one_u=0, epsilon_two_u=0, epsilon_three_u=0, q_u=0; // for fixed parameters
double total_out_flow_t[U];
int day = (int)((t-t_start_d)*365);
int u=0, v=0;

initialize_transformed_params_TO_BE_REPLACED;
initialize_transition_order_TO_BE_REPLACED;

// Manually set upper limit for introduction dates
int day_BA_one_int = (int)(3+42*day_BA_one[0]); // range between 3 and 45, i.e. between Sep-18 to Oct-25
int day_BA_two_int = (int)(10+65*day_BA_two[0]); // range between 10 and 75, i.e. between Sep-25 to Nov-25

for (u = 0 ; u < U ; u++) {
  total_out_flow_t[u]=0;
  for (v = 0; v < U ; v++) {
    total_out_flow_t[u] += v_by_g[day][u][v];
  }
  if(total_out_flow_t[u]==0){
    total_out_flow_t[u] = 1; // avoid dividing by zero
  }
}

for (u = 0 ; u < U ; u++) {
  // Reset rate and dN to be zero
  reset_rate_and_dN_TO_BE_REPLACED;

  // Manually transform parameters, from 0-1 to appropriate limits
  code_transform_parameters_TO_BE_REPLACED; // these are for estimating parameters
  eta_one_u = eta_one[u*eta_one_unit]; // this is a fixed parameter
  epsilon_one_u = epsilon_one[u*epsilon_one_unit]; // this is a fixed parameter
  epsilon_two_u = epsilon_two[u*epsilon_two_unit]; // this is a fixed parameter
  epsilon_three_u = epsilon_three[u*epsilon_three_unit]; // this is a fixed parameter
  q_u = q[u*q_unit]; // this is a fixed parameter

  // white noise (extrademographic stochasticity)
  dw = rgammawn(sigma_u,dt);
  
  // Update the aggregated E_k_i and I_k_i for each unit
  code_update_E_k_i_TO_BE_REPLACED;
  code_update_I_k_i_TO_BE_REPLACED;

  // Introduce Omicron BA.1 and BA.2 in South Africa
  // Needs to be changed when runnning new model with different number of spatial units
  if(day==day_BA_one_int){ // BA.1 introduction
    if(u==ZA_code){ // zero index number for South Africa, will be replaced later in R code
      E_2_c[u] += nearbyint(N_BA_one_introduced); // later replaced by R code
      I_2_c[u] += nearbyint(N_BA_one_introduced);
      //E_2_i[u] += nearbyint(N_BA_one_introduced/10);
      //I_2_i[u] += nearbyint(N_BA_one_introduced/10);

      // E_1_c[u] -= nearbyint(N_BA_one_introduced/2);
      // I_1_c[u] -= nearbyint(N_BA_one_introduced/2);
      // E_4_c[u] -= nearbyint(N_BA_one_introduced/2);
      // I_4_c[u] -= nearbyint(N_BA_one_introduced/2);
    }
  }
  if(day==day_BA_two_int){ // BA.2 
    if(u==ZA_code){ // zero index number for South Africa, will be replaced later in R code
      E_3_c[u] += nearbyint(N_BA_one_introduced);
      I_3_c[u] += nearbyint(N_BA_one_introduced);
      //E_3_i[u] += nearbyint(N_BA_one_introduced/10);
      //I_3_i[u] += nearbyint(N_BA_one_introduced/10);

      // E_1_c[u] -= nearbyint(N_BA_one_introduced/3);
      // I_1_c[u] -= nearbyint(N_BA_one_introduced/3);
      // E_4_c[u] -= nearbyint(N_BA_one_introduced/3);
      // I_4_c[u] -= nearbyint(N_BA_one_introduced/3);
      // E_2_c[u] -= nearbyint(N_BA_one_introduced/3);
      // I_2_c[u] -= nearbyint(N_BA_one_introduced/3);
    }
  }

  // Community related transitions, illustrated as below:
  //     > V        > T      > T > D
  //    /   \      /        /   /
  //   /     \>   /        /   /
  //  S  -->  E_k_c  -->  I_k_c  -->  R

  // S to V
  rate[S_V] = vaxx_rate[u];

  // S to E
  if(travel_control_level[u]==1){
    iota_u = iota_one_u;
  } else if(travel_control_level[u]==2){
    iota_u = iota_two_u;
  } else if(travel_control_level[u]==3){
    iota_u = iota_three_u;
  } else if(travel_control_level[u]==4){
    iota_u = iota_four_u;
  } else {
    iota_u = 1;
  }

  I_tilde1 = I_1_i[u]*iota_u + I_1_c[u];
  I_tilde2 = I_2_i[u]*iota_u + I_2_c[u];
  I_tilde3 = I_3_i[u]*iota_u + I_3_c[u];
  I_tilde4 = I_4_i[u]*iota_u + I_4_c[u];
  
  // stochastic force of infection
  rate[S_E1c] = eta_one_u*beta_naught_u/population[u]*I_tilde1*dw/dt;
  rate[S_E2c] = eta_two_u*beta_naught_u/population[u]*I_tilde2*dw/dt;
  rate[S_E3c] = eta_three_u*beta_naught_u/population[u]*I_tilde3*dw/dt;
  rate[S_E4c] = eta_four_u*beta_naught_u/population[u]*I_tilde4*dw/dt;

  reulermultinom(5,S[u],&rate[S_V],dt,&dN[S_V]);

  // V to E
  rate[V_E1c] = lambda_one_u*eta_one_u*beta_naught_u/population[u]*I_tilde1*dw/dt;
  rate[V_E2c] = lambda_two_u*eta_two_u*beta_naught_u/population[u]*I_tilde2*dw/dt;
  rate[V_E3c] = lambda_three_u*eta_three_u*beta_naught_u/population[u]*I_tilde3*dw/dt;
  rate[V_E4c] = lambda_four_u*eta_four_u*beta_naught_u/population[u]*I_tilde4*dw/dt;

  reulermultinom(4,V[u],&rate[V_E1c],dt,&dN[V_E1c]);

  // Community E to T, and community I to T
  rate[E1c_T] = total_out_flow_t[u] * E_1_c[u] / population[u];
  dN[E1c_T] = rpois(rate[E1c_T]*dt);
  if(dN[E1c_T] > E_1_c[u]){
    dN[E1c_T] = E_1_c[u];
  }
  rate[E2c_T] = total_out_flow_t[u] * E_2_c[u] / population[u];
  dN[E2c_T] = rpois(rate[E2c_T]*dt);
  if(dN[E2c_T] > E_2_c[u]){
    dN[E2c_T] = E_2_c[u];
  }
  rate[E3c_T] = total_out_flow_t[u] * E_3_c[u] / population[u];
  dN[E3c_T] = rpois(rate[E3c_T]*dt);
  if(dN[E3c_T] > E_3_c[u]){
    dN[E3c_T] = E_3_c[u];
  }
  rate[E4c_T] = total_out_flow_t[u] * E_4_c[u] / population[u];
  dN[E4c_T] = rpois(rate[E4c_T]*dt);
  if(dN[E4c_T] > E_4_c[u]){
    dN[E4c_T] = E_4_c[u];
  }

  rate[I1c_T] = total_out_flow_t[u] * I_1_c[u] / population[u];
  dN[I1c_T] = rpois(rate[I1c_T]*dt);
  if(dN[I1c_T] > I_1_c[u]){
    dN[I1c_T] = I_1_c[u];
  }
  rate[I2c_T] = total_out_flow_t[u] * I_2_c[u] / population[u];
  dN[I2c_T] = rpois(rate[I2c_T]*dt);
  if(dN[I2c_T] > I_2_c[u]){
    dN[I2c_T] = I_2_c[u];
  }
  rate[I3c_T] = total_out_flow_t[u] * I_3_c[u] / population[u];
  dN[I3c_T] = rpois(rate[I3c_T]*dt);
  if(dN[I3c_T] > I_3_c[u]){
    dN[I3c_T] = I_3_c[u];
  }
  rate[I4c_T] = total_out_flow_t[u] * I_4_c[u] / population[u];
  dN[I4c_T] = rpois(rate[I4c_T]*dt);
  if(dN[I4c_T] > I_4_c[u]){
    dN[I4c_T] = I_4_c[u];
  }

  // Community E to I
  rate[E1c_I1c] = epsilon_one_u;
  rate[E2c_I2c] = epsilon_two_u;
  rate[E3c_I3c] = epsilon_three_u;
  rate[E4c_I4c] = epsilon_four_u;

  reulermultinom(1,E_1_c[u]-dN[E1c_T],&rate[E1c_I1c],dt,&dN[E1c_I1c]);
  reulermultinom(1,E_2_c[u]-dN[E2c_T],&rate[E2c_I2c],dt,&dN[E2c_I2c]);
  reulermultinom(1,E_3_c[u]-dN[E3c_T],&rate[E3c_I3c],dt,&dN[E3c_I3c]);
  reulermultinom(1,E_4_c[u]-dN[E4c_T],&rate[E4c_I4c],dt,&dN[E4c_I4c]);

  // Community I to R/D
  rate[I1c_D] = infection_fatality[u]*gamma_u;
  rate[I2c_D] = infection_fatality[u]*gamma_u;
  rate[I3c_D] = infection_fatality[u]*gamma_u;
  rate[I4c_D] = infection_fatality[u]*gamma_u;
  rate[I1c_R] = (1-infection_fatality[u])*gamma_u;
  rate[I2c_R] = (1-infection_fatality[u])*gamma_u;
  rate[I3c_R] = (1-infection_fatality[u])*gamma_u;
  rate[I4c_R] = (1-infection_fatality[u])*gamma_u;

  reulermultinom(2,I_1_c[u]-dN[I1c_T],&rate[I1c_D],dt,&dN[I1c_D]);
  reulermultinom(2,I_2_c[u]-dN[I2c_T],&rate[I2c_D],dt,&dN[I2c_D]);
  reulermultinom(2,I_3_c[u]-dN[I3c_T],&rate[I3c_D],dt,&dN[I3c_D]);
  reulermultinom(2,I_4_c[u]-dN[I4c_T],&rate[I4c_D],dt,&dN[I4c_D]);

  S[u] += nearbyint(-dN[S_V] - dN[S_E1c] - dN[S_E2c] - dN[S_E3c] - dN[S_E4c]);
  V[u] += nearbyint(dN[S_V] - dN[V_E1c] - dN[V_E2c] - dN[V_E3c] - dN[V_E4c]);
  E_1_c[u] += nearbyint(dN[S_E1c] + dN[V_E1c] - dN[E1c_T] - dN[E1c_I1c]);
  E_2_c[u] += nearbyint(dN[S_E2c] + dN[V_E2c] - dN[E2c_T] - dN[E2c_I2c]);
  E_3_c[u] += nearbyint(dN[S_E3c] + dN[V_E3c] - dN[E3c_T] - dN[E3c_I3c]);
  E_4_c[u] += nearbyint(dN[S_E4c] + dN[V_E4c] - dN[E4c_T] - dN[E4c_I4c]);
  I_1_c[u] += nearbyint(dN[E1c_I1c] - dN[I1c_T] - dN[I1c_R] - dN[I1c_D]);
  I_2_c[u] += nearbyint(dN[E2c_I2c] - dN[I2c_T] - dN[I2c_R] - dN[I2c_D]);
  I_3_c[u] += nearbyint(dN[E3c_I3c] - dN[I3c_T] - dN[I3c_R] - dN[I3c_D]);
  I_4_c[u] += nearbyint(dN[E4c_I4c] - dN[I4c_T] - dN[I4c_R] - dN[I4c_D]);
  // Reset negative values to zero
  E_1_c[u] = (E_1_c[u] < 0) ? 0 : E_1_c[u];
  E_2_c[u] = (E_2_c[u] < 0) ? 0 : E_2_c[u];
  E_3_c[u] = (E_3_c[u] < 0) ? 0 : E_3_c[u];
  E_4_c[u] = (E_4_c[u] < 0) ? 0 : E_4_c[u];
  I_1_c[u] = (I_1_c[u] < 0) ? 0 : I_1_c[u];
  I_2_c[u] = (I_2_c[u] < 0) ? 0 : I_2_c[u];
  I_3_c[u] = (I_3_c[u] < 0) ? 0 : I_3_c[u];
  I_4_c[u] = (I_4_c[u] < 0) ? 0 : I_4_c[u];

  // Travel related transitions, illustrated as below:
  // 1. Direct travel //
  // T                             T              > D
  //  \(T_Eki0)                     \(T_Iki0)    /(Iki0_D)
  //   \>             (Eki0_Iki0)    \>         /     (Iki0_R)
  //  E_k_i_direct_in      -->      I_k_i_direct_in      -->      R
  //     \                             \
  //      \(Eki0_P)                     \(Iki0_P)
  //       \>P                           \>P
  // 2. Transit first //
  // P                                  P                 > D
  //  \(P_Eki1)                          \(P_Iki1)       /(Iki1_D)
  //   \>                  (Eki1_Iki1)    \>            /      (Iki1_R)
  //   E_k_i_transit_first      -->      I_k_i_transit_first      -->      R
  
  // Travel related transitions for direct travelers
  for (v=0; v < U ; v++) { // travelers from v to u
    if(v != u) {
      code_travel_related_transitions_direct_TO_BE_REPLACED;
    }
  }

  // Travel related transitions for transit first travelers
  for (v=0; v < U ; v++) { // travelers from v to u
    if(v != u) {
      code_travel_related_transitions_transit_first_TO_BE_REPLACED;
    }
  }

  // Transitions in R and D compartments
  R[u] += nearbyint(
    dN[I1c_R] + dN[I2c_R] + dN[I3c_R] + dN[I4c_R] +
    code_dN_Ikin_R_o_TO_BE_REPLACED
    );
  D[u] += nearbyint(
    dN[I1c_D] + dN[I2c_D] + dN[I3c_D] + dN[I4c_D] +
    code_dN_Ikin_D_o_TO_BE_REPLACED
    );

  // new cases tracker, reset at obs times
  C_1_c_infected_new[u] = dN[I1c_R] + dN[I1c_D];
  C_2_c_infected_new[u] = dN[I2c_R] + dN[I2c_D];
  C_3_c_infected_new[u] = dN[I3c_R] + dN[I3c_D];
  C_4_c_infected_new[u] = dN[I4c_R] + dN[I4c_D];
  
  code_C_k_i_infected_origin_o_new_TO_BE_REPLACED;
