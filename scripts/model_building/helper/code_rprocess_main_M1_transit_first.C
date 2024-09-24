// routes_E_transit_first_out
// // strain 1
rate[E1i1_I1i1] = epsilon_one_u;
reulermultinom(1,E_1_i_transit_first[u],&rate[E1i1_I1i1],dt,&dN[E1i1_I1i1]);
// // strain 2
rate[E2i1_I2i1] = epsilon_two_u;
reulermultinom(1,E_2_i_transit_first[u],&rate[E2i1_I2i1],dt,&dN[E2i1_I2i1]);
// // strain 3
rate[E3i1_I3i1] = epsilon_three_u;
reulermultinom(1,E_3_i_transit_first[u],&rate[E3i1_I3i1],dt,&dN[E3i1_I3i1]);
// // strain 4
rate[E4i1_I4i1] = epsilon_four_u;
reulermultinom(1,E_4_i_transit_first[u],&rate[E4i1_I4i1],dt,&dN[E4i1_I4i1]);

// routes_I_transit_first_out
// // strain 1
rate[I1i1_D] = infection_fatality[u] * gamma_u;
rate[I1i1_R] = (1-infection_fatality[u]) * gamma_u;
reulermultinom(2,I_1_i_transit_first[u],&rate[I1i1_D],dt,&dN[I1i1_D]);
// // strain 2
rate[I2i1_D] = infection_fatality[u] * gamma_u;
rate[I2i1_R] = (1-infection_fatality[u]) * gamma_u;
reulermultinom(2,I_2_i_transit_first[u],&rate[I2i1_D],dt,&dN[I2i1_D]);
// // strain 3
rate[I3i1_D] = infection_fatality[u] * gamma_u;
rate[I3i1_R] = (1-infection_fatality[u]) * gamma_u;
reulermultinom(2,I_3_i_transit_first[u],&rate[I3i1_D],dt,&dN[I3i1_D]);
// // strain 4
rate[I4i1_D] = infection_fatality[u] * gamma_u;
rate[I4i1_R] = (1-infection_fatality[u]) * gamma_u;
reulermultinom(2,I_4_i_transit_first[u],&rate[I4i1_D],dt,&dN[I4i1_D]);

// assigning dNs
// // strain 1
E_1_i_transit_first[u] += nearbyint(-dN[E1i1_I1i1]);
E_1_i_transit_first[u] = (E_1_i_transit_first[u]<0) ? 0 : E_1_i_transit_first[u];
I_1_i_transit_first[u] += nearbyint(dN[E1i1_I1i1] - dN[I1i1_D] - dN[I1i1_R]);
I_1_i_transit_first[u] = (I_1_i_transit_first[u]<0) ? 0 : I_1_i_transit_first[u];
// // strain 2
E_2_i_transit_first[u] += nearbyint(-dN[E2i1_I2i1]);
E_2_i_transit_first[u] = (E_2_i_transit_first[u]<0) ? 0 : E_2_i_transit_first[u];
I_2_i_transit_first[u] += nearbyint(dN[E2i1_I2i1] - dN[I2i1_D] - dN[I2i1_R]);
I_2_i_transit_first[u] = (I_2_i_transit_first[u]<0) ? 0 : I_2_i_transit_first[u];
// // strain 3
E_3_i_transit_first[u] += nearbyint(-dN[E3i1_I3i1]);
E_3_i_transit_first[u] = (E_3_i_transit_first[u]<0) ? 0 : E_3_i_transit_first[u];
I_3_i_transit_first[u] += nearbyint(dN[E3i1_I3i1] - dN[I3i1_D] - dN[I3i1_R]);
I_3_i_transit_first[u] = (I_3_i_transit_first[u]<0) ? 0 : I_3_i_transit_first[u];
// // strain 4
E_4_i_transit_first[u] += nearbyint(-dN[E4i1_I4i1]);
E_4_i_transit_first[u] = (E_4_i_transit_first[u]<0) ? 0 : E_4_i_transit_first[u];
I_4_i_transit_first[u] += nearbyint(dN[E4i1_I4i1] - dN[I4i1_D] - dN[I4i1_R]);
I_4_i_transit_first[u] = (I_4_i_transit_first[u]<0) ? 0 : I_4_i_transit_first[u];
