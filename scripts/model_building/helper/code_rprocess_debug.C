// code for debug
if(day<2){
  if(u==3){
  //if (!R_FINITE(C_new[u])) {
    //print S, E1c, E2c, E3c, E4c
    Rprintf("day=%i, U = %i, S=%lg, E1c=%lg, E2c=%lg, E3c=%lg, E4c=%lg \\n ",day, u, S[u], E_1_c[u], E_2_c[u], E_3_c[u], E_4_c[u]);
    // print different dN
    //Rprintf("dN[R_S]=%lg, dN[S_E1c]=%lg, dN[S_E2c]=%lg, dN[S_E3c]=%lg, dN[S_E4c]=%lg, dN[E1i_in]=%lg, dN[E1i_out]=%lg, dN[E1i_I1i]=%lg, dN[E2i_in]=%lg, dN[E2i_out]=%lg, dN[E2i_I2i]=%lg, dN[E3i_in]=%lg, dN[E3i_out]=%lg, dN[E3i_I3i]=%lg, dN[E4i_in]=%lg, dN[E4i_out]=%lg, dN[E4i_I4i]=%lg, dN[E1c_out]=%lg, dN[E1c_I1c]=%lg, dN[E2c_out]=%lg, dN[E2c_I2c]=%lg, dN[E3c_out]=%lg, dN[E3c_I3c]=%lg, dN[E4c_out]=%lg, dN[E4c_I4c]=%lg, dN[I1i_in]=%lg, dN[I1i_out]=%lg, dN[I1i_R]=%lg, dN[I1i_D]=%lg, dN[I2i_in]=%lg, dN[I2i_out]=%lg, dN[I2i_R]=%lg, dN[I2i_D]=%lg, dN[I3i_in]=%lg, dN[I3i_out]=%lg, dN[I3i_R]=%lg, dN[I3i_D]=%lg, dN[I4i_in]=%lg, dN[I4i_out]=%lg, dN[I4i_R]=%lg, dN[I4i_D]=%lg, dN[I1c_out]=%lg, dN[I1c_R]=%lg, dN[I1c_D]=%lg, dN[I2c_out]=%lg, dN[I2c_R]=%lg, dN[I2c_D]=%lg, dN[I3c_out]=%lg, dN[I3c_R]=%lg \\n", dN[R_S], dN[S_E1c], dN[S_E2c], dN[S_E3c], dN[S_E4c], dN[E1i_in], dN[E1i_out], dN[E1i_I1i], dN[E2i_in], dN[E2i_out], dN[E2i_I2i], dN[E3i_in], dN[E3i_out], dN[E3i_I3i], dN[E4i_in], dN[E4i_out], dN[E4i_I4i], dN[E1c_out], dN[E1c_I1c], dN[E2c_out], dN[E2c_I2c], dN[E3c_out], dN[E3c_I3c], dN[E4c_out], dN[E4c_I4c], dN[I1i_in], dN[I1i_out], dN[I1i_R], dN[I1i_D], dN[I2i_in], dN[I2i_out], dN[I2i_R], dN[I2i_D], dN[I3i_in], dN[I3i_out], dN[I3i_R], dN[I3i_D], dN[I4i_in], dN[I4i_out], dN[I4i_R], dN[I4i_D], dN[I1c_out], dN[I1c_R], dN[I1c_D], dN[I2c_out], dN[I2c_R], dN[I2c_D], dN[I3c_out], dN[I3c_R]);
    // print rate[S_E1c] etc
    Rprintf("day=%i, U = %i, S_E1c=%lg, S_E2c=%lg, S_E3c=%lg, S_E4c=%lg \\n",day, u, rate[S_E1c], rate[S_E2c], rate[S_E3c], rate[S_E4c]);
    // print omega1 etc
    //Rprintf("day=%i, U = %i, omega1=%lg, omega2=%lg, omega3=%lg, omega4=%lg \\n",day, u, omega1, omega2, omega3, omega4);
    //print Itilde1 etc
    Rprintf("day=%i, U = %i, I_tilde1=%lg, I_tilde2=%lg, I_tilde3=%lg, I_tilde4=%lg \\n",day, u, I_tilde1, I_tilde2, I_tilde3, I_tilde4);
    // print I1i, I1c etc
    Rprintf("day=%i, U = %i, I_1_i=%lg, I_2_i=%lg, I_3_i=%lg, I_4_i=%lg \\n",day, u, I_1_i[u], I_2_i[u], I_3_i[u], I_4_i[u]);
    Rprintf("day=%i, U = %i, I_1_c=%lg, I_2_c=%lg, I_3_c=%lg, I_4_c=%lg \\n",day, u, I_1_c[u], I_2_c[u], I_3_c[u], I_4_c[u]);
    // print rate[E1i_I1i] etc
    Rprintf("day=%i, U = %i, E1i_I1i=%lg, E2i_I2i=%lg, E3i_I3i=%lg, E4i_I4i=%lg \\n",day, u, rate[E1i_I1i], rate[E2i_I2i], rate[E3i_I3i], rate[E4i_I4i]);
    //print tau_in_E1i etc
    Rprintf("day=%i, U = %i, tau_in_E1i=%lg, tau_in_E2i=%lg, tau_in_E3i=%lg, tau_in_E4i=%lg \\n",day, u, tau_in_E1i, tau_in_E2i, tau_in_E3i, tau_in_E4i);
    // print tau_out_E1c etc
    Rprintf("day=%i, U = %i, tau_out_E1c=%lg, tau_out_E2c=%lg, tau_out_E3c=%lg, tau_out_E4c=%lg \\n",day, u, tau_out_E1c, tau_out_E2c, tau_out_E3c, tau_out_E4c);
    // print tau_in_I1i etc
    Rprintf("day=%i, U = %i, tau_in_I1i=%lg, tau_in_I2i=%lg, tau_in_I3i=%lg, tau_in_I4i=%lg \\n",day, u, tau_in_I1i, tau_in_I2i, tau_in_I3i, tau_in_I4i);
    // print tau_out_I1c etc
    Rprintf("day=%i, U = %i, tau_out_I1c=%lg, tau_out_I2c=%lg, tau_out_I3c=%lg, tau_out_I4c=%lg \\n",day, u, tau_out_I1c, tau_out_I2c, tau_out_I3c, tau_out_I4c);
  }
}