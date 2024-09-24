// community detection
C_1_c_detected_new[u] = C_1_c_infected_new[u]*infection_detection[u];
C_2_c_detected_new[u] = C_2_c_infected_new[u]*infection_detection[u];
C_3_c_detected_new[u] = C_3_c_infected_new[u]*infection_detection[u];
C_4_c_detected_new[u] = C_4_c_infected_new[u]*infection_detection[u];

// community sequencing, will be fitted, N=4
C_1_c_sequenced_new[u] = C_1_c_detected_new[u]*detection_sequencing[u];
C_2_c_sequenced_new[u] = C_2_c_detected_new[u]*detection_sequencing[u];
C_3_c_sequenced_new[u] = C_3_c_detected_new[u]*detection_sequencing[u];
C_4_c_sequenced_new[u] = C_4_c_detected_new[u]*detection_sequencing[u];

// imported detection
code_C_k_i_detected_origin_o_new_TO_BE_REPLACED;
// imported sequencing, will be fitted, N=116
code_C_k_i_sequenced_origin_o_new_TO_BE_REPLACED;

// All sequenced, will be fitted, N=4
code_C_1_all_sequenced_new_TO_BE_REPLACED;
code_C_2_all_sequenced_new_TO_BE_REPLACED;
code_C_3_all_sequenced_new_TO_BE_REPLACED;
code_C_4_all_sequenced_new_TO_BE_REPLACED;

// All deceased, will be fitted, N=1
D_new[u] = dN[I1c_D] + dN[I2c_D] + dN[I3c_D] + dN[I4c_D] + code_dN_Ikin_D_o_TO_BE_REPLACED;
// All reported, will be fitted, N=1
C_reported_new[u] = C_1_c_detected_new[u] + C_2_c_detected_new[u] + C_3_c_detected_new[u] + C_4_c_detected_new[u] + code_C_k_i_detected_all_TO_BE_REPLACED - D_new[u];
if(C_reported_new[u] < 0) C_reported_new[u] = 0;

} // intentionally added bracket to close the loop
