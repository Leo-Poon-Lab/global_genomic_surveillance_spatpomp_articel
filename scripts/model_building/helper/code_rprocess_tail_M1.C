// Detection
C_1_c_detected_new[u] = C_1_c_infected_new[u]*infection_detection[u];
C_2_c_detected_new[u] = C_2_c_infected_new[u]*infection_detection[u];
C_3_c_detected_new[u] = C_3_c_infected_new[u]*infection_detection[u];
C_4_c_detected_new[u] = C_4_c_infected_new[u]*infection_detection[u];

C_1_i_detected_new[u] = C_1_i_infected_new[u]*infection_detection[u];
C_2_i_detected_new[u] = C_2_i_infected_new[u]*infection_detection[u];
C_3_i_detected_new[u] = C_3_i_infected_new[u]*infection_detection[u];
C_4_i_detected_new[u] = C_4_i_infected_new[u]*infection_detection[u];

C_1_all_detected_new[u] = C_1_c_detected_new[u] + C_1_i_detected_new[u];
C_2_all_detected_new[u] = C_2_c_detected_new[u] + C_2_i_detected_new[u];
C_3_all_detected_new[u] = C_3_c_detected_new[u] + C_3_i_detected_new[u];
C_4_all_detected_new[u] = C_4_c_detected_new[u] + C_4_i_detected_new[u];

// Sequencing
C_1_c_sequenced_new[u] = C_1_c_detected_new[u]*detection_sequencing[u];
C_2_c_sequenced_new[u] = C_2_c_detected_new[u]*detection_sequencing[u];
C_3_c_sequenced_new[u] = C_3_c_detected_new[u]*detection_sequencing[u];
C_4_c_sequenced_new[u] = C_4_c_detected_new[u]*detection_sequencing[u];

C_1_i_sequenced_new[u] = C_1_i_detected_new[u]*detection_sequencing[u];
C_2_i_sequenced_new[u] = C_2_i_detected_new[u]*detection_sequencing[u];
C_3_i_sequenced_new[u] = C_3_i_detected_new[u]*detection_sequencing[u];
C_4_i_sequenced_new[u] = C_4_i_detected_new[u]*detection_sequencing[u];

C_1_all_sequenced_new[u] = C_1_c_sequenced_new[u] + C_1_i_sequenced_new[u]; // will be fitted, N=1
C_2_all_sequenced_new[u] = C_2_c_sequenced_new[u] + C_2_i_sequenced_new[u]; // will be fitted, N=1
C_3_all_sequenced_new[u] = C_3_c_sequenced_new[u] + C_3_i_sequenced_new[u]; // will be fitted, N=1
C_4_all_sequenced_new[u] = C_4_c_sequenced_new[u] + C_4_i_sequenced_new[u]; // will be fitted, N=1

// All deceased, will be fitted, N=1
D_new[u] += dN[I1c_D] + dN[I2c_D] + dN[I3c_D] + dN[I4c_D] + dN[I1i0_D] + dN[I2i0_D] + dN[I3i0_D] + dN[I4i0_D] + dN[I1i1_D] + dN[I2i1_D] + dN[I3i1_D] + dN[I4i1_D];
// All reported, will be fitted, N=1
C_reported_new[u] += C_1_all_detected_new[u] + C_2_all_detected_new[u] + C_3_all_detected_new[u] + C_4_all_detected_new[u] - D_new[u];
if(C_reported_new[u] < 0) C_reported_new[u] = 0;

} // intentionally added bracket to close the loop
