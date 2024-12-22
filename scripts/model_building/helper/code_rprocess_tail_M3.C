C_k_c_infected_u = C_1_c_infected_new[u] + C_2_c_infected_new[u] + C_3_c_infected_new[u] + C_4_c_infected_new[u];
C_k_i_infected_u = code_sum_C_k_i_infected_TO_BE_REPLACED;

// whether ther is a specified traveler weight
if(traveler_weight[u]<=1.0){
  is_traveler_weight_u = 1;
} else {
  is_traveler_weight_u = 0;
}

infection_detection_new_u = infection_detection[u]*change_IDR[u];
if((C_k_c_infected_u+C_k_i_infected_u)*infection_detection_new_u > diagnostic_capacity[u]*change_IDR[u]){ // if the number of diagnostic cases is larger than the diagnostic capacity, then the number total diagnostic cases will be capped by the diagnostic capacity.
  infection_detection_new_u = diagnostic_capacity[u]*change_IDR[u]/C_k_i_infected_u;
}
if(infection_detection_new_u>1){
  // community detection
  C_1_c_detected_new[u] = C_1_c_infected_new[u];
  C_2_c_detected_new[u] = C_2_c_infected_new[u];
  C_3_c_detected_new[u] = C_3_c_infected_new[u];
  C_4_c_detected_new[u] = C_4_c_infected_new[u];
  C_k_c_detected_u = C_1_c_detected_new[u] + C_2_c_detected_new[u] + C_3_c_detected_new[u] + C_4_c_detected_new[u];
  // imported detection
  code_C_k_i_detected_origin_o_new_IDRover1_TO_BE_REPLACED;
  C_k_i_detected_u = code_sum_C_k_i_detected_TO_BE_REPLACED;
} else {
  // community detection
  C_1_c_detected_new[u] = C_1_c_infected_new[u]*infection_detection_new_u;
  C_2_c_detected_new[u] = C_2_c_infected_new[u]*infection_detection_new_u;
  C_3_c_detected_new[u] = C_3_c_infected_new[u]*infection_detection_new_u;
  C_4_c_detected_new[u] = C_4_c_infected_new[u]*infection_detection_new_u;
  C_k_c_detected_u = C_1_c_detected_new[u] + C_2_c_detected_new[u] + C_3_c_detected_new[u] + C_4_c_detected_new[u];
  // imported detection
  code_C_k_i_detected_origin_o_new_IDRlower1_TO_BE_REPLACED;
  C_k_i_detected_u = code_sum_C_k_i_detected_TO_BE_REPLACED;

  if(is_traveler_weight_u==1){
    // traveler weight means the proportion of detection/sequencing capacity is left for the imported cases.
    // 1. calculate the total number of detected cases
    C_k_all_detected_u = C_k_c_detected_u + C_k_i_detected_u;
    // 2. calculate the proportion of detected cases that left for the imported
    num_C_k_i_propensity_u = C_k_all_detected_u * traveler_weight[u];
    num_C_k_c_propensity_u = C_k_all_detected_u * (1-traveler_weight[u]);
    // update the C_k_c_detected_new and C_k_i_detected_origin_o_new
    if(C_k_i_infected_u==0){
      new_IDR_i_u = 0;
    } else {
      new_IDR_i_u = num_C_k_i_propensity_u / C_k_i_infected_u;
    }
    if(C_k_c_infected_u==0){
      new_IDR_c_u = 0;
    } else {
      new_IDR_c_u = num_C_k_c_propensity_u / C_k_c_infected_u;
    }
    if(new_IDR_i_u > 1){ // all the imported cases can be detected
      code_C_k_i_detected_origin_o_new_IDRover1_TO_BE_REPLACED;
      C_k_i_detected_u = code_sum_C_k_i_detected_TO_BE_REPLACED;
    } else { // only a proportion of imported cases can be detected
      code_C_k_i_detected_origin_o_new_newIDRlower1_TO_BE_REPLACED;
      C_k_i_detected_u = code_sum_C_k_i_detected_TO_BE_REPLACED;
    }
    if(new_IDR_c_u > 1){ // all the community cases can be detected
      C_1_c_detected_new[u] = C_1_c_infected_new[u];
      C_2_c_detected_new[u] = C_2_c_infected_new[u];
      C_3_c_detected_new[u] = C_3_c_infected_new[u];
      C_4_c_detected_new[u] = C_4_c_infected_new[u];
      C_k_c_detected_u = C_1_c_detected_new[u] + C_2_c_detected_new[u] + C_3_c_detected_new[u] + C_4_c_detected_new[u];
    } else { // only a proportion of community cases can be detected
      C_1_c_detected_new[u] = C_1_c_infected_new[u]*new_IDR_c_u;
      C_2_c_detected_new[u] = C_2_c_infected_new[u]*new_IDR_c_u;
      C_3_c_detected_new[u] = C_3_c_infected_new[u]*new_IDR_c_u;
      C_4_c_detected_new[u] = C_4_c_infected_new[u]*new_IDR_c_u;
      C_k_c_detected_u = C_1_c_detected_new[u] + C_2_c_detected_new[u] + C_3_c_detected_new[u] + C_4_c_detected_new[u];
    }
  }
}

detection_sequencing_new_u = ((C_k_i_infected_u + C_k_c_infected_u) * infection_detection[u] * detection_sequencing[u])/(C_k_i_detected_u + C_k_c_detected_u)*change_DSR[u]; // Note that '(C_k_i_infected_u + C_k_c_infected_u) * infection_detection[u] * detection_sequencing[u]' is the original sequenced number, '(C_k_i_detected_u + C_k_c_detected_u)' is the updated detection number. The detection sequencing ratio is therefore not dependent on the changes of infection detection ratio.
if((C_k_i_detected_u + C_k_c_detected_u)*detection_sequencing_new_u > sequencing_capacity[u]*change_DSR[u]){ // if the number of sequenced cases is larger than the sequencing capacity, then the number total sequenced cases will be capped by the sequencing capacity.
  detection_sequencing_new_u = sequencing_capacity[u]*change_DSR[u]/(C_k_i_detected_u + C_k_c_detected_u);
}

if(detection_sequencing_new_u>1){
  // community sequencing, will be fitted, N=4
  C_1_c_sequenced_new[u] = C_1_c_detected_new[u];
  C_2_c_sequenced_new[u] = C_2_c_detected_new[u];
  C_3_c_sequenced_new[u] = C_3_c_detected_new[u];
  C_4_c_sequenced_new[u] = C_4_c_detected_new[u];
  C_k_c_sequenced_u = C_1_c_sequenced_new[u] + C_2_c_sequenced_new[u] + C_3_c_sequenced_new[u] + C_4_c_sequenced_new[u];
  // imported sequencing, will be fitted, N=116
  code_C_k_i_sequenced_origin_o_new_DSRover1_TO_BE_REPLACED;
  C_k_i_sequenced_u = code_sum_C_k_i_sequenced_TO_BE_REPLACED;
} else {
  // community sequencing, will be fitted, N=4
  C_1_c_sequenced_new[u] = C_1_c_detected_new[u]*detection_sequencing_new_u;
  C_2_c_sequenced_new[u] = C_2_c_detected_new[u]*detection_sequencing_new_u;
  C_3_c_sequenced_new[u] = C_3_c_detected_new[u]*detection_sequencing_new_u;
  C_4_c_sequenced_new[u] = C_4_c_detected_new[u]*detection_sequencing_new_u;
  C_k_c_sequenced_u = C_1_c_sequenced_new[u] + C_2_c_sequenced_new[u] + C_3_c_sequenced_new[u] + C_4_c_sequenced_new[u];
  // imported sequencing, will be fitted, N=116
  code_C_k_i_sequenced_origin_o_new_DSRlower1_TO_BE_REPLACED;
  C_k_i_sequenced_u = code_sum_C_k_i_sequenced_TO_BE_REPLACED;

  if(is_traveler_weight_u==1){
    // traveler weight means the proportion of detection/sequencing capacity is left for the imported cases.
    // 1. calculate the total number of sequenced cases
    C_k_all_sequenced_u = C_k_c_sequenced_u + C_k_i_sequenced_u;
    // 2. calculate the proportion of sequenced cases that left for the imported
    num_C_k_i_propensity_u = C_k_all_sequenced_u * traveler_weight[u];
    num_C_k_c_propensity_u = C_k_all_sequenced_u * (1-traveler_weight[u]);
    // update the C_k_c_sequenced_new and C_k_i_sequenced_origin_o_new
    if(C_k_i_detected_u==0){
      new_DSR_i_u = 0;
    } else {
      new_DSR_i_u = num_C_k_i_propensity_u / C_k_i_detected_u;
    }
    if(C_k_c_detected_u==0){
      new_DSR_c_u = 0;
    } else {
      new_DSR_c_u = num_C_k_c_propensity_u / C_k_c_detected_u;
    }
    if(new_DSR_i_u > 1){ // all the imported cases can be sequenced
      code_C_k_i_sequenced_origin_o_new_DSRover1_TO_BE_REPLACED;
    } else { // only a proportion of imported cases can be sequenced
      code_C_k_i_sequenced_origin_o_new_newDSRlower1_TO_BE_REPLACED;
    }
    if(new_DSR_c_u > 1){ // all the community cases can be sequenced
      C_1_c_sequenced_new[u] = C_1_c_detected_new[u];
      C_2_c_sequenced_new[u] = C_2_c_detected_new[u];
      C_3_c_sequenced_new[u] = C_3_c_detected_new[u];
      C_4_c_sequenced_new[u] = C_4_c_detected_new[u];
    } else { // only a proportion of community cases can be sequenced
      C_1_c_sequenced_new[u] = C_1_c_detected_new[u]*new_DSR_c_u;
      C_2_c_sequenced_new[u] = C_2_c_detected_new[u]*new_DSR_c_u;
      C_3_c_sequenced_new[u] = C_3_c_detected_new[u]*new_DSR_c_u;
      C_4_c_sequenced_new[u] = C_4_c_detected_new[u]*new_DSR_c_u;
    }
  }
}

// All sequenced, will be fitted, N=4
code_C_1_all_sequenced_new_TO_BE_REPLACED;
code_C_2_all_sequenced_new_TO_BE_REPLACED;
code_C_3_all_sequenced_new_TO_BE_REPLACED;
code_C_4_all_sequenced_new_TO_BE_REPLACED;

// All deceased, will be fitted, N=1
D_new[u] += dN[I1c_D] + dN[I2c_D] + dN[I3c_D] + dN[I4c_D] + code_dN_Ikin_D_o_TO_BE_REPLACED;
// All reported, will be fitted, N=1
C_reported_new[u] += C_k_c_detected_u + C_k_i_detected_u - D_new[u];
if(C_reported_new[u] < 0) C_reported_new[u] = 0;

} // intentionally added bracket to close the loop
