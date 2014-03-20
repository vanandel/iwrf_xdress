/////////////////////////////////////////////////////////////////
// iwrf_functions.c
//
// Utility routines for iwrf_data structs
//
// Mike Dixon, RAL, NCAR, POBox 3000, Boulder, CO, 80307-3000
// Feb 2009
// ported for C compiler -- D. Brunkow CSU-CHILL  Sep 2009

#include "dataport/swap.h"
// #include <toolsa/DateTime.hh>
// #include <toolsa/mem.h>
// #include <toolsa/str.h>
#include <stdlib.h>
#include <stdio.h>
#include <sys/time.h>
#include <time.h>
#include "iwrf_functions.h"
#include "iwrf_user_interface.h"
#include <math.h>
#ifndef NAN
#define NAN (0.0/0.0)
#endif

#define MEM_zero(a) memset(a, 0, sizeof(*a))

//using namespace std;

/////////////////////////////
// check for a missing values

int iwrf_int_is_missing(si32 val)
{
  if (val == IWRF_MISSING_INT) {
    return true;
  } else {
    return false;
  }
}

int iwrf_float_is_missing(fl32 val)
{
  if (isnan(val)) {
    return true;
  } else if (fabs(val - IWRF_MISSING_FLOAT) < 0.001)  {
    return true;
  } else {
    return false;
  }
}


//////////////////////////////////////////////////////////////////
// packet initialization
// sets values to missing as appropriate

//////////////////////////////////////////////////////
// init sync struct

void iwrf_sync_init(iwrf_sync_t *val)

{

  MEM_zero(val);
  val->packet.id = IWRF_SYNC_ID;
  val->packet.len_bytes = sizeof(*val);
  val->packet.version_num = 1;
  iwrf_set_packet_time_to_now(&val->packet);

  val->magik[0] = IWRF_SYNC_VAL_00;
  val->magik[1] = IWRF_SYNC_VAL_01;

}

//////////////////////////////////////////////////////
// init version struct

void iwrf_version_init(iwrf_version_t *val)

{

  MEM_zero(val);

  val->packet.id = IWRF_VERSION_ID;
  val->packet.len_bytes = sizeof(*val);
  val->packet.version_num = 1;
  iwrf_set_packet_time_to_now(&val->packet);
  
}

//////////////////////////////////////////////////////
// init radar_info struct

void iwrf_radar_info_init(iwrf_radar_info_t *val)

{

  MEM_zero(val);

  val->packet.id = IWRF_RADAR_INFO_ID;
  val->packet.len_bytes = sizeof(*val);
  val->packet.version_num = 1;
  iwrf_set_packet_time_to_now(&val->packet);

  val->latitude_deg = IWRF_MISSING_FLOAT;
  val->longitude_deg = IWRF_MISSING_FLOAT;
  val->altitude_m = IWRF_MISSING_FLOAT;
  val->platform_type = IWRF_RADAR_PLATFORM_NOT_SET;

  val->beamwidth_deg_h = IWRF_MISSING_FLOAT;
  val->beamwidth_deg_v = IWRF_MISSING_FLOAT;
  val->wavelength_cm = IWRF_MISSING_FLOAT;
  
  val->nominal_gain_ant_db_h = IWRF_MISSING_FLOAT;
  val->nominal_gain_ant_db_v = IWRF_MISSING_FLOAT;

}

//////////////////////////////////////////////////////
// init scan_segment struct

void iwrf_scan_segment_init(iwrf_scan_segment_t *val)

{

  MEM_zero(val);

  val->packet.id = IWRF_SCAN_SEGMENT_ID;
  val->packet.len_bytes = sizeof(*val);
  val->packet.version_num = 1;
  iwrf_set_packet_time_to_now(&val->packet);

  val->scan_mode = IWRF_SCAN_MODE_NOT_SET;
  val->follow_mode = IWRF_FOLLOW_MODE_NOT_SET;
  val->volume_num = IWRF_MISSING_INT;
  val->sweep_num = IWRF_MISSING_INT;
  val->time_limit = 0;

  val->az_manual = IWRF_MISSING_FLOAT;
  val->el_manual = IWRF_MISSING_FLOAT;
  val->az_start = IWRF_MISSING_FLOAT;
  val->el_start = IWRF_MISSING_FLOAT;
  val->scan_rate = IWRF_MISSING_FLOAT;
  val->left_limit = IWRF_MISSING_FLOAT;
  val->right_limit = IWRF_MISSING_FLOAT;
  val->up_limit = IWRF_MISSING_FLOAT;
  val->down_limit = IWRF_MISSING_FLOAT;
  val->step = IWRF_MISSING_FLOAT;

  val->current_fixed_angle = IWRF_MISSING_FLOAT;
  val->init_direction_cw = 1;
  val->init_direction_up = 1;

  val->n_sweeps = 0;
  
  val->optimizer_rmax_km = IWRF_MISSING_FLOAT;
  val->optimizer_htmax_km = IWRF_MISSING_FLOAT;
  val->optimizer_res_m = IWRF_MISSING_FLOAT;  

}

//////////////////////////////////////////////////////
// init antenna_correction struct

void iwrf_antenna_correction_init(iwrf_antenna_correction_t *val)

{

  MEM_zero(val);

  val->packet.id = IWRF_ANTENNA_CORRECTION_ID;
  val->packet.len_bytes = sizeof(*val);
  val->packet.version_num = 1;
  iwrf_set_packet_time_to_now(&val->packet);

  val->az_correction = IWRF_MISSING_FLOAT;
  val->el_correction = IWRF_MISSING_FLOAT;

}

//////////////////////////////////////////////////////
// init ts_processing struct

void iwrf_ts_processing_init(iwrf_ts_processing_t *val)

{

  MEM_zero(val);

  val->packet.id = IWRF_TS_PROCESSING_ID;
  val->packet.len_bytes = sizeof(*val);
  val->packet.version_num = 1;
  iwrf_set_packet_time_to_now(&val->packet);

  val->xmit_rcv_mode = IWRF_XMIT_RCV_MODE_NOT_SET;
  val->xmit_phase_mode = IWRF_XMIT_PHASE_MODE_NOT_SET;
  val->prf_mode = IWRF_PRF_MODE_NOT_SET;
  val->pulse_type = IWRF_PULSE_TYPE_NOT_SET;

  val->prt_usec = IWRF_MISSING_FLOAT;
  val->prt2_usec = IWRF_MISSING_FLOAT;

  val->cal_type = IWRF_CAL_TYPE_NOT_SET;
  
  val->burst_range_offset_m = IWRF_MISSING_FLOAT;
  val->pulse_width_us = IWRF_MISSING_FLOAT;
  val->start_range_m = IWRF_MISSING_FLOAT;
  val->gate_spacing_m = IWRF_MISSING_FLOAT;

  val->integration_cycle_pulses = IWRF_MISSING_INT;
  val->clutter_filter_number = IWRF_MISSING_INT;
  val->range_gate_averaging = 1;
  val->max_gate = 0;

  val->test_power_dbm = IWRF_MISSING_FLOAT;
  val->test_pulse_range_km = IWRF_MISSING_FLOAT;
  val->test_pulse_length_usec = IWRF_MISSING_FLOAT;

  val->num_prts = IWRF_MISSING_INT;
  val->prt2_usec = IWRF_MISSING_FLOAT;
  val->prt4_usec = IWRF_MISSING_FLOAT;
  val->block_mode_prt2_pulses = IWRF_MISSING_INT;
  val->block_mode_prt3_pulses = IWRF_MISSING_INT;
  val->block_mode_prt4_pulses = IWRF_MISSING_INT;
  val->pol_sync_mode = 0;  // No synchronization by default
  
}

//////////////////////////////////////////////////////
// init xmit_power struct

void iwrf_xmit_power_init(iwrf_xmit_power_t *val)

{

  MEM_zero(val);

  val->packet.id = IWRF_XMIT_POWER_ID;
  val->packet.len_bytes = sizeof(*val);
  val->packet.version_num = 1;
  iwrf_set_packet_time_to_now(&val->packet);

  val->power_dbm_h = IWRF_MISSING_FLOAT;
  val->power_dbm_v = IWRF_MISSING_FLOAT;

}

//////////////////////////////////////////////////////
// init xmit_sample struct

void iwrf_xmit_sample_init(iwrf_xmit_sample_t *val)

{

  MEM_zero(val);

  val->packet.id = IWRF_XMIT_SAMPLE_ID;
  val->packet.len_bytes = sizeof(*val);
  val->packet.version_num = 1;
  iwrf_set_packet_time_to_now(&val->packet);

  val->power_dbm_h = IWRF_MISSING_FLOAT;
  val->power_dbm_v = IWRF_MISSING_FLOAT;
  val->offset = 0;
  val->n_samples = 0;

  val->sampling_freq = IWRF_MISSING_FLOAT;

  val->scale_h = IWRF_MISSING_FLOAT;
  val->offset_h = IWRF_MISSING_FLOAT;

  val->scale_v = IWRF_MISSING_FLOAT;
  val->offset_v = IWRF_MISSING_FLOAT;

}

//////////////////////////////////////////////////////
// init xmit_info struct

void iwrf_xmit_info_init(iwrf_xmit_info_t *val)

{

  MEM_zero(val);

  val->packet.id = IWRF_XMIT_INFO_ID;
  val->packet.len_bytes = sizeof(*val);
  val->packet.version_num = 1;
  iwrf_set_packet_time_to_now(&val->packet);

  val->xmit_0_enabled = 0;
  val->xmit_1_enabled = 0;

  val->xmit_rcv_mode = IWRF_XMIT_RCV_MODE_NOT_SET;
  val->xmit_phase_mode = IWRF_XMIT_PHASE_MODE_NOT_SET;
  val->prf_mode = IWRF_PRF_MODE_NOT_SET;
  val->pulse_type = IWRF_PULSE_TYPE_NOT_SET;

  val->prt_usec = IWRF_MISSING_FLOAT;
  val->prt2_usec = IWRF_MISSING_FLOAT;

}

//////////////////////////////////////////////////////
// init calibration struct

void iwrf_calibration_init(iwrf_calibration_t *val)

{

  MEM_zero(val);

  val->packet.id = IWRF_CALIBRATION_ID;
  val->packet.len_bytes = sizeof(*val);
  val->packet.version_num = 1;
  iwrf_set_packet_time_to_now(&val->packet);

  val->wavelength_cm = IWRF_MISSING_FLOAT;
  val->beamwidth_deg_h = IWRF_MISSING_FLOAT;
  val->beamwidth_deg_v = IWRF_MISSING_FLOAT;
  val->gain_ant_db_h = IWRF_MISSING_FLOAT;
  val->gain_ant_db_v = IWRF_MISSING_FLOAT;
  val->pulse_width_us = IWRF_MISSING_FLOAT;
  val->xmit_power_dbm_h = IWRF_MISSING_FLOAT;
  val->xmit_power_dbm_v = IWRF_MISSING_FLOAT;
  val->two_way_waveguide_loss_db_h = IWRF_MISSING_FLOAT;
  val->two_way_waveguide_loss_db_v = IWRF_MISSING_FLOAT;
  val->two_way_radome_loss_db_h = IWRF_MISSING_FLOAT;
  val->two_way_radome_loss_db_v = IWRF_MISSING_FLOAT;
  val->receiver_mismatch_loss_db = IWRF_MISSING_FLOAT;
  val->radar_constant_h = IWRF_MISSING_FLOAT;
  val->radar_constant_v = IWRF_MISSING_FLOAT;
  val->noise_dbm_hc = IWRF_MISSING_FLOAT;
  val->noise_dbm_hx = IWRF_MISSING_FLOAT;
  val->noise_dbm_vc = IWRF_MISSING_FLOAT;
  val->noise_dbm_vx = IWRF_MISSING_FLOAT;
  val->receiver_gain_db_hc = IWRF_MISSING_FLOAT;
  val->receiver_gain_db_hx = IWRF_MISSING_FLOAT;
  val->receiver_gain_db_vc = IWRF_MISSING_FLOAT;
  val->receiver_gain_db_vx = IWRF_MISSING_FLOAT;
  val->base_dbz_1km_hc = IWRF_MISSING_FLOAT;
  val->base_dbz_1km_hx = IWRF_MISSING_FLOAT;
  val->base_dbz_1km_vc = IWRF_MISSING_FLOAT;
  val->base_dbz_1km_vx = IWRF_MISSING_FLOAT;
  val->sun_power_dbm_hc = IWRF_MISSING_FLOAT;
  val->sun_power_dbm_hx = IWRF_MISSING_FLOAT;
  val->sun_power_dbm_vc = IWRF_MISSING_FLOAT;
  val->sun_power_dbm_vx = IWRF_MISSING_FLOAT;
  val->noise_source_power_dbm_h = IWRF_MISSING_FLOAT;
  val->noise_source_power_dbm_v = IWRF_MISSING_FLOAT;
  val->power_meas_loss_db_h = IWRF_MISSING_FLOAT;
  val->power_meas_loss_db_v = IWRF_MISSING_FLOAT;
  val->coupler_forward_loss_db_h = IWRF_MISSING_FLOAT;
  val->coupler_forward_loss_db_v = IWRF_MISSING_FLOAT;
  val->test_power_dbm_h = IWRF_MISSING_FLOAT;
  val->test_power_dbm_v = IWRF_MISSING_FLOAT;
  val->zdr_correction_db = IWRF_MISSING_FLOAT;
  val->ldr_correction_db_h = IWRF_MISSING_FLOAT;
  val->ldr_correction_db_v = IWRF_MISSING_FLOAT;
  val->phidp_rot_deg = IWRF_MISSING_FLOAT;
  
}

//////////////////////////////////////////////////////
// init event_notice struct

void iwrf_event_notice_init(iwrf_event_notice_t *val)

{

  MEM_zero(val);

  val->packet.id = IWRF_EVENT_NOTICE_ID;
  val->packet.len_bytes = sizeof(*val);
  val->packet.version_num = 1;
  iwrf_set_packet_time_to_now(&val->packet);

  val->start_of_sweep = 0;
  val->end_of_sweep = 0;

  val->start_of_volume = 0;
  val->end_of_volume = 0;
  
  val->scan_mode = IWRF_SCAN_MODE_NOT_SET;
  val->follow_mode = IWRF_FOLLOW_MODE_NOT_SET;
  val->volume_num = IWRF_MISSING_INT;
  val->sweep_num = IWRF_MISSING_INT;
  
  val->cause = IWRF_EVENT_CAUSE_NOT_SET;

}

//////////////////////////////////////////////////////
// init phasecode struct

void iwrf_phasecode_init(iwrf_phasecode_t *val)

{

  MEM_zero(val);

  val->packet.id = IWRF_PHASECODE_ID;
  val->packet.len_bytes = sizeof(*val);
  val->packet.version_num = 1;
  iwrf_set_packet_time_to_now(&val->packet);

  val->seq_length = 0;

}

//////////////////////////////////////////////////////
// init pulse_header struct

void iwrf_pulse_header_init(iwrf_pulse_header_t *val)

{

  MEM_zero(val);

  val->packet.id = IWRF_PULSE_HEADER_ID;
  val->packet.len_bytes = sizeof(*val);
  val->packet.version_num = 1;
  iwrf_set_packet_time_to_now(&val->packet);

  val->pulse_seq_num = 0;

  val->scan_mode = IWRF_SCAN_MODE_NOT_SET;
  val->follow_mode = IWRF_FOLLOW_MODE_NOT_SET;
  val->sweep_num = IWRF_MISSING_INT;
  val->volume_num = IWRF_MISSING_INT;

  val->fixed_el = IWRF_MISSING_FLOAT;
  val->fixed_az = IWRF_MISSING_FLOAT;
  val->elevation = IWRF_MISSING_FLOAT;
  val->azimuth = IWRF_MISSING_FLOAT;
  
  val->prt = IWRF_MISSING_FLOAT;
  val->prt_next = IWRF_MISSING_FLOAT;
  
  val->pulse_width_us = IWRF_MISSING_FLOAT;

  val->n_gates = 0;

  val->n_channels = 1;
  val->iq_encoding = IWRF_IQ_ENCODING_NOT_SET;
  val->hv_flag = IWRF_MISSING_INT;

  val->antenna_transition = 0;
  val->phase_cohered = IWRF_MISSING_INT;
  val->status = 0;
  val->n_data = 0;
  
  int ichan = 0;
  for (; ichan < IWRF_MAX_CHAN; ichan++) {
    val->iq_offset[ichan] = 0;
    val->burst_mag[ichan] = IWRF_MISSING_FLOAT;
    val->burst_arg[ichan] = IWRF_MISSING_FLOAT;
    val->burst_arg_diff[ichan] = IWRF_MISSING_FLOAT;
  }

  val->scale = 1.0;
  val->offset = 0.0;
  val->n_gates_burst = 0;

}

//////////////////////////////////////////////////////
// init rvp8_ops_info struct

void iwrf_rvp8_ops_info_init(iwrf_rvp8_ops_info_t *val)

{

  MEM_zero(val);

  val->packet.id = IWRF_RVP8_OPS_INFO_ID;
  val->packet.len_bytes = sizeof(*val);
  val->packet.version_num = 1;
  iwrf_set_packet_time_to_now(&val->packet);

  val->i_version = 0;

  val->i_major_mode = IWRF_MISSING_INT;
  val->i_polarization = IWRF_MISSING_INT;
  val->i_phase_mode_seq = IWRF_MISSING_INT;
  
  val->i_task_sweep = IWRF_MISSING_INT;
  val->i_task_aux_num = IWRF_MISSING_INT;
  val->i_task_scan_type = IWRF_MISSING_INT;

  val->i_aq_mode = IWRF_MISSING_INT;
  val->i_unfold_mode = IWRF_MISSING_INT;
  
  val->i_pwidth_code = IWRF_MISSING_INT;
  val->f_pwidth_usec = IWRF_MISSING_FLOAT;

  val->f_dbz_calib = IWRF_MISSING_FLOAT;
  val->i_sample_size = IWRF_MISSING_INT;
  val->i_mean_angle_sync = IWRF_MISSING_INT;

  val->i_flags = 0;

  val->i_playback_version = IWRF_MISSING_INT;

  val->f_sy_clk_mhz = IWRF_MISSING_FLOAT;
  val->f_wavelength_cm = IWRF_MISSING_FLOAT;
  val->f_saturation_dbm = 6;
  
  val->f_range_mask_res = IWRF_MISSING_FLOAT;

  int ichan = 0;
  for (; ichan < IWRF_MAX_CHAN; ichan++) {
    val->f_noise_dbm[ichan] = IWRF_MISSING_FLOAT;
    val->f_noise_stdv_db[ichan] = IWRF_MISSING_FLOAT;
  }
  val->f_noise_range_km = IWRF_MISSING_FLOAT;
  val->f_noise_prf_hz = IWRF_MISSING_FLOAT;

  int ii = 0;
  for (ii = 0; ii < 2; ii++) {
    val->i_gparm_latch_sts[ii] = IWRF_MISSING_INT;
  }
  for (ii = 0; ii < 6; ii++) {
    val->i_gparm_immed_sts[ii] = IWRF_MISSING_INT;
  }
  for (ii = 0; ii < 4; ii++) {
    val->i_gparm_diag_bits[ii] = IWRF_MISSING_INT;
  }

}

//////////////////////////////////////////////////////
// init rvp8_pulse_header struct

void iwrf_rvp8_pulse_header_init(iwrf_rvp8_pulse_header_t *val)

{

  MEM_zero(val);

  val->packet.id = IWRF_RVP8_PULSE_HEADER_ID;
  val->packet.len_bytes = sizeof(*val);
  val->packet.version_num = 1;
  iwrf_set_packet_time_to_now(&val->packet);

  val->i_version = 0;
  val->i_flags = 0;
  val->i_aq_mode = 0;
  val->i_polar_bits = 0;
  val->i_viq_per_bin = 1;
  val->i_tg_bank = 1;

  val->i_tx_phase = IWRF_MISSING_INT;
  val->i_az = IWRF_MISSING_INT;
  val->i_el = IWRF_MISSING_INT;
  val->i_num_vecs = IWRF_MISSING_INT;
  val->i_max_vecs = IWRF_MISSING_INT;
  val->i_tg_wave = IWRF_MISSING_INT;
  
  val->i_btime_api = IWRF_MISSING_INT;
  val->i_sys_time = IWRF_MISSING_INT;
  val->i_prev_prt = IWRF_MISSING_INT;
  val->i_next_prt = IWRF_MISSING_INT;
  val->i_seq_num = IWRF_MISSING_INT;

  int ii = 0;
  for (; ii < 2; ii++) {
    val->uiq_perm[ii] = IWRF_MISSING_INT;
    val->uiq_once[ii] = IWRF_MISSING_INT;
  }

  int ichan = 0;
  for (; ichan < IWRF_MAX_CHAN; ichan++) {
    val->i_data_off[ichan] = 0;
    val->f_burst_mag[ichan] = IWRF_MISSING_FLOAT;
    val->i_burst_arg[ichan] = IWRF_MISSING_INT;
    val->i_wrap_iq[ichan] = IWRF_MISSING_INT;
  }

}

////////////////////////////////////////////////////////////
// set packet sequence number

void iwrf_set_packet_seq_num(iwrf_packet_info_t *packet, si64 seq_num) {
  packet->seq_num = seq_num;
}

////////////////////////////////////////////////////////////
// set packet time

void iwrf_set_packet_time(iwrf_packet_info_t *packet,
			  double dtime) {
  time_t secs = (time_t) dtime;
  int nano_secs = (int) ((dtime - secs) * 1.0e9 + 0.5);
  packet->time_secs_utc = secs;
  packet->time_nano_secs = nano_secs;
}

void iwrf_set_packet_time2(iwrf_packet_info_t *packet,
			  time_t secs, int nano_secs) {
  packet->time_secs_utc = secs;
  packet->time_nano_secs = nano_secs;
}

void iwrf_set_packet_time_to_now(iwrf_packet_info_t *packet) {
  struct timeval time;
  gettimeofday(&time, NULL);
  packet->time_secs_utc = time.tv_sec;
  packet->time_nano_secs = time.tv_usec * 1000;
}

//////////////////////////////////////////////////////////////////
// check packet id for validity, swapping as required.
// returns 0 on success, -1 on failure

int iwrf_check_packet_id(si32 packetId)

{

  si32 id = packetId;
  if (id >= 0x77770001 && id < 0x777700ff) {
    // no swapping needed
  } else {
    // swap
    SWAP_array_32(&id, sizeof(si32));
  }

  switch (id) {
    case IWRF_SYNC_ID:
    case IWRF_RADAR_INFO_ID:
    case IWRF_SCAN_SEGMENT_ID:
    case IWRF_ANTENNA_CORRECTION_ID:
    case IWRF_TS_PROCESSING_ID:
    case IWRF_XMIT_POWER_ID:
    case IWRF_XMIT_SAMPLE_ID:
    case IWRF_CALIBRATION_ID:
    case IWRF_EVENT_NOTICE_ID:
    case IWRF_PHASECODE_ID:
    case IWRF_XMIT_INFO_ID:
    case IWRF_PULSE_HEADER_ID:
    case IWRF_RVP8_PULSE_HEADER_ID:
    case IWRF_RVP8_OPS_INFO_ID:
      return 0;
  }

  return -1;

}

//////////////////////////////////////////////////////////////////
// check packet id for validity, swapping in-place as required.
// also swaps the packet_len argument.
// returns 0 on success, -1 on failure

int iwrf_check_packet_id_swp(si32 *packetId, si32 *packetLen)

{

  if (*packetId >= 0x77770001 && *packetId < 0x777700ff) {
    // no swapping needed
  } else {
    // swap
    SWAP_array_32(packetId, sizeof(si32));
    SWAP_array_32(packetLen, sizeof(si32));
  }

  return iwrf_check_packet_id(*packetId);

}

//////////////////////////////////////////////////////////////////
// get packet id, check for validity of this packet
// checks the packet length
// prints error in debug mode
// returns 0 on success, -1 on failure

int iwrf_get_packet_id(const void* buf, int len, int *packet_id)

{

  if (len < (int) sizeof(si32)) {
    return -1;
  }

  // get packet ID
  
  si32 id;
  memcpy(&id, buf, sizeof(si32));
  if (id >= 0x77770001 && id < 0x777700ff) {
    // no swapping needed
  } else {
    // swap
    printf("Warning swapping needed\n");
    SWAP_array_32(&id, sizeof(si32));
  }

  *packet_id = id;

  int iret = 0;
  switch (id) {

    case IWRF_SYNC_ID:
      if (len < (int) sizeof(iwrf_sync_t)) {
	iret = -1;
      } break;

    case IWRF_RADAR_INFO_ID:
      if (len < (int) sizeof(iwrf_radar_info_t)) {
	iret = -1;
      } break;

    case IWRF_SCAN_SEGMENT_ID:
      if (len < (int) sizeof(iwrf_scan_segment_t)) {
	iret = -1;
      } break;

    case IWRF_ANTENNA_CORRECTION_ID:
      if (len < (int) sizeof(iwrf_antenna_correction_t)) {
	iret = -1;
      } break;

    case IWRF_TS_PROCESSING_ID:
      if (len < (int) sizeof(iwrf_ts_processing_t)) {
	iret = -1;
      } break;

    case IWRF_XMIT_POWER_ID:
      if (len < (int) sizeof(iwrf_xmit_power_t)) {
	iret = -1;
      } break;

    case IWRF_XMIT_SAMPLE_ID:
      if (len < (int) sizeof(iwrf_xmit_sample_t)) {
	iret = -1;
      } break;

    case IWRF_CALIBRATION_ID:
      if (len < (int) sizeof(iwrf_calibration_t)) {
	iret = -1;
      } break;

    case IWRF_EVENT_NOTICE_ID:
      if (len < (int) sizeof(iwrf_event_notice_t)) {
	iret = -1;
      } break;

    case IWRF_PHASECODE_ID:
      if (len < (int) sizeof(iwrf_phasecode_t)) {
	iret = -1;
      } break;

    case IWRF_XMIT_INFO_ID:
      if (len < (int) sizeof(iwrf_xmit_info_t)) {
	iret = -1;
      } break;

    case IWRF_PULSE_HEADER_ID:
      if (len < (int) sizeof(iwrf_pulse_header_t)) {
	iret = -1;
      } break;

    case IWRF_RVP8_PULSE_HEADER_ID:
      if (len < (int) sizeof(iwrf_rvp8_pulse_header_t)) {
	iret = -1;
      } break;

    case IWRF_RVP8_OPS_INFO_ID:
      if (len < (int) sizeof(iwrf_rvp8_ops_info_t)) {
	iret = -1;
      } break;
  default:
      printf("unknown id %d\n",  id);
      iret = -1;
      break;
  }

  return iret;

}

//////////////////////////////////////////////////////////////////
// get packet time as a double

double iwrf_get_packet_time_as_double(const iwrf_packet_info_t *packet)

{
  return (packet->time_secs_utc + packet->time_nano_secs / 1.0e9);
}

//////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////
// swapping routines
//
// swap to native as required
// swapping is the responsibility of the user, data is always
// written in native

////////////////////////////////
// swap depending on packet type
// returns 0 on success, -1 on failure

int iwrf_packet_swap(void *buf, int len)
{

  int packet_id;
  if (iwrf_get_packet_id(buf, len, &packet_id)) {
    return -1;
  }

  switch (packet_id) {

    case IWRF_SYNC_ID:
      iwrf_sync_swap((iwrf_sync_t *) buf);
      break;

    case IWRF_RADAR_INFO_ID:
      iwrf_radar_info_swap((iwrf_radar_info_t *) buf);
      break;

    case IWRF_SCAN_SEGMENT_ID:
      iwrf_scan_segment_swap((iwrf_scan_segment_t *) buf);
      break;

    case IWRF_ANTENNA_CORRECTION_ID:
      iwrf_antenna_correction_swap((iwrf_antenna_correction_t *) buf);
      break;

    case IWRF_TS_PROCESSING_ID:
      iwrf_ts_processing_swap((iwrf_ts_processing_t *) buf);
      break;

    case IWRF_XMIT_POWER_ID:
      iwrf_xmit_power_swap((iwrf_xmit_power_t *) buf);
      break;

    case IWRF_XMIT_SAMPLE_ID:
      iwrf_xmit_sample_swap((iwrf_xmit_sample_t *) buf);
      break;

    case IWRF_CALIBRATION_ID:
      iwrf_calibration_swap((iwrf_calibration_t *) buf);
      break;

    case IWRF_EVENT_NOTICE_ID:
      iwrf_event_notice_swap((iwrf_event_notice_t *) buf);
      break;

    case IWRF_PHASECODE_ID:
      iwrf_phasecode_swap((iwrf_phasecode_t *) buf);
      break;

    case IWRF_XMIT_INFO_ID:
      iwrf_xmit_info_swap((iwrf_xmit_info_t *) buf);
      break;

    case IWRF_PULSE_HEADER_ID:
      iwrf_pulse_header_swap((iwrf_pulse_header_t *) buf);
      break;

    case IWRF_RVP8_PULSE_HEADER_ID:
      iwrf_rvp8_pulse_header_swap((iwrf_rvp8_pulse_header_t *) buf);
      break;

    case IWRF_RVP8_OPS_INFO_ID:
      iwrf_rvp8_ops_info_swap((iwrf_rvp8_ops_info_t *) buf);
      break;
      
  }

  return 0;

}

//////////////////////////////////////////////////////
// swap packet header
// returns true is swapped, false if already in native

int iwrf_packet_info_swap(iwrf_packet_info_t *packet)

{
  si32 id = packet->id;
  if (id >= 0x77770001 && id < 0x777700ff) {
    // no swapping needed
    return false;
  }
  SWAP_array_32(&packet->id, 4 * sizeof(si32));
  SWAP_array_64(&packet->seq_num, sizeof(si64));
  SWAP_array_64(&packet->time_secs_utc, sizeof(si64));
  SWAP_array_32(&packet->time_nano_secs, 4 * sizeof(si32));
  return true;
}

//////////////////////////////////////////////////////
// swap sync
// returns true is swapped, false if already in native

int iwrf_sync_swap(iwrf_sync_t *sync)

{
  // only swap the header
  // no data to be swapped since all bytes are identical
  return iwrf_packet_info_swap(&sync->packet);
}

//////////////////////////////////////////////////////
// swap radar_info
// returns true is swapped, false if already in native

int iwrf_radar_info_swap(iwrf_radar_info_t *radar_info)

{
  int swap = iwrf_packet_info_swap(&radar_info->packet);
  if (swap) {
    ui08 *start = (ui08 *) radar_info + sizeof(iwrf_packet_info_t);
    ui08 *end = (ui08 *) &radar_info->radar_name;
    int nbytes = end - start;
    SWAP_array_32(start, nbytes);
  }
  return swap;
}

//////////////////////////////////////////////////////
// swap scan_segment
// returns true is swapped, false if already in native

int iwrf_scan_segment_swap(iwrf_scan_segment_t *segment)

{
  int swap = iwrf_packet_info_swap(&segment->packet);
  if (swap) {
    ui08 *start = (ui08 *) segment + sizeof(iwrf_packet_info_t);
    ui08 *end = (ui08 *) &segment->segment_name;
    int nbytes = end - start;
    SWAP_array_32(start, nbytes);
  }
  return swap;
}

//////////////////////////////////////////////////////
// swap antenna_correction
// returns true is swapped, false if already in native

int iwrf_antenna_correction_swap(iwrf_antenna_correction_t *correction)

{
  int swap = iwrf_packet_info_swap(&correction->packet);
  if (swap) {
    ui08 *start = (ui08 *) correction + sizeof(iwrf_packet_info_t);
    int nbytes = sizeof(iwrf_antenna_correction_t) - sizeof(iwrf_packet_info_t);
    SWAP_array_32(start, nbytes);
  }
  return swap;
}

//////////////////////////////////////////////////////
// swap ts_processing
// returns true is swapped, false if already in native

int iwrf_ts_processing_swap(iwrf_ts_processing_t *processing)

{
  int swap = iwrf_packet_info_swap(&processing->packet);
  if (swap) {
    ui08 *start = (ui08 *) processing + sizeof(iwrf_packet_info_t);
    int nbytes = sizeof(iwrf_ts_processing_t) - sizeof(iwrf_packet_info_t);
    SWAP_array_32(start, nbytes);
  }
  return swap;
}

//////////////////////////////////////////////////////
// swap xmit_power
// returns true is swapped, false if already in native

int iwrf_xmit_power_swap(iwrf_xmit_power_t *power)

{
  int swap = iwrf_packet_info_swap(&power->packet);
  if (swap) {
    ui08 *start = (ui08 *) power + sizeof(iwrf_packet_info_t);
    int nbytes = sizeof(iwrf_xmit_power_t) - sizeof(iwrf_packet_info_t);
    SWAP_array_32(start, nbytes);
  }
  return swap;
}

//////////////////////////////////////////////////////
// swap xmit_sample
// returns true is swapped, false if already in native

int iwrf_xmit_sample_swap(iwrf_xmit_sample_t *sample)

{
  int swap = iwrf_packet_info_swap(&sample->packet);
  if (swap) {
    ui08 *start = (ui08 *) sample + sizeof(iwrf_packet_info_t);
    int nbytes = sizeof(iwrf_xmit_sample_t) - sizeof(iwrf_packet_info_t);
    SWAP_array_32(start, nbytes);
  }
  return swap;
}

//////////////////////////////////////////////////////
// swap xmit_info
// returns true is swapped, false if already in native

int iwrf_xmit_info_swap(iwrf_xmit_info_t *info)

{
  int swap = iwrf_packet_info_swap(&info->packet);
  if (swap) {
    ui08 *start = (ui08 *) info + sizeof(iwrf_packet_info_t);
    int nbytes = sizeof(iwrf_xmit_info_t) - sizeof(iwrf_packet_info_t);
    SWAP_array_32(start, nbytes);
  }
  return swap;
}

//////////////////////////////////////////////////////
// swap calibration
// returns true is swapped, false if already in native

int iwrf_calibration_swap(iwrf_calibration_t *calib)

{
  int swap = iwrf_packet_info_swap(&calib->packet);
  if (swap) {
    ui08 *start = (ui08 *) &calib + sizeof(iwrf_packet_info_t);
    int nbytes = sizeof(iwrf_calibration_t) - sizeof(iwrf_packet_info_t);
    SWAP_array_32(start, nbytes);
  }
  return swap;
}

//////////////////////////////////////////////////////
// swap event_notice
// returns true is swapped, false if already in native

int iwrf_event_notice_swap(iwrf_event_notice_t *notice)

{
  int swap = iwrf_packet_info_swap(&notice->packet);
  if (swap) {
    ui08 *start = (ui08 *) notice + sizeof(iwrf_packet_info_t);
    int nbytes = sizeof(iwrf_event_notice_t) - sizeof(iwrf_packet_info_t);
    SWAP_array_32(start, nbytes);
  }
  return swap;
}

//////////////////////////////////////////////////////
// swap phasecode
// returns true is swapped, false if already in native

int iwrf_phasecode_swap(iwrf_phasecode_t *code)

{
  int swap = iwrf_packet_info_swap(&code->packet);
  if (swap) {
    ui08 *start = (ui08 *) code + sizeof(iwrf_packet_info_t);
    int nbytes = sizeof(iwrf_phasecode_t) - sizeof(iwrf_packet_info_t);
    SWAP_array_32(start, nbytes);
  }
  return swap;
}

//////////////////////////////////////////////////////
// swap pulse_header
// returns true is swapped, false if already in native

int iwrf_pulse_header_swap(iwrf_pulse_header_t *pulse)

{
  int swap = iwrf_packet_info_swap(&pulse->packet);
  if (swap) {
    SWAP_array_64(&pulse->pulse_seq_num, sizeof(si64));
    int nn = sizeof(iwrf_packet_info_t) + sizeof(si64);
    ui08 *start = (ui08 *) pulse + nn;
    int nbytes = sizeof(iwrf_pulse_header_t) - nn;
    SWAP_array_32(start, nbytes);
  }
  return swap;
}

//////////////////////////////////////////////////////
// swap rvp8_pulse_header
// returns true is swapped, false if already in native

int iwrf_rvp8_pulse_header_swap(iwrf_rvp8_pulse_header_t *pulse)

{
  int swap = iwrf_packet_info_swap(&pulse->packet);

  if (swap) {

    pulse->i_tx_phase = SWAP_ui16(pulse->i_tx_phase);
    pulse->i_az = SWAP_ui16(pulse->i_az);
    pulse->i_el = SWAP_ui16(pulse->i_el);
    pulse->i_num_vecs = SWAP_si16(pulse->i_num_vecs);
    pulse->i_max_vecs = SWAP_si16(pulse->i_max_vecs);
    pulse->i_tg_wave = SWAP_ui16(pulse->i_tg_wave);

    pulse->i_btime_api = SWAP_ui32(pulse->i_btime_api);
    pulse->i_sys_time = SWAP_ui32(pulse->i_sys_time);
    pulse->i_prev_prt = SWAP_ui32(pulse->i_prev_prt);
    pulse->i_next_prt = SWAP_ui32(pulse->i_next_prt);
    pulse->i_seq_num = SWAP_ui32(pulse->i_seq_num);

    SWAP_array_32(&pulse->uiq_perm, 2 * sizeof(ui32));
    SWAP_array_32(&pulse->uiq_once, 2 * sizeof(ui32));

    SWAP_array_32(&pulse->i_data_off, IWRF_MAX_CHAN * sizeof(si32));
    SWAP_array_32(&pulse->f_burst_mag, IWRF_MAX_CHAN * sizeof(fl32));
    SWAP_array_16(&pulse->i_burst_arg, IWRF_MAX_CHAN * sizeof(si16));
    SWAP_array_16(&pulse->i_wrap_iq, IWRF_MAX_CHAN * sizeof(si16));

    SWAP_array_32(&pulse->unused2, 29 * sizeof(si32));

  }
  return swap;
}

//////////////////////////////////////////////////////
// swap rvp8_ops_info
// returns true is swapped, false if already in native

int iwrf_rvp8_ops_info_swap(iwrf_rvp8_ops_info_t *info)

{
  int swap = iwrf_packet_info_swap(&info->packet);

  if (swap) {

    info->i_version = SWAP_si32(info->i_version);


    info->i_major_mode = SWAP_ui32(info->i_major_mode);
    info->i_polarization = SWAP_ui32(info->i_polarization);
    info->i_phase_mode_seq = SWAP_ui32(info->i_phase_mode_seq);

    info->i_task_sweep = SWAP_ui16(info->i_task_sweep);
    info->i_task_aux_num = SWAP_ui16(info->i_task_aux_num);

    info->i_task_scan_type = SWAP_si32(info->i_task_scan_type);
    SWAP_array_32(&info->unused1, 3 * sizeof(si32));

    info->i_aq_mode = SWAP_ui32(info->i_aq_mode);
    info->i_unfold_mode = SWAP_ui32(info->i_unfold_mode);

    info->i_pwidth_code = SWAP_ui32(info->i_pwidth_code);
    info->f_pwidth_usec = SWAP_fl32(info->f_pwidth_usec);

    info->f_dbz_calib = SWAP_fl32(info->f_dbz_calib);

    info->i_sample_size = SWAP_si32(info->i_sample_size);

    info->i_mean_angle_sync = SWAP_ui32(info->i_mean_angle_sync);
    info->i_flags = SWAP_ui32(info->i_flags);

    info->i_playback_version = SWAP_si32(info->i_playback_version);

    info->f_sy_clk_mhz = SWAP_fl32(info->f_sy_clk_mhz);
    info->f_wavelength_cm = SWAP_fl32(info->f_wavelength_cm);
    info->f_saturation_dbm = SWAP_fl32(info->f_saturation_dbm);
    info->f_range_mask_res = SWAP_fl32(info->f_range_mask_res);

    SWAP_array_16(&info->i_range_mask, IWRF_RVP8_GATE_MASK_LEN * sizeof(si16));

    SWAP_array_32(&info->f_noise_dbm, IWRF_MAX_CHAN * sizeof(fl32));
    SWAP_array_32(&info->f_noise_stdv_db, IWRF_MAX_CHAN * sizeof(fl32));

    info->f_noise_range_km = SWAP_fl32(info->f_noise_range_km);
    info->f_noise_prf_hz = SWAP_fl32(info->f_noise_prf_hz);

    SWAP_array_16(&info->i_gparm_latch_sts, 2 * sizeof(si16));
    SWAP_array_16(&info->i_gparm_immed_sts, 6 * sizeof(si16));
    SWAP_array_16(&info->i_gparm_diag_bits, 4 * sizeof(si16));

    SWAP_array_32(&info->unused2, 189 * sizeof(si32));

  }
  return swap;
}

//////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////
// string representation of enums

// string representation of packet_id 

char *iwrf_packet_id_to_str(int packet_id)

{
  
  switch (packet_id) {
    case IWRF_SYNC_ID:
      return "IWRF_SYNC_ID";
    case IWRF_RADAR_INFO_ID:
      return "IWRF_RADAR_INFO_ID";
    case IWRF_SCAN_SEGMENT_ID:
      return "IWRF_SCAN_SEGMENT_ID";
    case IWRF_ANTENNA_CORRECTION_ID:
      return "IWRF_ANTENNA_CORRECTION_ID";
    case IWRF_TS_PROCESSING_ID:
      return "IWRF_TS_PROCESSING_ID";
    case IWRF_XMIT_POWER_ID:
      return "IWRF_XMIT_POWER_ID";
    case IWRF_XMIT_SAMPLE_ID:
      return "IWRF_XMIT_SAMPLE_ID";
    case IWRF_CALIBRATION_ID:
      return "IWRF_CALIBRATION_ID";
    case IWRF_EVENT_NOTICE_ID:
      return "IWRF_EVENT_NOTICE_ID";
    case IWRF_PHASECODE_ID:
      return "IWRF_PHASECODE_ID";
    case IWRF_XMIT_INFO_ID:
      return "IWRF_XMIT_INFO_ID";
    case IWRF_PULSE_HEADER_ID:
      return "IWRF_PULSE_HEADER_ID";
    case IWRF_RVP8_PULSE_HEADER_ID:
      return "IWRF_RVP8_PULSE_HEADER_ID";
    case IWRF_RVP8_OPS_INFO_ID:
      return "IWRF_RVP8_OPS_INFO_ID";
    default:
      return "IWRF_UNKNOWN_ID";
  }

}

// string representation of xmit_rcv_mode

char *iwrf_xmit_rcv_mode_to_str(int xmit_rcv_mode)

{
  
  switch (xmit_rcv_mode) {
    case IWRF_SINGLE_POL:
      return "IWRF_SINGLE_POL";
    case IWRF_ALT_HV_CO_ONLY:
      return "IWRF_ALT_HV_CO_ONLY";
    case IWRF_ALT_HV_CO_CROSS:
      return "IWRF_ALT_HV_CO_CROSS";
    case IWRF_ALT_HV_FIXED_HV:
      return "IWRF_ALT_HV_FIXED_HV";
    case IWRF_SIM_HV_FIXED_HV:
      return "IWRF_SIM_HV_FIXED_HV";
    case IWRF_SIM_HV_SWITCHED_HV:
      return "IWRF_SIM_HV_SWITCHED_HV";
    case IWRF_H_ONLY_FIXED_HV:
      return "IWRF_H_ONLY_FIXED_HV";
    case IWRF_V_ONLY_FIXED_HV:
      return "IWRF_V_ONLY_FIXED_HV";
    default:
      return "IWRF_XMIT_RCV_MODE unknown";
  }

}

// string representation of iwrf_xmit_phase_mode

char *iwrf_xmit_phase_mode_to_str(int xmit_phase_mode)

{
  
  switch (xmit_phase_mode) {
    case IWRF_XMIT_PHASE_MODE_FIXED:
      return "IWRF_XMIT_PHASE_MODE_FIXED";
    case IWRF_XMIT_PHASE_MODE_RANDOM:
      return "IWRF_XMIT_PHASE_MODE_RANDOM";
    case IWRF_XMIT_PHASE_MODE_SZ864:
      return "IWRF_XMIT_PHASE_MODE_SZ864";
    default:
      return "UNKNOWN";
  }

}

// string representation of prf_mode

char *iwrf_prf_mode_to_str(int prf_mode)

{
  
  switch (prf_mode) {
    case IWRF_PRF_MODE_FIXED:
      return "IWRF_PRF_MODE_FIXED";
    case IWRF_PRF_MODE_STAGGERED_2_3:
      return "IWRF_PRF_MODE_STAGGERED_2_3";
    case IWRF_PRF_MODE_STAGGERED_3_4:
      return "IWRF_PRF_MODE_STAGGERED_3_4";
    case IWRF_PRF_MODE_STAGGERED_4_5:
      return "IWRF_PRF_MODE_STAGGERED_4_5";
    case IWRF_PRF_MODE_MULTI_PRT:
      return "IWRF_PRF_MODE_MULTI_PRT";
    case IWRF_PRF_MODE_BLOCK_MODE:
      return "IWRF_PRF_MODE_BLOCK_MODE";
    default:
      return "UNKNOWN";
  }

}

// Scale float pulse data

void iwrf_pulse_scale_data(void * buf,int len, double scale, double bias)
  
{
   iwrf_pulse_header_t *pHdr = (iwrf_pulse_header_t *) buf;
   ui08 *uptr = (ui08 *) buf;  // for pointer math
   float *fptr = (float *) (uptr + sizeof(iwrf_pulse_header_t));
   int i=0;
   for(i=0; i < pHdr->n_data; i++) fptr[i] = fptr[i] * scale + bias;
}


// string representation of pulse_type

char *iwrf_pulse_type_to_str(int pulse_type)
  
{
  
  switch (pulse_type) {
    case IWRF_PULSE_TYPE_RECT:
      return "IWRF_PULSE_TYPE_RECT";
    case IWRF_PULSE_TYPE_GAUSSIAN:
      return "IWRF_PULSE_TYPE_GAUSSIAN";
    default:
      return "UNKNOWN";
  }

}

// string representation of pulse_polarization

char *iwrf_pol_mode_to_str(int pol_mode)

{
  
  switch (pol_mode) {
    case IWRF_POL_MODE_H:
      return "IWRF_POL_MODE_H";
    case IWRF_POL_MODE_V:
      return "IWRF_POL_MODE_V";
    case IWRF_POL_MODE_HV_ALT:
      return "IWRF_POL_MODE_HV_ALT";
    case IWRF_POL_MODE_HV_SIM:
      return "IWRF_POL_MODE_HV_SIM";
    default:
      return "UNKNOWN";
  }

}

// string representation of scan_mode

char *iwrf_scan_mode_to_str(int scan_mode)

{
  
  switch (scan_mode) {
    case IWRF_SCAN_MODE_SECTOR:
      return "IWRF_SCAN_MODE_SECTOR";
    case IWRF_SCAN_MODE_COPLANE:
      return "IWRF_SCAN_MODE_COPLANE";
    case IWRF_SCAN_MODE_RHI:
      return "IWRF_SCAN_MODE_RHI";
    case IWRF_SCAN_MODE_VERTICAL_POINTING:
      return "IWRF_SCAN_MODE_VERTICAL_POINTING";
    case IWRF_SCAN_MODE_IDLE:
      return "IWRF_SCAN_MODE_IDLE";
    case IWRF_SCAN_MODE_AZ_SUR_360:
      return "IWRF_SCAN_MODE_AZ_SUR_360";
    case IWRF_SCAN_MODE_EL_SUR_360:
      return "IWRF_SCAN_MODE_EL_SUR_360";
    case IWRF_SCAN_MODE_SUNSCAN:
      return "IWRF_SCAN_MODE_SUNSCAN";
    case IWRF_SCAN_MODE_POINTING:
      return "IWRF_SCAN_MODE_POINTING";
    case IWRF_SCAN_MODE_MANPPI:
      return "IWRF_SCAN_MODE_MANPPI";
    case IWRF_SCAN_MODE_MANRHI:
      return "IWRF_SCAN_MODE_MANRHI";
    default:
      return "UNKNOWN";
  }

}

// string representation of follow_mode

char *iwrf_follow_mode_to_str(int follow_mode)

{
  
  switch (follow_mode) {
    case IWRF_FOLLOW_MODE_NONE:
      return "IWRF_FOLLOW_MODE_NONE";
    case IWRF_FOLLOW_MODE_SUN:
      return "IWRF_FOLLOW_MODE_SUN";
    case IWRF_FOLLOW_MODE_VEHICLE:
      return "IWRF_FOLLOW_MODE_VEHICLE";
    case IWRF_FOLLOW_MODE_AIRCRAFT:
      return "IWRF_FOLLOW_MODE_AIRCRAFT";
    case IWRF_FOLLOW_MODE_TARGET:
      return "IWRF_FOLLOW_MODE_TARGET";
    case IWRF_FOLLOW_MODE_MANUAL:
      return "IWRF_FOLLOW_MODE_MANUAL";
    default:
      return "UNKNOWN";
  }

}

// string representation of radar_platform

char *iwrf_radar_platform_to_str(int radar_platform)

{
  
  switch (radar_platform) {
    case IWRF_RADAR_PLATFORM_FIXED:
      return "IWRF_RADAR_PLATFORM_FIXED";
    case IWRF_RADAR_PLATFORM_VEHICLE:
      return "IWRF_RADAR_PLATFORM_VEHICLE";
    case IWRF_RADAR_PLATFORM_SHIP:
      return "IWRF_RADAR_PLATFORM_SHIP";
    case IWRF_RADAR_PLATFORM_AIRCRAFT:
      return "IWRF_RADAR_PLATFORM_AIRCRAFT";
    default:
      return "UNKNOWN";
  }

}

// string representation of cal_type
  
char *iwrf_cal_type_to_str(int cal_type)

{
  
  switch (cal_type) {
    case IWRF_CAL_TYPE_NONE:
      return "IWRF_CAL_TYPE_NONE";
    case IWRF_CAL_TYPE_CW_CAL:
      return "IWRF_CAL_TYPE_CW_CAL";
    case IWRF_CAL_TYPE_SOLAR_CAL_FIXED:
      return "IWRF_CAL_TYPE_SOLAR_CAL_FIXED";
    case IWRF_CAL_TYPE_SOLAR_CAL_SCAN:
      return "IWRF_CAL_TYPE_SOLAR_CAL_SCAN";
    case IWRF_CAL_TYPE_NOISE_SOURCE_H:
      return "IWRF_CAL_TYPE_NOISE_SOURCE_H";
    case IWRF_CAL_TYPE_NOISE_SOURCE_V:
      return "IWRF_CAL_TYPE_NOISE_SOURCE_V";
    case IWRF_CAL_TYPE_NOISE_SOURCE_HV:
      return "IWRF_CAL_TYPE_NOISE_SOURCE_HV";
    case IWRF_CAL_TYPE_BLUESKY:
      return "IWRF_CAL_TYPE_BLUESKY";
    case IWRF_CAL_TYPE_SAVEPARAMS:
      return "IWRF_CAL_TYPE_SAVEPARAMS";
    default:
      return "UNKNOWN";
  }

}

// string representation of event_cause

char *iwrf_event_cause_to_str(int event_cause)

{
  
  switch (event_cause) {
    case IWRF_EVENT_CAUSE_DONE:
      return "IWRF_EVENT_CAUSE_DONE";
    case IWRF_EVENT_CAUSE_TIMEOUT:
      return "IWRF_EVENT_CAUSE_TIMEOUT";
    case IWRF_EVENT_CAUSE_TIMER:
      return "IWRF_EVENT_CAUSE_TIMER";
    case IWRF_EVENT_CAUSE_ABORT:
      return "IWRF_EVENT_CAUSE_ABORT";
    case IWRF_EVENT_CAUSE_SCAN_ABORT:
      return "IWRF_EVENT_CAUSE_SCAN_ABORT";
    case IWRF_EVENT_CAUSE_RESTART:
      return "IWRF_EVENT_CAUSE_RESTART";
    case IWRF_EVENT_CAUSE_SCAN_STATE_TIMEOUT:
      return "IWRF_EVENT_CAUSE_SCAN_STATE_TIMEOUT";
    case IWRF_EVENT_CAUSE_ANTENNA_FAULT:
      return "IWRF_EVENT_CAUSE_ANTENNA_FAULT";
    default:
      return "UNKNOWN";
  }

}

// string representation of iq_encoding/

char *iwrf_iq_encoding_to_str(int iq_encoding)

{
  
  switch (iq_encoding) {
    case IWRF_IQ_ENCODING_FL32:
      return "IWRF_IQ_ENCODING_FL32";
    case IWRF_IQ_ENCODING_SCALED_SI16:
      return "IWRF_IQ_ENCODING_SCALED_SI16";
    case IWRF_IQ_ENCODING_DBM_PHASE_SI16:
      return "IWRF_IQ_ENCODING_DBM_PHASE_SI16";
    case IWRF_IQ_ENCODING_SIGMET_FL16:
      return "IWRF_IQ_ENCODING_SIGMET_FL16";
    default:
      return "UNKNOWN";
  }

}

// string representation of (user interface) ui_error

char *iwrf_ui_error_to_str(int error)

{
  switch(error) {
    case IWRF_UI_ERROR_DELETE_FAILED:
      return "IWRF_UI_ERROR_DELETE_FAILED task was running";
    case IWRF_UI_ERROR_TASKLIST_SIZE:
        return "IWRF_UI_ERROR_TASKLIST_SIZE";
    case IWRF_UI_ERROR_MISSING_TASK_NAME:
        return "IWRF_UI_ERROR_MISSING_TASK_NAME";
    case IWRF_UI_ERROR_TASKLIST_INDEX_RANGE:
        return "IWRF_UI_ERROR_TASKLIST_INDEX_RANGE";
    case IWRF_UI_ERROR_TASKLIST_FULL:
        return "IWRF_UI_ERROR_TASKLIST_FULL";
    case IWRF_UI_ERROR_APPEND_TASK_UNDEFINED:
        return "IWRF_UI_ERROR_APPEND_TASK_UNDEFINED";
    case IWRF_UI_ERROR_APPEND_TASK_NONAME:
        return "IWRF_UI_ERROR_APPEND_TASK_NONAME";
    case IWRF_UI_ERROR_CANT_SCHEDULE:
        return "IWRF_UI_ERROR_CANT_SCHEDULE";
    case IWRF_UI_ERROR_REPEAT_CYCLE_RANGE:
        return "IWRF_UI_ERROR_REPEAT_CYCLE_RANGE";
    case IWRF_UI_ERROR_TASKLIST_ENTRY_UNKNOWN:
		return "IWRF_UI_ERROR_TASKLIST_ENTRY_UNKNOWN";
    case IWRF_UI_ERROR_TASKLIST_ENTRY_UNDEFINED:
		return "IWRF_UI_ERROR_TASKLIST_ENTRY_UNDEFINED";
    case IWRF_UI_ERROR_CANT_SCHEDULE_TASK_NOT_IN_TASKLIST:
		return "IWRF_UI_ERROR_CANT_SCHEDULE_TASK_NOT_IN_TASKLIST";
    case IWRF_ANTCON_NOT_CONNECTED:
		return "IWRF_ANTCON_NOT_CONNECTED";
    case IWRF_ANTENA_FAULTED:
		return "IWRF_ANTENA_FAULTED";
    case IWRF_TXCTRL_NOT_CONNECTED:
    	return "IWRF_TXCTRL_NOT_CONNECTED";
    case IWRF_SCAN_SEGMENT_LIST_ERROR:
		return "IWRF_SCAN_SEGMENT_LIST_ERROR";
    case IWRF_UI_ERROR_BAD_PACKET_ID:
		return "IWRF_UI_ERROR_BAD_PACKET_ID";
    case IWRF_UI_WARN_START_TIME_WARN:
		return "IWRF_UI_WARN_START_TIME_WARN";
	case IWRF_UI_DATA_SYSTEM_NOT_CONNECTED:
		return "IWRF_UI_DATA_SYSTEM_NOT_CONNECTED";
    default:
      return "UNKNOWN";
  }

}

//////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////
// printing routines

////////////////////////////////
// print depending on packet type

void iwrf_packet_print(FILE *out, const void *buf, int len)

{
  
  int packet_id;
  if (iwrf_get_packet_id(buf, len, &packet_id)) {
    return;
  }

  switch (packet_id) {

    case IWRF_SYNC_ID:
      iwrf_sync_print(out, (iwrf_sync_t *) buf);
      break;

    case IWRF_RADAR_INFO_ID:
      iwrf_radar_info_print(out, (iwrf_radar_info_t *) buf);
      break;

    case IWRF_SCAN_SEGMENT_ID:
      iwrf_scan_segment_print(out, (iwrf_scan_segment_t *) buf);
      break;

    case IWRF_ANTENNA_CORRECTION_ID:
      iwrf_antenna_correction_print(out, (iwrf_antenna_correction_t *) buf);
      break;

    case IWRF_TS_PROCESSING_ID:
      iwrf_ts_processing_print(out, (iwrf_ts_processing_t *) buf);
      break;

    case IWRF_XMIT_POWER_ID:
      iwrf_xmit_power_print(out, (iwrf_xmit_power_t *) buf);
      break;

    case IWRF_XMIT_SAMPLE_ID:
      iwrf_xmit_sample_print(out, (iwrf_xmit_sample_t *) buf);
      break;

    case IWRF_CALIBRATION_ID:
      iwrf_calibration_print(out, (iwrf_calibration_t *) buf);
      break;

    case IWRF_EVENT_NOTICE_ID:
      iwrf_event_notice_print(out, (iwrf_event_notice_t *) buf);
      break;

    case IWRF_PHASECODE_ID:
      iwrf_phasecode_print(out, (iwrf_phasecode_t *) buf);
      break;

    case IWRF_XMIT_INFO_ID:
      iwrf_xmit_info_print(out, (iwrf_xmit_info_t *) buf);
      break;

    case IWRF_PULSE_HEADER_ID:
      iwrf_pulse_header_print(out, (iwrf_pulse_header_t *) buf);
      break;

    case IWRF_RVP8_PULSE_HEADER_ID:
      iwrf_rvp8_pulse_header_print(out, (iwrf_rvp8_pulse_header_t *) buf);
      break;

    case IWRF_RVP8_OPS_INFO_ID:
      iwrf_rvp8_ops_info_print(out, (iwrf_rvp8_ops_info_t *) buf);
      break;
    default:
        fprintf(out, "unknown id %d\n", packet_id);
        break;
  }

}

char *iwrf_time_str(const si64 *ptime, si32 *nano_secs)
{
	static char str1[30]={"                             "};
	static char str2[50];
	const time_t *ptime2= (const time_t *)ptime;

	char *asc = asctime(gmtime(ptime2));
	if (asc) {
	    strncpy(str1,asc, 29);
	    strncpy(str2, str1+4, 6);	// Month Day	
	    strncpy(str2+6, str1+19, 5);	// Year
	    *(str2+11)= ' ';
	    strncpy(str2+12, str1+11, 8);	// HH:MM:SS
	    *(str2+20)= '.';
	    *(str2+21)= 0;
	    if(nano_secs) sprintf(str2+21, "%.9d", *nano_secs);	// fractions of sec
	} else {
	    strcpy(str2, "BAD_TIME");
        }
	return(&str2[0]);
}
void iwrf_time_str2(const si64 *ptime, si32 *nano_secs, char *str2)
{
	char str1[30]={"                             "};

	char *asc = asctime(gmtime( (time_t *)ptime));
	if (asc) {
	    strncpy(str1,asc, 29);
	    strncpy(str2, str1+4, 6);	// Month Day	
	    strncpy(str2+6, str1+19, 5);	// Year
	    *(str2+11)= ' ';
	    strncpy(str2+12, str1+11, 8);	// HH:MM:SS
	    *(str2+20)= '.';
	    *(str2+21)= 0;
	    if(nano_secs) sprintf(str2+21, "%.9d", *nano_secs);	// fractions of sec
	} else {
	    strcpy(str2, "BAD_TIME");
        }
}


//////////////////////////////////////////////////////
// print UI schedule header
void iwrf_ui_schedule_print(FILE *out,
			    const iwrf_ui_schedule_info_t *schedule)
{
	fprintf(out, "begin time: %s\n", iwrf_time_str(&schedule->begin_time_utc,NULL));
	fprintf(out, "repeat cycle %d secs\n", schedule->repeat_cycle_secs);
	fprintf(out, "priority: %d\n", schedule->priority);
	if(schedule->last_run_time_utc==0) fprintf(out, "last run at: never run\n");
	else fprintf(out, "last run at: %s\n", 
			iwrf_time_str(&schedule->last_run_time_utc,NULL));
}

//////////////////////////////////////////////////////
// print packet header

void iwrf_packet_info_print(FILE *out,
			    const iwrf_packet_info_t *packet)
{

  fprintf(out, "  id: 0x%x\n", packet->id);
  fprintf(out, "  len_bytes: %d\n", packet->len_bytes);
  fprintf(out, "  seq_num: %lld\n", (long long) packet->seq_num);
  fprintf(out, "  version_num: %d\n", packet->version_num);
  fprintf(out, "  radar_id: %d\n", packet->radar_id);
  
  fprintf(out, "  time UTC: %s\n", iwrf_time_str(&packet->time_secs_utc, 
					(si32 *)&packet->time_nano_secs));

}

//////////////////////////////////////////////////////
// print radar_info

void iwrf_radar_info_print(FILE *out,
			   const iwrf_radar_info_t *info)

{
  char *cp;

  fprintf(out, "==================== iwrf_radar_info ============================\n");
  iwrf_packet_info_print(out, &info->packet);

  fprintf(out, "  latitude_deg: %g\n", info->latitude_deg);
  fprintf(out, "  longitude_deg: %g\n", info->longitude_deg);
  fprintf(out, "  altitude_m: %g\n", info->altitude_m);
  fprintf(out, "  platform_type: %s\n",
	  iwrf_radar_platform_to_str(info->platform_type));
  fprintf(out, "  beamwidth_deg_h: %g\n", info->beamwidth_deg_h);
  fprintf(out, "  beamwidth_deg_v: %g\n", info->beamwidth_deg_v);
  fprintf(out, "  wavelength_cm: %g\n", info->wavelength_cm);
  fprintf(out, "  nominal_gain_ant_db_h: %g\n", info->nominal_gain_ant_db_h);
  fprintf(out, "  nominal_gain_ant_db_v: %g\n", info->nominal_gain_ant_db_v);
  cp= iwrf_safe_str(info->radar_name, IWRF_MAX_RADAR_NAME);
  fprintf(out, "  radar_name: %s\n", cp);
  free(cp);
  cp= iwrf_safe_str(info->site_name, IWRF_MAX_SITE_NAME);
  fprintf(out, "  site_name: %s\n", cp);
  free(cp);
  fprintf(out, "=================================================================\n");

}

//////////////////////////////////////////////////////
// print sync packet

void iwrf_sync_print(FILE *out,
		     const iwrf_sync_t *sync)

{

  fprintf(out, "==================== iwrf_sync ==================================\n");
  iwrf_packet_info_print(out, &sync->packet);
  fprintf(out, "  magik[0]: 0x%x\n", sync->magik[0]);
  fprintf(out, "  magik[1]: 0x%x\n", sync->magik[1]);
  fprintf(out, "=================================================================\n");

}

//////////////////////////////////////////////////////
// print version packet

void iwrf_version_print(FILE *out,
		     const iwrf_version_t *version)

{

  fprintf(out, "==================== iwrf_version ===============================\n");
  iwrf_packet_info_print(out, &version->packet);
  fprintf(out, "  major_version_num: %d\n", version->major_version_num);
  fprintf(out, "  minor_version_num: %d\n", version->minor_version_num);
  fprintf(out, "  version_name: %s\n", version->version_name);
  fprintf(out, "=================================================================\n");

}

//////////////////////////////////////////////////////
// print scan_segment

void iwrf_scan_segment_print(FILE *out,
			   const iwrf_scan_segment_t *seg)

{
  char *cp;

  fprintf(out, "==================== iwrf_scan_segment ==========================\n");
  iwrf_packet_info_print(out, &seg->packet);

  fprintf(out, "  scan_mode: %s\n", iwrf_scan_mode_to_str(seg->scan_mode));
  fprintf(out, "  follow_mode: %s\n", iwrf_follow_mode_to_str(seg->follow_mode));
  fprintf(out, "  volume_num: %d\n", seg->volume_num);
  fprintf(out, "  sweep_num: %d\n", seg->sweep_num);
  fprintf(out, "  time_limit: %d\n", seg->time_limit);
  fprintf(out, "  az_manual: %g\n", seg->az_manual);
  fprintf(out, "  el_manual: %g\n", seg->el_manual);
  fprintf(out, "  az_start: %g\n", seg->az_start);
  fprintf(out, "  el_start: %g\n", seg->el_start);
  fprintf(out, "  scan_rate: %g\n", seg->scan_rate);
  fprintf(out, "  left_limit: %g\n", seg->left_limit);
  fprintf(out, "  right_limit: %g\n", seg->right_limit);
  fprintf(out, "  up_limit: %g\n", seg->up_limit);
  fprintf(out, "  down_limit: %g\n", seg->down_limit);
  fprintf(out, "  step: %g\n", seg->step);
  fprintf(out, "  current_fixed_angle: %g\n", seg->current_fixed_angle);
  fprintf(out, "  init_direction_cw: %d\n", seg->init_direction_cw);
  fprintf(out, "  init_direction_up: %d\n", seg->init_direction_up);
  fprintf(out, "  n_sweeps: %d\n", seg->n_sweeps);

  fprintf(out, "  fixed_angles:");
  int ii = 0;
  for (; ii < seg->n_sweeps; ii++) {
    fprintf(out, " %g", seg->fixed_angles[ii]);
  }
  fprintf(out, "\n");

  cp= iwrf_safe_str(seg->segment_name, IWRF_MAX_SEGMENT_NAME);
  fprintf(out, "  segment_name: %s\n", cp);
  free(cp);
  cp= iwrf_safe_str(seg->project_name, IWRF_MAX_PROJECT_NAME);
  fprintf(out, "  project_name: %s\n", cp);
  free(cp);
  
  fprintf(out, "=================================================================\n");

}

//////////////////////////////////////////////////////
// print antenna_correction

void iwrf_antenna_correction_print(FILE *out,
				   const iwrf_antenna_correction_t *corr)

{

  fprintf(out, "==================== iwrf_antenna_correction ====================\n");
  iwrf_packet_info_print(out, &corr->packet);
  
  fprintf(out, "  az_correction: %g\n", corr->az_correction);
  fprintf(out, "  el_correction: %g\n", corr->el_correction);
  fprintf(out, "=================================================================\n");

}

//////////////////////////////////////////////////////
// print ts_processing

void iwrf_ts_processing_print(FILE *out,
			      const iwrf_ts_processing_t *proc)

{

  fprintf(out, "==================== iwrf_ts_processing =========================\n");
  iwrf_packet_info_print(out, &proc->packet);

  fprintf(out, "  xmit_rcv_mode: %s\n",
	  iwrf_xmit_rcv_mode_to_str(proc->xmit_rcv_mode));
  fprintf(out, "  pol_mode: %s\n", 
	  iwrf_pol_mode_to_str(proc->pol_mode));
  fprintf(out, "  xmit_phase_mode: %s\n",
	  iwrf_xmit_phase_mode_to_str(proc->xmit_phase_mode));
  fprintf(out, "  pulse_type: %s\n",
          iwrf_pulse_type_to_str(proc->pulse_type));

  fprintf(out, "  prf_mode: %s\n",
          iwrf_prf_mode_to_str(proc->prf_mode));
  fprintf(out, "  num_prts: %d\n", proc->num_prts);
  fprintf(out, "  prt_usec: %g\n", proc->prt_usec);
  fprintf(out, "  prt2_usec: %g\n", proc->prt2_usec);
  fprintf(out, "  prt3_usec: %g\n", proc->prt3_usec);
  fprintf(out, "  prt4_usec: %g\n", proc->prt4_usec);
  fprintf(out, "  integration_cycle_pulses: %d\n",
          proc->integration_cycle_pulses);
  fprintf(out, "  block_mode_prt2_pulses: %d\n", proc->block_mode_prt2_pulses);
  fprintf(out, "  block_mode_prt3_pulses: %d\n", proc->block_mode_prt3_pulses);
  fprintf(out, "  block_mode_prt4_pulses: %d\n", proc->block_mode_prt4_pulses);

  fprintf(out, "  cal_type: %s\n",
          iwrf_cal_type_to_str(proc->cal_type));
  
  fprintf(out, "  burst_range_offset_m: %g\n", proc->burst_range_offset_m);
  fprintf(out, "  pulse_width_us: %g\n", proc->pulse_width_us);
  fprintf(out, "  start_range_m: %g\n", proc->start_range_m);
  fprintf(out, "  gate_spacing_m: %g\n", proc->gate_spacing_m);

  fprintf(out, "  indexed_beam_spacing: %.2f  beams are indexed=%d\n", 
	proc->indexed_beam_spacing_deg, proc->beams_are_indexed);
  fprintf(out, "  indexed_beam_width: %.2f  use indexed_beam_width=%d\n", 
	proc->indexed_beam_width_deg, proc->specify_dwell_width);

  fprintf(out, "  clutter_filter_number: %d\n", proc->clutter_filter_number);
  fprintf(out, "  range_gate_averaging: %d\n", proc->range_gate_averaging);
  fprintf(out, "  max_gate: %d\n", proc->max_gate);

  fprintf(out, "  test_power_dbm: %g\n", proc->test_power_dbm);

  fprintf(out, "  test_pulse_range_km: %g\n", proc->test_pulse_range_km);
  fprintf(out, "  test_pulse_length_usec: %g\n", proc->test_pulse_length_usec);
			   
}

//////////////////////////////////////////////////////
// print xmit_power

void iwrf_xmit_power_print(FILE *out,
			   const iwrf_xmit_power_t *pwr)

{

  fprintf(out, "==================== iwrf_xmit_power ============================\n");
  iwrf_packet_info_print(out, &pwr->packet);
  
  fprintf(out, "  power_dbm_h: %g\n", pwr->power_dbm_h);
  fprintf(out, "  power_dbm_v: %g\n", pwr->power_dbm_v);
  fprintf(out, "=================================================================\n");

}

//////////////////////////////////////////////////////
// print xmit_sample

void iwrf_xmit_sample_print(FILE *out,
			    const iwrf_xmit_sample_t *samp)

{

  fprintf(out, "==================== iwrf_xmit_sample ===========================\n");
  iwrf_packet_info_print(out, &samp->packet);
  
  fprintf(out, "  power_dbm_h: %g\n", samp->power_dbm_h);
  fprintf(out, "  power_dbm_v: %g\n", samp->power_dbm_v);
  fprintf(out, "  offset: %d\n", samp->offset);
  fprintf(out, "  n_samples: %d\n", samp->n_samples);
  fprintf(out, "  sampling_freq: %g\n", samp->sampling_freq);
  fprintf(out, "  scale_h: %g\n", samp->scale_h);
  fprintf(out, "  offset_h: %g\n", samp->offset_h);
  fprintf(out, "  scale_v: %g\n", samp->scale_v);
  fprintf(out, "  offset_v: %g\n", samp->offset_v);

  int n_samples = samp->n_samples;
  if (n_samples > IWRF_N_TXSAMP) {
    n_samples = IWRF_N_TXSAMP;
  }
 
  fprintf(out, "  samples_h:\n");
  int ii = 0;
  for (; ii < n_samples; ii++) {
    fprintf(out, "    index, val: %d, %d\n", ii,
            samp->samples_h[ii]);
  }

  fprintf(out, "  samples_v:\n");
  for (ii = 0; ii < n_samples; ii++) {
    fprintf(out, "    index, val: %d, %d\n", ii,
            samp->samples_v[ii]);
  }
  fprintf(out, "=================================================================\n");

}

//////////////////////////////////////////////////////
// print xmit_info

void iwrf_xmit_info_print(FILE *out,
			  const iwrf_xmit_info_t *info)

{
  
  fprintf(out, "==================== iwrf_xmit_info =============================\n");
  iwrf_packet_info_print(out, &info->packet);
  
  fprintf(out, "  xmit_0_enabled: %d\n", info->xmit_0_enabled);
  fprintf(out, "  xmit_1_enabled: %d\n", info->xmit_1_enabled);
  fprintf(out, "  xmit_rcv_mode: %d\n", info->xmit_rcv_mode);
  fprintf(out, "  xmit_phase_mode: %s\n",
	  iwrf_xmit_phase_mode_to_str(info->xmit_phase_mode));
  fprintf(out, "  prf_mode: %s\n", iwrf_prf_mode_to_str(info->prf_mode));
  fprintf(out, "  pulse_type: %s\n", iwrf_pulse_type_to_str(info->pulse_type));
  fprintf(out, "  prt_usec: %g\n", info->prt_usec);
  fprintf(out, "  prt2_usec: %g\n", info->prt2_usec);
  fprintf(out, "=================================================================\n");

}

//////////////////////////////////////////////////////
// print calibration

void iwrf_calibration_print(FILE *out,
			    const iwrf_calibration_t *calib)
  
{

  fprintf(out, "==================== iwrf_calibration ===========================\n");
  iwrf_packet_info_print(out, &calib->packet);
  
  fprintf(out, "  wavelength_cm: %g\n", calib->wavelength_cm);
  fprintf(out, "  beamwidth_deg_h: %g\n", calib->beamwidth_deg_h);
  fprintf(out, "  beamwidth_deg_v: %g\n", calib->beamwidth_deg_v);
  fprintf(out, "  gain_ant_db_h: %g\n", calib->gain_ant_db_h);
  fprintf(out, "  gain_ant_db_v: %g\n", calib->gain_ant_db_v);
  fprintf(out, "  pulse_width_us: %g\n", calib->pulse_width_us);
  fprintf(out, "  xmit_power_dbm_h: %g\n", calib->xmit_power_dbm_h);
  fprintf(out, "  xmit_power_dbm_v: %g\n", calib->xmit_power_dbm_v);
  fprintf(out, "  two_way_waveguide_loss_db_h: %g\n", calib->two_way_waveguide_loss_db_h);
  fprintf(out, "  two_way_waveguide_loss_db_v: %g\n", calib->two_way_waveguide_loss_db_v);
  fprintf(out, "  two_way_radome_loss_db_h: %g\n", calib->two_way_radome_loss_db_h);
  fprintf(out, "  two_way_radome_loss_db_v: %g\n", calib->two_way_radome_loss_db_v);
  fprintf(out, "  receiver_mismatch_loss_db: %g\n", calib->receiver_mismatch_loss_db);
  fprintf(out, "  radar_constant_h: %g\n", calib->radar_constant_h);
  fprintf(out, "  radar_constant_v: %g\n", calib->radar_constant_v);
  fprintf(out, "  noise_dbm_hc: %g\n", calib->noise_dbm_hc);
  fprintf(out, "  noise_dbm_hx: %g\n", calib->noise_dbm_hx);
  fprintf(out, "  noise_dbm_vc: %g\n", calib->noise_dbm_vc);
  fprintf(out, "  noise_dbm_vx: %g\n", calib->noise_dbm_vx);
  fprintf(out, "  receiver_gain_db_hc: %g\n", calib->receiver_gain_db_hc);
  fprintf(out, "  receiver_gain_db_hx: %g\n", calib->receiver_gain_db_hx);
  fprintf(out, "  receiver_gain_db_vc: %g\n", calib->receiver_gain_db_vc);
  fprintf(out, "  receiver_gain_db_vx: %g\n", calib->receiver_gain_db_vx);
  fprintf(out, "  base_dbz_1km_hc: %g\n", calib->base_dbz_1km_hc);
  fprintf(out, "  base_dbz_1km_hx: %g\n", calib->base_dbz_1km_hx);
  fprintf(out, "  base_dbz_1km_vc: %g\n", calib->base_dbz_1km_vc);
  fprintf(out, "  base_dbz_1km_vx: %g\n", calib->base_dbz_1km_vx);
  fprintf(out, "  sun_power_dbm_hc: %g\n", calib->sun_power_dbm_hc);
  fprintf(out, "  sun_power_dbm_hx: %g\n", calib->sun_power_dbm_hx);
  fprintf(out, "  sun_power_dbm_vc: %g\n", calib->sun_power_dbm_vc);
  fprintf(out, "  sun_power_dbm_vx: %g\n", calib->sun_power_dbm_vx);
  fprintf(out, "  noise_source_power_dbm_h: %g\n", calib->noise_source_power_dbm_h);
  fprintf(out, "  noise_source_power_dbm_v: %g\n", calib->noise_source_power_dbm_v);
  fprintf(out, "  power_meas_loss_db_h: %g\n", calib->power_meas_loss_db_h);
  fprintf(out, "  power_meas_loss_db_v: %g\n", calib->power_meas_loss_db_v);
  fprintf(out, "  coupler_forward_loss_db_h: %g\n", calib->coupler_forward_loss_db_h);
  fprintf(out, "  coupler_forward_loss_db_v: %g\n", calib->coupler_forward_loss_db_v);
  fprintf(out, "  test_power_dbm_h: %g\n", calib->test_power_dbm_h);
  fprintf(out, "  test_power_dbm_v: %g\n", calib->test_power_dbm_v);
  fprintf(out, "  zdr_correction_db: %g\n", calib->zdr_correction_db);
  fprintf(out, "  ldr_correction_db_h: %g\n", calib->ldr_correction_db_h);
  fprintf(out, "  ldr_correction_db_v: %g\n", calib->ldr_correction_db_v);
  fprintf(out, "  phidp_rot_deg: %g\n", calib->phidp_rot_deg);
  fprintf(out, "=================================================================\n");
  
}

//////////////////////////////////////////////////////
// print event_notice

void iwrf_event_notice_print(FILE *out,
			     const iwrf_event_notice_t *note)

{

  fprintf(out, "==================== iwrf_event_notice ==========================\n");
  iwrf_packet_info_print(out, &note->packet);
  
  fprintf(out, "  start_of_sweep: %d\n", note->start_of_sweep);
  fprintf(out, "  end_of_sweep: %d\n", note->end_of_sweep);
  fprintf(out, "  start_of_volume: %d\n", note->start_of_volume);
  fprintf(out, "  end_of_volume: %d\n", note->end_of_volume);
  fprintf(out, "  scan_mode: %s\n", iwrf_scan_mode_to_str(note->scan_mode));
  fprintf(out, "  follow_mode: %s\n", iwrf_follow_mode_to_str(note->follow_mode));
  fprintf(out, "  volume_num: %d\n", note->volume_num);
  fprintf(out, "  sweep_num: %d\n", note->sweep_num);
  fprintf(out, "  cause: %s\n", iwrf_event_cause_to_str(note->cause));
  fprintf(out, "=================================================================\n");

}

//////////////////////////////////////////////////////
// print iwrf_phasecode

void iwrf_phasecode_print(FILE *out,
			  const iwrf_phasecode_t *code)
  
{

  fprintf(out, "==================== iwrf_phasecode =============================\n");
  iwrf_packet_info_print(out, &code->packet);
  
  int seq_length = code->seq_length;
  if (seq_length > IWRF_MAX_PHASE_SEQ_LEN) {
    seq_length = IWRF_MAX_PHASE_SEQ_LEN;
  }
  fprintf(out, "  seq_length: %d\n", seq_length);
  
  int ii = 0;
  for (; ii < seq_length; ii++) {
    fprintf(out, "  Sequence[%d]: phase_deg_h, phase_deg_v: %g, %g",
	    ii,
	    code->phase[ii].phase_deg_h,
	    code->phase[ii].phase_deg_v);
  }
  fprintf(out, "=================================================================\n");

}

//////////////////////////////////////////////////////
// print pulse_header

void iwrf_pulse_header_print(FILE *out,
			     const iwrf_pulse_header_t *pulse)

{
  
  fprintf(out, "==================== iwrf_pulse_header ==========================\n");
  iwrf_packet_info_print(out, &pulse->packet);
  
  fprintf(out, "  pulse_seq_num: %lld\n", (long long) pulse->pulse_seq_num);
  fprintf(out, "  scan_mode: %s\n", iwrf_scan_mode_to_str(pulse->scan_mode));
  fprintf(out, "  follow_mode: %s\n", iwrf_follow_mode_to_str(pulse->follow_mode));
  fprintf(out, "  sweep_num: %d\n", pulse->sweep_num);
  fprintf(out, "  volume_num: %d\n", pulse->volume_num);
  fprintf(out, "  fixed_el: %g\n", pulse->fixed_el);
  fprintf(out, "  fixed_az: %g\n", pulse->fixed_az);
  fprintf(out, "  elevation: %g\n", pulse->elevation);
  fprintf(out, "  azimuth: %g\n", pulse->azimuth);
  fprintf(out, "  prt: %g\n", pulse->prt);
  fprintf(out, "  prt_next: %g\n", pulse->prt_next);
  fprintf(out, "  pulse_width_us: %g\n", pulse->pulse_width_us);
  fprintf(out, "  n_gates: %d\n", pulse->n_gates);
  fprintf(out, "  n_channels: %d\n", pulse->n_channels);
  fprintf(out, "  iq_encoding: %s\n",
          iwrf_iq_encoding_to_str(pulse->iq_encoding));  
  fprintf(out, "  hv_flag: %d\n", pulse->hv_flag);
  fprintf(out, "  antenna_transition: %d\n", pulse->antenna_transition);
  fprintf(out, "  phase_cohered: %d\n", pulse->phase_cohered);
  fprintf(out, "  status: %d\n", pulse->status);
  fprintf(out, "  n_data: %d\n", pulse->n_data);

  fprintf(out, "  iq_offset:");
  int ii = 0;
  for (; ii < IWRF_MAX_CHAN; ii++) {
    fprintf(out, " %d", pulse->iq_offset[ii]);
  }
  fprintf(out, "\n");

  fprintf(out, "  burst_mag:");
  for (ii = 0; ii < IWRF_MAX_CHAN; ii++) {
    fprintf(out, " %g", pulse->burst_mag[ii]);
  }
  fprintf(out, "\n");

  fprintf(out, "  burst_arg:");
  for (ii = 0; ii < IWRF_MAX_CHAN; ii++) {
    fprintf(out, " %g", pulse->burst_arg[ii]);
  }
  fprintf(out, "\n");

  fprintf(out, "  burst_arg_diff:");
  for (ii = 0; ii < IWRF_MAX_CHAN; ii++) {
    fprintf(out, " %g", pulse->burst_arg_diff[ii]);
  }
  fprintf(out, "\n");

  fprintf(out, "  scale: %g\n", pulse->scale);
  fprintf(out, "  offset: %g\n", pulse->offset);
  fprintf(out, "  n_gates_burst: %d\n", pulse->n_gates_burst);

  fprintf(out, "  start_range_m: %g\n", pulse->start_range_m);
  fprintf(out, "  gate_spacing_m: %g\n", pulse->gate_spacing_m);

  if (pulse->event_flags & IWRF_END_OF_SWEEP) {
    fprintf(out, "  event: end_of_sweep\n");
  }
  if (pulse->event_flags & IWRF_START_OF_SWEEP) {
    fprintf(out, "  event: start_of_sweep\n");
  }
  if (pulse->event_flags & IWRF_END_OF_VOLUME) {
    fprintf(out, "  event: end_of_volume\n");
  }
  if (pulse->event_flags & IWRF_START_OF_VOLUME) {
    fprintf(out, "  event: start_of_volume\n");
  }

  fprintf(out, "=================================================================\n");

}

//////////////////////////////////////////////////////
// print rvp8_pulse_header

void iwrf_rvp8_pulse_header_print(FILE *out,
				  const iwrf_rvp8_pulse_header_t *pulse)

{
  
  fprintf(out, "==================== iwrf_rvp8_pulse_header =====================\n");
  iwrf_packet_info_print(out, &pulse->packet);
  
  fprintf(out, "  i_flags: %d\n", pulse->i_flags);
  fprintf(out, "  i_aq_mode: %d\n", pulse->i_aq_mode);
  fprintf(out, "  i_polar_bits: %d\n", pulse->i_polar_bits);
  fprintf(out, "  i_viq_per_bin: %d\n", pulse->i_viq_per_bin);
  fprintf(out, "  i_tg_bank: %d\n", pulse->i_tg_bank);
  fprintf(out, "  i_tx_phase: %d\n", pulse->i_tx_phase);
  fprintf(out, "  i_az: %d\n", pulse->i_az);
  fprintf(out, "  i_el: %d\n", pulse->i_el);
  fprintf(out, "  i_num_vecs: %d\n", pulse->i_num_vecs);
  fprintf(out, "  i_max_vecs: %d\n", pulse->i_max_vecs);
  fprintf(out, "  i_tg_wave: %d\n", pulse->i_tg_wave);
  fprintf(out, "  i_btime_api: %d\n", pulse->i_btime_api);
  fprintf(out, "  i_sys_time: %d\n", pulse->i_sys_time);
  fprintf(out, "  i_prev_prt: %d\n", pulse->i_prev_prt);
  fprintf(out, "  i_next_prt: %d\n", pulse->i_next_prt);
  fprintf(out, "  i_seq_num: %d\n", pulse->i_seq_num);

  fprintf(out, "  uiq_perm:");
  int ii = 0;
  for (ii = 0; ii < 2; ii++) {
    fprintf(out, " %d", pulse->uiq_perm[ii]);
  }
  fprintf(out, "\n");

  fprintf(out, "  i_data_off:");
  for (ii = 0; ii < IWRF_MAX_CHAN; ii++) {
    fprintf(out, " %d", pulse->i_data_off[ii]);
  }
  fprintf(out, "\n");

  fprintf(out, "  f_burst_mag:");
  for (ii = 0; ii < IWRF_MAX_CHAN; ii++) {
    fprintf(out, " %g", pulse->f_burst_mag[ii]);
  }
  fprintf(out, "\n");

  fprintf(out, "  i_burst_arg:");
  for (ii = 0; ii < IWRF_MAX_CHAN; ii++) {
    fprintf(out, " %d", pulse->i_burst_arg[ii]);
  }
  fprintf(out, "\n");

  fprintf(out, "  i_wrap_iq:");
  for (ii = 0; ii < IWRF_MAX_CHAN; ii++) {
    fprintf(out, " %d", pulse->i_wrap_iq[ii]);
  }
  fprintf(out, "\n");
  fprintf(out, "=================================================================\n");
}

//////////////////////////////////////////////////////
// print rvp8_ops_info

void iwrf_rvp8_ops_info_print(FILE *out,
			      const iwrf_rvp8_ops_info_t *info)

{
  char *cp;
  
  fprintf(out, "==================== iwrf_rvp8_ops_info =========================\n");
  iwrf_packet_info_print(out, &info->packet);
  
  fprintf(out, "  i_version: %d\n", info->i_version);
  fprintf(out, "  i_major_mode: %d\n", info->i_major_mode);
  fprintf(out, "  i_polarization: %d\n", info->i_polarization);
  fprintf(out, "  i_phase_mode_seq: %d\n", info->i_phase_mode_seq);
  fprintf(out, "  i_task_sweep: %d\n", info->i_task_sweep);
  fprintf(out, "  i_task_aux_num: %d\n", info->i_task_aux_num);
  fprintf(out, "  i_task_scan_type: %d\n", info->i_task_scan_type);

  
  cp= iwrf_safe_str(info->s_task_name, 32);
  fprintf(out, "  task_name: %s\n", cp);
  free(cp);
  cp= iwrf_safe_str(info->s_site_name, 32);
  fprintf(out, "  site_name: %s\n", cp);
  free(cp);

  fprintf(out, "  i_aq_mode: %d\n", info->i_aq_mode);
  fprintf(out, "  i_unfold_mode: %d\n", info->i_unfold_mode);
  fprintf(out, "  i_pwidth_code: %d\n", info->i_pwidth_code);
  fprintf(out, "  f_pwidth_usec: %g\n", info->f_pwidth_usec);
  fprintf(out, "  f_dbz_calib: %g\n", info->f_dbz_calib);
  fprintf(out, "  i_sample_size: %d\n", info->i_sample_size);
  fprintf(out, "  i_mean_angle_sync: %d\n", info->i_mean_angle_sync);
  fprintf(out, "  i_flags: %d\n", info->i_flags);
  fprintf(out, "  i_playback_version: %d\n", info->i_playback_version);
  fprintf(out, "  f_sy_clk_mhz: %g\n", info->f_sy_clk_mhz);
  fprintf(out, "  f_wavelength_cm: %g\n", info->f_wavelength_cm);
  fprintf(out, "  f_saturation_dbm: %g\n", info->f_saturation_dbm);
  fprintf(out, "  f_range_mask_res: %g\n", info->f_range_mask_res);

  fprintf(out, "  i_range_mask:");
  int ii = 0;
  for (ii = 0; ii < IWRF_RVP8_GATE_MASK_LEN; ii++) {
    fprintf(out, " %d", info->i_range_mask[ii]);
  }
  fprintf(out, "\n");

  fprintf(out, "  f_noise_dbm:");
  for (ii = 0; ii < IWRF_MAX_CHAN; ii++) {
    fprintf(out, " %g", info->f_noise_dbm[ii]);
  }
  fprintf(out, "\n");

  fprintf(out, "  f_noise_stdv_db:");
  for (ii = 0; ii < IWRF_MAX_CHAN; ii++) {
    fprintf(out, " %g", info->f_noise_stdv_db[ii]);
  }
  fprintf(out, "\n");

  fprintf(out, "  f_noise_range_km: %g\n", info->f_noise_range_km);
  fprintf(out, "  f_noise_prf_hz: %g\n", info->f_noise_prf_hz);

  fprintf(out, "  i_gparm_latch_sts:");
  for (ii = 0; ii < 2; ii++) {
    fprintf(out, " %d", info->i_gparm_latch_sts[ii]);
  }
  fprintf(out, "\n");

  fprintf(out, "  i_gparm_immed_sts:");
  for (ii = 0; ii < 6; ii++) {
    fprintf(out, " %d", info->i_gparm_immed_sts[ii]);
  }
  fprintf(out, "\n");

  fprintf(out, "  i_gparm_diag_bits:");
  for (ii = 0; ii < 4; ii++) {
    fprintf(out, " %d", info->i_gparm_diag_bits[ii]);
  }
  fprintf(out, "\n");

  cp= iwrf_safe_str(info->s_version_string, 12);
  fprintf(out, "  version_string: %s\n", cp);
  free(cp);

  fprintf(out, "=================================================================\n");

}
void iwrf_ui_tasklist_print(FILE *out, iwrf_ui_tasklist_full_t *tlp)
{	int i;
	fprintf(out, "Tasklist Contents:\n");
	for(i=0; i<tlp->num_list_items; ++i) fprintf(out,"%3d %s\n", i, tlp->tasks[i].name);
	fprintf(out, "Tasklist Schedule:\n");
	iwrf_ui_schedule_print(out, &tlp->tasklist_schedule);
	fprintf(out, "\n");
}
//////////////////////////////////////////////////////
// Return a string formed safely from a char* array
// Null-termination of the input string is guaranteed.
// calling program should free returned string

char *iwrf_safe_str(const char *str, int maxLen)

{

  char *safechar = malloc(maxLen + 1);
  strncpy(safechar, str, maxLen);
  safechar[maxLen] = '\0';
  return safechar;

}
