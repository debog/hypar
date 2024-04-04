# Defaut values of all the input parameters for HyPar

# sparse_grids.inp
export p_imin=3
export p_interp_order=6
export p_write_sg_solutions="no"
export p_write_sg_errors="yes"

# solver.inp
export p_ndims=1
export p_nvars=1
export p_size=100
export p_size_exact=
export p_iproc=1
export p_ghost=3
export p_niter=200
export p_ts="rk"
export p_tstype="ssprk3"
export p_hyp_scheme="weno5"
export p_hyp_flux_split="no"
export p_hyp_int_type="components"
export p_par_type="nonconservative-1.5stage"
export p_par_scheme="4"
export p_dt="0.005"
export p_cons_check="yes"
export p_screen_op_iter=1
export p_file_op_iter=99999999
export p_op_format="none"
export p_ip_type="binary"
export p_input_mode="serial"
export p_output_mode="serial"
export p_n_io=
export p_op_overwrite="yes"
export p_model="linear-advection-diffusion-reaction"

# boundary.inp
export p_nb=2
BC_COUNTER=0
while [ $BC_COUNTER -lt $p_nb ]; do
  varname="p_bctype_${BC_COUNTER}"
  eval export "$varname"="periodic"
  varname="p_dim_${BC_COUNTER}"
  eval export "$varname"="0"
  varname="p_limits_${BC_COUNTER}"
  eval 'export "$varname"="0 0"'
  let BC_COUNTER=BC_COUNTER+1
done
BC_COUNTER=0
while [ $BC_COUNTER -lt $p_nb ]; do
  varname="p_face_${BC_COUNTER}"
  eval export "$varname"="1"
  let BC_COUNTER=BC_COUNTER+2
done
BC_COUNTER=1
while [ $BC_COUNTER -lt $p_nb ]; do
  varname="p_face_${BC_COUNTER}"
  eval export "$varname"="-1"
  let BC_COUNTER=BC_COUNTER+2
done

# weno.inp
export p_mapped=0
export p_borges=0
export p_yc=1
export p_nl=0
export p_eps="0.000001"
export p_p=2.0
export p_rc="0.3"
export p_xi="0.001"
export p_tol="1e-16"

# lusolver.inp
export p_lutype="jacobi"
export p_norm=1
export p_maxiter=10
export p_atol="1e-12"
export p_rtol="1e-10"
export p_verbose=0

# linear advection
export p_adv_speed="1.0"
export p_adv_file="none"
export p_linadv_centered_flux="no"

# Euler/Navier-Stokes equations
export p_gamma="1.4"
export p_gravity="0.0"
export p_pr="0.72"
export p_re="-1"
export p_minf="1.0"
export p_rho_ref="1.0"
export p_p_ref="1.0"
export p_HB=0
export p_R="1.0"
export p_upwinding="roe"
