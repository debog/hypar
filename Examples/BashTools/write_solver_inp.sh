FILE="solver.inp"

if [ "x$p_size_exact" == "x" ]; then
  export p_size_exact=$p_size
fi

/bin/cat <<EOM >$FILE
begin
  ndims              $p_ndims
  nvars              $p_nvars
  size               $p_size
  size_exact         $p_size_exact
  iproc              $p_iproc
  ghost              $p_ghost
  n_iter             $p_niter
  restart_iter       0
  time_scheme        $p_ts
  time_scheme_type   $p_tstype
  hyp_space_scheme   $p_hyp_scheme
  hyp_flux_split     $p_hyp_flux_split
  hyp_interp_type    $p_hyp_int_type
  par_space_type     $p_par_type
  par_space_scheme   $p_par_scheme
  dt                 $p_dt
  conservation_check $p_cons_check
  screen_op_iter     $p_screen_op_iter
  file_op_iter       $p_file_op_iter
  op_file_format     $p_op_format
  ip_file_type       $p_ip_type
  input_mode         $p_input_mode $p_n_io
  output_mode        $p_output_mode $p_n_io
  op_overwrite       $p_op_overwrite
  model              $p_model
end
EOM
