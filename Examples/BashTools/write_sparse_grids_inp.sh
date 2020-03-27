FILE="sparse_grids.inp"

/bin/cat <<EOM > $FILE
begin
  log2_imin           $p_imin
  interp_order        $p_interp_order
  write_sg_solutions  $p_write_sg_solutions
  write_sg_errors     $p_write_sg_errors
end
EOM
