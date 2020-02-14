FILE="physics.inp"

/bin/cat <<EOM > $FILE
begin
  gamma     ${p_gamma}
  upwinding ${p_upwinding}
  Pr        ${p_pr}
  Minf      ${p_minf}
  gravity   ${p_gravity}
  rho_ref   ${p_rho_ref}
  p_ref     ${p_p_ref}
  HB        ${p_HB}
  R         ${p_R}
end
EOM
