FILE="physics.inp"

/bin/cat <<EOM > $FILE
begin
  gamma     ${p_gamma}
  gravity   ${p_gravity}
  upwinding ${p_upwinding}
end
EOM
