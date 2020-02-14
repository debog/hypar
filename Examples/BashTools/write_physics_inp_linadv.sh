FILE="physics.inp"

/bin/cat <<EOM > $FILE
begin
  advection ${p_adv_speed}
end
EOM
