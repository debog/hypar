FILE="physics.inp"

if [ "$p_adv_file" == "none" ]; then
/bin/cat <<EOM > $FILE
begin
  advection       ${p_adv_speed}
  centered_flux   ${p_linadv_centered_flux}
end
EOM
else
/bin/cat <<EOM > $FILE
begin
  advection_filename  ${p_adv_file}
  centered_flux       ${p_linadv_centered_flux}
end
EOM
fi
