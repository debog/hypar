FILE="physics.inp"

if [ "$p_adv_file" == "none" ]; then
/bin/cat <<EOM > $FILE
begin
  advection ${p_adv_speed}
end
EOM
else
/bin/cat <<EOM > $FILE
begin
  advection_filename ${p_adv_file}
end
EOM
fi
