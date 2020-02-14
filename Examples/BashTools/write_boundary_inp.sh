FILE="boundary.inp"

/bin/cat <<EOM >$FILE
$p_nb
EOM

BC_COUNTER=0
while [ $BC_COUNTER -lt $p_nb ]; do
  varname_bctype=p_bctype_${BC_COUNTER}
  varname_dim=p_dim_${BC_COUNTER}
  varname_face=p_face_${BC_COUNTER}
  varname_limits=p_limits_${BC_COUNTER}
/bin/cat <<EOM >>$FILE
${!varname_bctype} ${!varname_dim} ${!varname_face} ${!varname_limits}
EOM
  let BC_COUNTER=BC_COUNTER+1
done
