FILE=( HMEC_unique BT549_unique HMEC_common BT549_common )
for f in ${FILE[@]}; do
  echo ${f%%_*}
done