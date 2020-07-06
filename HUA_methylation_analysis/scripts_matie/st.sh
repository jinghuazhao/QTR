ls *.R | parallel -C' ' 'dos2unix {}'
ls *.R | parallel -C' ' 'sed -i "s/rrbs_clean_data/rrbs_clean_data_matie/g" {}'
