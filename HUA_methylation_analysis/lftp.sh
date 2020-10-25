#!/bin/bash

HOST=
USER=
PASS=

cd rrbs_clean_data_lmer
lftp -u ${USER},${PASS} sftp://${HOST} <<EOF
cd /home/sftp-epic-omics/HUA_EWAS/rrbs_clean_data_lmer;
mirror --parallel=15 --verbose
bye;
EOF
