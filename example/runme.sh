#!/usr/bin/bash
deactivate
eval "$(conda shell.bash hook)"
conda activate cgcms-poem
cd ..
bash bin/run_poem_cgcms.sh -f /example -a n -p pro >>poem.log
conda deactivate