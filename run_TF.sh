export OMP_NUM_THREADS=4
python CSHMM_TF_train_release.py \
--data_file treutlein2014 \
-tf tfDNA_predicted_100.txt.update \
--structure_file init_cluster_treutlein2014.txt \
--n_split 100 -ng 16000 \
--n_iteration 2 \
--cross_validation 0 \
--random_seed 5 \
-k 10 \
-opt genlasso \
--model_name lung_developmental_TF_v7 \
>>  lung_developmental_TF_v7.log 2>&1 
