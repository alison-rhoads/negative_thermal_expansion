pheasy -s --dim 3 3 3 -w 3 --c3 5 --nbody 2 3 --ndata 50
echo "step 1 done"
pheasy -s --read_scell --dim 3 3 3 -w 3 --c3 5 --nbody 2 3 --ndata 50
echo "step 2 done"
pheasy -c --read_scell --dim 3 3 3 -w 3 --ndata 50
echo "step 3 done"
pheasy -d --read_scell --dim 3 3 3 --ndata 50 --disp_file
echo "step 4 done"
pheasy -f --read_scell --dim 3 3 3 -w 3 --model LASSO --fix_fc2 --ndata 50 --hdf5 --full_ifc > fit.out
echo "step 5 done"
echo "FITTING COMPLETE!"
