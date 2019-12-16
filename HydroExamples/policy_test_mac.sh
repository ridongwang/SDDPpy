
python ./Hydro_ARx_SP.py --T=12 --R=10 --max_iter=1000 --max_time=500 --sim_iter=1000 --lines_freq=10 --dro_r=0 --lag=1 --N=10 --dynamic_sampling=False --multicut=True
for c in {-2,-1,0,1,2}
do
		for b in {1,5}
        do
			r=$(bc -l <<< "$b*(10^$c)")
			## python ./Hydro_ARx_DW_Phlipott.py --T=12 --R=10 --max_iter=100001 --max_time=100 --sim_iter=1000 --lines_freq=10 --dro_r=$r --lag=1 --N=10 --dynamic_sampling=False --multicut=True
			python ./Hydro_ARx_DW_DUAL.py --T=12 --R=10 --max_iter=1000 --max_time=500 --sim_iter=1000 --lines_freq=20 --dro_r=$r --lag=1 --N=10 --dynamic_sampling=False --multicut=True

        done
done
