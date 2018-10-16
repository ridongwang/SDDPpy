for c in {-4,-3,-2}
do
		for b in {1,1.5,2,3,4,5,6,7,8,9}
        do
			r=$(bc -l <<< "$b*(10^$c)")
			python ./Hydro_ARx_DW_DUAL.py --T=12 --R=2 --max_iter=10000 --max_time=60 --sim_iter=1000 --lines_freq=10 --dro_r=$r --lag=1 --N=5 --dynamic_sampling=True --multicut=True --DW_extended=3 --DW_sampling=cvx_hull2
        done
done
