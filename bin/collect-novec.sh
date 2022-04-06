rm intel-baseline.txt
touch intel-baseline.txt
names=("s123" "s124" "s161" "s1161" "s253" "s258" "s271" "s273" "s274" "s277" "s278" "s2711" "s2712" "s314" "s315" "s316" "s3111" "s3113" "s341" "s342" "s343" "s443" "vif")
for name in "${names[@]}"
do
	  for i in `seq 1 10`
		    do
			      $1/tsvc_novec_default $name
        done
done
