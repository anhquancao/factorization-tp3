rank=10
etas=(0.00001 0.001 0.1 1)
iter=20

for eta in "${etas[@]}"; do    
    ./sgd.sh $rank $eta $iter;
done
