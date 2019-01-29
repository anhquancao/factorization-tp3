ranks=(10 20 30 40 50 60)
eta=0.00001
iter=20

for rank in "${ranks[@]}"; do
    ./sgd.sh $rank $eta $iter;
done
