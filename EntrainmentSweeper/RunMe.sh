rm -r Job* 
module load intel/19.4
counter=$(python3 Python.py 2>&1)
for ((i = 1; i<= $counter;i++)); do
    cp -r ./Template/* ./Job"$i"_*
    cd ./Job"$i"_*
    make
    #./Runme.sh
    cd ../
done

