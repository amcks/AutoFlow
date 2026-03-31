#!/bin/bash

input_file="cleaned_smiles.txt"
counter=0

# Path to your script
autoflow_script="autoflow_v0_9.sh"

# Ensure the autoflow_v0_8.sh script exists in the current directory
if [[ ! -f "$autoflow_script" ]]; then
    echo "Error: $autoflow_script not found in the current directory."
    exit 1
fi

while IFS= read -r smiles || [[ -n "$smiles" ]]; do
    # Skip empty lines
    if [[ -z "$smiles" ]]; then
        ((counter++))
        continue
    fi

    # Create a numbered directory
    dir="$counter"
    mkdir -p "$dir"

    # If the SMILES string is the dummy atom, skip further processing
    if [[ "$smiles" == "[*]" ]]; then
        echo "Skipping dummy atom at line $counter"
        ((counter++))
        continue
    fi

    # Copy the autoflow_v0_8.sh script into the directory
    cp "$autoflow_script" "$dir/"
    
    cd "$dir"
    # Create the SLURM submit.sh script inside the directory
    cat > "submit.sh" <<EOF
#!/bin/bash

#SBATCH -p hawkcpu,haswell
#SBATCH -t 10:00:00
#SBATCH -N 1
#SBATCH --ntasks-per-node=1
#SBATCH -J vasp
#SBATCH --qos=nogpu
#SBATCH --job-name="EO_RING_$counter"

source /etc/profile.d/zlmod.sh

ulimit -s unlimited

./autoflow_v0_9.sh -s Ag -m 1,1,1 -a "$smiles" -l 4.13 -p fcc

exit
EOF

    # Submit the SLURM script
    sbatch "submit.sh"
    cd ..

    ((counter++))
done < "$input_file"
