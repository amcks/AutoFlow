from collections import Counter
import csv
import json

def parse_reaction(step):
    reactants, products = step.split(">>")
    return {
        "reactants": Counter(reactants.split(".")),
        "products": Counter(products.split("."))
    }


def parse_mechanisms(filepath):
    with open(filepath, "r") as f:
        steps = []
        rules = {}
        length = None
        in_rules = False

        for line in f:
            line = line.strip()

            if not line:
                continue

            # Reaction step
            if ">>" in line:
                if in_rules:
                    yield {
                        "steps": steps,
                        "length": length,
                        "rules": rules
                    }
                    steps, rules, length = [], {}, None
                    in_rules = False

                steps.append(parse_reaction(line))
                continue

            # Length line
            if line.startswith("the above pathway is of length"):
                length = int(line.split()[6])
                in_rules = True
                continue

            # Rule lines
            if line.startswith("rule") and in_rules:
                parts = line.split()
                rules[parts[1]] = int(parts[-2])
                continue

        if steps:
            yield {
                "steps": steps,
                "length": length,
                "rules": rules
            }


def load_energy_map(csv_path):
    energy_map = {}

    with open(csv_path, newline="") as f:
        reader = csv.DictReader(f)
        for row in reader:
            smi = row["original"]
            energy = row.get("energy")

            if smi and energy:
                energy_map[smi] = float(energy)

    # hardcode
    energy_map.setdefault("[{M}]", -171.1818084716797) #MACE-MP energy

    return energy_map


def compute_step_energy(step, energy_map):
    def total(counter):
        return sum(energy_map[smi] * count for smi, count in counter.items() if smi!="[{M}]")
    
    n_reactants = sum(count for smi, count in step["reactants"].items() if smi!="[{M}]")
    n_products = sum(count for smi, count in step["products"].items() if smi!="[{M}]")
    Surf_corr = (n_products - n_reactants)*energy_map["[{M}]"]
    
    return total(step["products"]) - total(step["reactants"]) - Surf_corr


def evaluate_mechanism(mech, energy_map, current_best=float("inf")):
    max_barrier = float("-inf")

    for step in mech["steps"]:
        try:
            dE = compute_step_energy(step, energy_map)
        except KeyError as e:
            # Skip mechanism for missing species
            return None

        if dE > max_barrier:
            max_barrier = dE

        # Early pruning
        if max_barrier >= current_best:
            return None

    return max_barrier


def find_best_mechanism(mech_file, energy_csv):
    energy_map = load_energy_map(energy_csv)

    best_mech = None
    best_barrier = float("inf")

    for mech in parse_mechanisms(mech_file):
        barrier = evaluate_mechanism(mech, energy_map, best_barrier)

        if barrier is None:
            continue

        if barrier < best_barrier:
            best_barrier = barrier
            best_mech = mech
    
    stepwise_barr = []
    for step in best_mech["steps"]:
        dE = compute_step_energy(step, energy_map)
        stepwise_barr.append(dE)

    return best_mech, best_barrier, stepwise_barr


if __name__ == "__main__":
    best_mech, barrier, stepwise = find_best_mechanism(
        "CO2pathways.txt",
        "output_MACE-MP.csv"
    )

    print("Best barrier:", barrier)
    if best_mech:
        print("Mechanism length:", best_mech["length"])
        print(stepwise)
        print(best_mech["steps"])
