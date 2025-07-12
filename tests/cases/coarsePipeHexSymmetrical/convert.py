import json
import sys

import pandas as pd

if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Usage: convert.py <foam_solution.csv>")
        sys.exit(1)

    csv_file = sys.argv[1]

    df = pd.read_csv(csv_file)
    p_list = df["p"].tolist()
    u_components = df[["U_0", "U_1"]].values.tolist()

    u_vectors = [u + [0.0] for u in u_components]

    output = {"U": u_vectors, "p": p_list}

    json_string = json.dumps(output, indent=2)

    with open("foam_fields.json", "w") as f:
        f.write(json_string)
