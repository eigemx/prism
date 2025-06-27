import pandas as pd
import json

if __name__ == "__main__":
    df = pd.read_csv("foam_solution_fine.csv")

    p_list = df['p'].tolist()
    u_components = df[['U_0', 'U_1']].values.tolist()

    u_vectors = [u + [0.0] for u in u_components]

    output = {
        "U": u_vectors,
        "p": p_list
    }

    print(len(p_list))
    print(len(u_vectors))

    json_string = json.dumps(output, indent=2)

    with open("foam_fields.json", "w") as f:
        f.write(json_string)
