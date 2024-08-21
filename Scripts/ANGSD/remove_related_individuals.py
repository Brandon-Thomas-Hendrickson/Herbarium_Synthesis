import pandas as pd

# Load the ngsrelate output
ngsrelate_output = "newres" 

# Read the ngsrelate output into a DataFrame
df = pd.read_csv(ngsrelate_output, delim_whitespace=True)

relatedness_threshold = 0.25

# Identify pairs with high relatedness
highly_related_pairs = df[df['relatedness'] > relatedness_threshold]

# Create a set to keep track of individuals to remove
individuals_to_remove = set()

# Iterate through the highly related pairs
for index, row in highly_related_pairs.iterrows():
    ind1 = row['ind1']
    ind2 = row['ind2']
    
    # If neither individual is already marked for removal, mark one for removal
    if ind1 not in individuals_to_remove and ind2 not in individuals_to_remove:
        individuals_to_remove.add(ind2)

# Output the list of individuals to remove
with open("individuals_to_remove.txt", "w") as f:
    for individual in individuals_to_remove:
        f.write(f"{individual}\n")

print("List of individuals to remove has been written to individuals_to_remove.txt")