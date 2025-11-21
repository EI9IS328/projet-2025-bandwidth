import csv

input_file = "data.txt"
output_file = "data.csv"

with open(input_file, "r") as txt_file, open(output_file, "w", newline="") as csv_file:
    writer = csv.writer(csv_file)
    
    for line in txt_file:
        # Skip empty lines
        if not line.strip():
            continue
        
        # Split by whitespace into columns
        row = line.split()
        writer.writerow(row)

print("Conversion complete! CSV saved as:", output_file)
