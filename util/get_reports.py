import re
import csv
import os

def extract_and_convert(line, data_type):
    # Extract the number after the colon or equal sign, possibly with a unit or percentage
    match = re.search(r'[:=]\s*([-+]?\d*\.?\d+)', line)
    if match:
        number = match.group(1)
        if '%' in line:
            # Convert percentage to decimal and round to four decimal places
            return round(float(number) / 100, 4)
        try:
            if data_type == 'float':
                return round(float(number), 4)
            elif data_type == 'int':
                return int(number)
            else:
                return number  # Default to string
        except ValueError:
            return number  # Default to string if conversion fails
    return None

def process_files(directory):
    # Keywords and their corresponding data types
    keywords = {
        'Ps3_factor': 'int',
        'Ps4_factor': 'int',
        'BCM1 Beam Cut Current': 'float',
        'BCM2 Beam Cut Current': 'float',
        'BCM4A Beam Cut Current': 'float',
        'BCM4C Beam Cut Current': 'float',
        'BCM1  Beam Cut Charge': 'float',
        'BCM2  Beam Cut Charge': 'float',
        'BCM4A Beam Cut Charge': 'float',
        'BCM4C Beam Cut Charge': 'float',
        'Pre-Scaled Ps3 HMS Computer Live Time': 'float',
        'Pre-Scaled Ps4 HMS Computer Live Time': 'float',
        'E SING FID TRACK EFFIC': 'float',
        'hTRIG3': 'int',
        'hTRIG4': 'int',
        'Run Num': 'int'
    }
    header = ['Run Num'] + list(keywords.keys())[:-1]  # Exclude 'Run Num' from parameters for the header
    data = []

    # Iterate through each file in the directory
    for filename in os.listdir(directory):
        if filename.endswith(".report"):  # Check for text files
            file_path = os.path.join(directory, filename)
            params = {key: None for key in keywords}  # Initialize dictionary to store values

            with open(file_path, 'r') as file:
                for line in file:
                    for key in keywords:
                        if key in line:
                            params[key] = extract_and_convert(line, keywords[key])

            # Collect the data from the dictionary into a list maintaining the order from 'header'
            if params['Run Num'] is not None:
                data.append([params[key] for key in header])

    # Write to CSV
    with open('output.csv', 'w', newline='') as csvfile:
        writer = csv.writer(csvfile)
        writer.writerow(header)
        writer.writerows(data)

# Usage example
process_files('data/reports2')
