import matplotlib.ticker as ticker
import matplotlib.dates as mdates
import matplotlib.pyplot as plt
import pandas as pd
import subprocess
import datetime
import sys
import os


def parse_fasta(fasta_file):
    """Parse the FASTA file to extract dates."""
    dates = []
    with open(fasta_file, 'r') as file:
        for line in file:
            if line.startswith('>'):
                # Extract the date from the header line
                parts = line.strip().split('|')
                date = parts[-1].strip()
                dates.append(date)
    return dates

def combine_data_with_dates(nc_files, fasta_file):
    dates = parse_fasta(fasta_file)
    combined_data = pd.DataFrame({'Date': pd.to_datetime(dates)})

    for file in nc_files:
        data = pd.read_csv(file, sep="\t", header=None, usecols=[1])
        column_name = os.path.splitext(os.path.basename(file))[0]  # Extract file name without extension
        combined_data[column_name] = data.iloc[:, 0]  # Add data as a new column
    
    combined_file = "combined_data.csv"
    combined_data.to_csv(combined_file, index=False)

    return combined_file


def run_altair(input_file, window_size, output_file):
    # Read the file using pandas
    data = pd.read_csv(input_file, sep='\t', header=None)
    temp_file = f"temp_{window_size}.txt"
    # Write the first two columns to a temporary file
    data.iloc[:, :2].to_csv(temp_file, sep='\t', index=False, header=False)
    # Run AltaiR
    cmd = f"./AltaiR average -v -c 2 -p -w {window_size} < {temp_file} > {output_file}"
    subprocess.run(cmd, shell=True, check=True)
    # Clean up the temporary file
    os.remove(temp_file)

def get_min_max_time(input_file):
    data = pd.read_csv(input_file, sep='\t', header=None)
    min_time = data.iloc[:, 0].min()
    max_time = data.iloc[:, 0].max()
    return min_time, max_time

def plot_data(combined_file, yrange, colors, ylabel, output_file):
    data = pd.read_csv(combined_file)
    data['Date'] = pd.to_datetime(data['Date'])  # Ensure dates are parsed correctly

    plt.figure(figsize=(16, 4))

    # Adjust this list to match the column names in your CSV
    columns_to_plot = ['.NCFiltered_5', '.NCFiltered_20', '.NCFiltered_100']

    for idx, column in enumerate(columns_to_plot):
        # Extract window size from the column name for the legend
        window_size = column.split('_')[-1]
        label = f'w={window_size}'
        plt.plot(data['Date'], data[column], color=colors[idx], label=label)

    # Set date format on x-axis
    locator = mdates.AutoDateLocator(minticks=3, maxticks=7)
    formatter = mdates.ConciseDateFormatter(locator)
    plt.gca().xaxis.set_major_locator(locator)
    plt.gca().xaxis.set_major_formatter(formatter)
    
    # Set larger font size and bold font weight for x-axis label
    plt.xlabel('Date', fontsize='xx-large', fontweight='bold')
    
    # Rotate x-axis ticks
    plt.xticks(rotation=45)

    # Formatting for y-axis
    plt.gca().yaxis.set_major_formatter(ticker.FormatStrFormatter('%.4f'))

    # Increase legend font size
    plt.legend(loc='upper right', fontsize='xx-large')

    # Set larger font size and bold font weight for y-axis label
    plt.ylabel(ylabel, fontsize='xx-large', fontweight='bold')

    # Set tick parameters
    plt.tick_params(axis='both', which='major', labelsize='xx-large')  # <-- Update here

    # Set x-axis range from January 2020 to January 2022
    plt.xlim(datetime.datetime(2020, 1, 1), datetime.datetime(2022, 1, 31))
    plt.ylim(yrange)
    plt.grid(True)
    plt.tight_layout()
    plt.savefig(output_file)
    plt.show()  # Optionally display the plot



    
def main():
    if len(sys.argv) != 5:
        print("Usage: python script.py input_file fasta_file y_min y_max")
        sys.exit(1)
    
    input_file = sys.argv[1]
    fasta_file = sys.argv[2]
    y_min = float(sys.argv[3])
    y_max = float(sys.argv[4])
    output_file = f"NCProfile{input_file}.pdf"

    window_sizes = [5, 20, 100]
    temp_files = [f".NCFiltered_{w}.txt" for w in window_sizes]
    colors = ['#322152', '#009900', '#CC0000']  # blue, green, red
    # Update labels to match column names
    labels = [os.path.splitext(os.path.basename(file))[0] for file in temp_files]

    for w, temp_file in zip(window_sizes, temp_files):
        run_altair(input_file, w, temp_file)

    combined_file = combine_data_with_dates(temp_files, fasta_file)
    plot_data(combined_file, (y_min, y_max), colors, "NC", output_file)



if __name__ == "__main__":
    main()
