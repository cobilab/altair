import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
import numpy as np
from matplotlib.ticker import FormatStrFormatter
import matplotlib.font_manager as font_manager


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

def read_frequency_data(freq_file):
    """Read the frequency data from the CSV file."""
    with open(freq_file, 'r') as file:
        return [line.strip() for line in file]

def merge_data(dates, freq_data):
    """Merge dates with frequency data."""
    combined_data = [f"{freq}\t{date}" for freq, date in zip(freq_data, dates)]
    return combined_data

def write_combined_data(combined_data, output_file):
    """Write the combined data to a new file."""
    with open(output_file, 'w') as file:
        for line in combined_data:
            file.write(line + '\n')

def read_data(file_path):
    return pd.read_csv(file_path, sep='\t', header=None, names=['A', 'C', 'T', 'G', 'Date'])



def plot_data(data, window_size=20):
    bases = ['A', 'G', 'C', 'T']
    colors = ['skyblue', 'salmon', 'lightgreen', 'khaki']
    fig, axs = plt.subplots(2, 2, figsize=(22, 8))  # Wider plot

    # Convert 'Date' to datetime for better x-axis handling
    data['Date'] = pd.to_datetime(data['Date'])
    data = data[(data['Date'] >= '2020-01-01') & (data['Date'] <= '2022-01-01')]

    # Determine the range of data to set plot limits
    date_min, date_max = data['Date'].min(), data['Date'].max()

    locator = mdates.AutoDateLocator(minticks=3, maxticks=7)
    formatter = mdates.ConciseDateFormatter(locator)

    for i, base in enumerate(bases):
        row, col = divmod(i, 2)

        # Calculate moving average
        smoothed_data = data[base].rolling(window=window_size, min_periods=1).mean()

        axs[row, col].plot(data['Date'], smoothed_data, label=f'Base {base}', color=colors[i], linewidth=2)
        
        # Numerical format with 4 decimal places
        axs[row, col].yaxis.set_major_formatter(FormatStrFormatter('%.4f'))
        
        # Increase legend font size
        axs[row, col].legend(loc='upper left', fontsize='xx-large')
        
        # Set larger font size for y-axis label
        axs[row, col].set_ylabel('Frequency', fontsize='xx-large', fontweight='bold')
        
        # Set larger font size for y-axis label
        axs[row, col].set_xlabel('Date', fontsize='xx-large', fontweight='bold')

        # Set plot limits
        axs[row, col].set_xlim([date_min, date_max])

        # Apply the locator and formatter
        axs[row, col].xaxis.set_major_locator(locator)
        axs[row, col].xaxis.set_major_formatter(formatter)

        # Increase tick visibility
        axs[row, col].tick_params(axis='x', rotation=45, labelsize='xx-large')
        axs[row, col].tick_params(axis='y', labelsize='xx-large')


    plt.tight_layout()
    plt.show()
    plt.savefig('base_frequencies_plot.pdf')


def check_base_pair_correlation(data_file):
    data = pd.read_csv(data_file, sep='\t', header=None, names=['A', 'C', 'T', 'G', 'Date'])
    
    at_correlation = data['A'].corr(data['T'])
    cg_correlation = data['C'].corr(data['G'])
    
    print(f"Correlation between A and T frequencies: {at_correlation}")
    print(f"Correlation between C and G frequencies: {cg_correlation}")

def main():
    # File paths
    fasta_file = 'sorted_output.fa'
    freq_file = 'freq-data.csv'
    output_file = 'freq-date-data.csv'

    # Process the files
    dates = parse_fasta(fasta_file)
    freq_data = read_frequency_data(freq_file)
    combined_data = merge_data(dates, freq_data)
    write_combined_data(combined_data, output_file)
    
    print("Combined data written to", output_file)
    file_path = 'freq-date-data.csv'
    data = read_data(file_path)
    plot_data(data)
    

main()



