import matplotlib.pyplot as plt
import sys
import numpy as np
import glob
import math
import matplotlib.dates as mdates
from datetime import datetime


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


def parse_file(file_path):
    """
    Parses the given file and returns two dictionaries:
    1. mRAW_count_per_sequence: {sequence_number: count_of_mRAWs}
    2. average_gc_content_per_sequence: {sequence_number: average_gc_content}
    """
    mraw_count_per_sequence = {}
    gc_content_sum_per_sequence = {}
    mraw_lengths_per_sequence = {}

    with open(file_path, 'r') as file:
        for line in file:
            sequence, _, mraw = line.strip().split('\t')
            sequence = int(sequence)

            # Count mRAWs per sequence
            mraw_count_per_sequence[sequence] = mraw_count_per_sequence.get(sequence, 0) + 1

            # Calculate GC content for the mRAW
            gc_count = sum(nucleotide in ['G', 'C'] for nucleotide in mraw)
            gc_content = gc_count / len(mraw)

            # Sum up GC content and lengths for calculating averages later
            gc_content_sum_per_sequence[sequence] = gc_content_sum_per_sequence.get(sequence, 0) + gc_content
            mraw_lengths_per_sequence[sequence] = mraw_lengths_per_sequence.get(sequence, 0) + 1

    # Calculating average GC content per sequence
    average_gc_content_per_sequence = {seq: gc_content_sum_per_sequence[seq] / mraw_lengths_per_sequence[seq]
                                       for seq in gc_content_sum_per_sequence}

    if len(mraw_count_per_sequence) < 10:
        return None  # Indicate insufficient data
    else:
        return mraw_count_per_sequence, average_gc_content_per_sequence



def moving_average(data, window_size):
    """Simple moving average of the data."""
    cumsum = np.cumsum(data, dtype=float)
    cumsum[window_size:] = cumsum[window_size:] - cumsum[:-window_size]
    return cumsum[window_size - 1:] / window_size


def plot_data(all_data, dates, window_size=100):
    num_k_values = len(all_data)
    num_cols = 2
    num_rows = math.ceil(num_k_values / (num_cols / 2))

    # Determine the split row index for the main and supplementary figures
    split_row_index = num_rows - 2

    # Prepare dates for plotting
    plot_dates = [datetime.strptime(date, '%Y-%m-%d') for date in dates]  # Adjust date format as needed
    sorted_keys = sorted(all_data.keys())

    # Main figure
    fig, main_axes = plt.subplots(nrows=split_row_index, ncols=num_cols, figsize=(22, 5 * split_row_index), squeeze=False)
    # Supplementary figure
    sup_fig, sup_axes = plt.subplots(nrows=2, ncols=num_cols, figsize=(22, 10), squeeze=False)

    for i, k in enumerate(sorted_keys):
        data = all_data[k]
        row = i // (num_cols // 2)
        col = (i % (num_cols // 2)) * 2

        mraw_count_per_sequence, average_gc_content_per_sequence = data
        mraw_counts = [mraw_count_per_sequence[seq] for seq in sorted(mraw_count_per_sequence.keys())]
        average_gc_contents = [average_gc_content_per_sequence[seq] * 100 for seq in sorted(average_gc_content_per_sequence.keys())]

        smoothed_mraw_counts = moving_average(mraw_counts, window_size)
        smoothed_gc_contents = moving_average(average_gc_contents, window_size)

        label = f'k={k}'

        if row < split_row_index:
            axes = main_axes
        else:
            axes = sup_axes
            row -= split_row_index  # Adjust row index for supplementary figure

        ax1 = axes[row, col]  # mRAW counts plot
        ax1.plot(plot_dates[window_size - 1:], smoothed_mraw_counts, color='royalblue', linewidth=2, label=label)
        ax1.set_xlabel('Date', fontsize='xx-large', fontweight='bold')
        ax1.set_ylabel('Count of mRAWs', fontsize='xx-large', fontweight='bold')
        ax1.grid(True, linestyle='--', alpha=0.7)
        ax1.tick_params(axis='both', which='major', labelsize='xx-large')
        ax1.tick_params(axis='both', which='minor', labelsize='xx-large')
        ax1.legend(loc='upper right', fontsize='xx-large')

        ax2 = axes[row, col + 1]  # GC content plot
        ax2.plot(plot_dates[window_size - 1:], smoothed_gc_contents, color='forestgreen', linewidth=2, label=label)
        ax2.set_xlabel('Date', fontsize='xx-large', fontweight='bold')
        ax2.set_ylabel('Average GC Content (%)', fontsize='xx-large', fontweight='bold')
        ax2.grid(True, linestyle='--', alpha=0.7)
        ax2.tick_params(axis='both', which='major', labelsize='xx-large')
        ax2.tick_params(axis='both', which='minor', labelsize='xx-large')
        ax2.legend(loc='upper right', fontsize='xx-large')

        for ax in [ax1, ax2]:
            ax.xaxis.set_major_locator(mdates.YearLocator())
            ax.xaxis.set_major_formatter(mdates.DateFormatter('%Y'))
            ax.xaxis.set_minor_locator(mdates.MonthLocator())
            ax.xaxis.set_minor_formatter(mdates.DateFormatter('%b'))
            plt.setp(ax.xaxis.get_majorticklabels(), rotation=70)
            plt.setp(ax.xaxis.get_minorticklabels(), rotation=70)

    fig.tight_layout()
    fig.savefig('relativeSingularityProfile.pdf')
    sup_fig.tight_layout()
    sup_fig.savefig('relativeSingularityProfile_suplementary.pdf')

    plt.show()


def main(base_file_name):
    dates = parse_fasta(base_file_name)
    file_pattern = f"{base_file_name}-k*.eg"
    files = glob.glob(file_pattern)

    all_data = {}
    for file_path in files:
        k = int(file_path.split('-k')[-1].split('.')[0])  # Extract k value and convert to int
        result = parse_file(file_path)
        if result is not None:
            all_data[k] = result

    plot_data(all_data, dates)

if __name__ == "__main__":
    if len(sys.argv) > 1:
        base_file_name = sys.argv[1]
        main(base_file_name)
    else:
        print("Please provide the base file name and the FASTA file.")

