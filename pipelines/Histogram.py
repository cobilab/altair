import matplotlib.pyplot as plt
import matplotlib.ticker as ticker

def read_fasta(file_path):
    with open(file_path, 'r') as file:
        sequence_lengths = []
        current_seq = ''
        for line in file:
            if line.startswith('>'):
                if current_seq:
                    sequence_lengths.append(len(current_seq))
                    current_seq = ''
            else:
                current_seq += line.strip()
        if current_seq:
            sequence_lengths.append(len(current_seq))
    return sequence_lengths

def plot_histogram(sequence_lengths, output_file, min_length=None, max_length=None):
    plt.figure(figsize=(10, 6))
    plt.hist(sequence_lengths, bins=10000, color='blue', edgecolor='black')

    # Customizing the plot to resemble the gnuplot style
    plt.title('Histogram of Sequence Lengths', color='black', fontsize=12)
    plt.xlabel('Sequence Length', fontsize=12)
    plt.ylabel('Frequency', fontsize=12)

    # Set grid and border styles
    plt.grid(True, which='major', linestyle='-', linewidth='0.5', color='grey')
    plt.minorticks_on()
    plt.grid(True, which='minor', linestyle=':', linewidth='0.5', color='grey')

    # Set range for x-axis and y-axis if provided
    if min_length is not None and max_length is not None:
        plt.xlim(min_length, max_length)

    # Add an arrow annotation if needed
    plt.annotate('', xy=(29903, 0), xytext=(29903, max(sequence_lengths)/3),  
                 arrowprops=dict(facecolor='#CC0000', shrink=0.05, lw=2), fontsize=12)

    # Save the plot as a PDF
    plt.savefig(output_file, format='pdf')

# Replace 'SARS-CoV-2.fa' with the path to your FASTA file
fasta_file = 'SARS-CoV-2.fa'  # Example file path
histogram_file = 'Histogram.pdf'  # Output file

# Read the FASTA file and plot the histogram
sequence_lengths = read_fasta(fasta_file)
plot_histogram(sequence_lengths, histogram_file, min_length=29500, max_length=30000)
