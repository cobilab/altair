import unittest
from datetime import datetime
import argparse 

def parse_fasta_file(file_content):
    sections = file_content.split('>')[1:]
    parsed_data = []

    for section in sections:
        header, *body = section.strip().split('\n')
        country, date_str = header.split('|')
        date_str = date_str.strip()
        date_format = "%Y-%m-%d" if len(date_str) > 7 else "%Y-%m" if len(date_str) > 4 else "%Y"
        date = datetime.strptime(date_str, date_format)
        standardized_date_str = date.strftime("%Y-%m-%d")
        standardized_header = f"> {country.strip()} | {standardized_date_str}"
        parsed_data.append((standardized_header, body))

    return parsed_data

def sort_and_reconstruct_fasta(parsed_data):
    sorted_data = sorted(parsed_data, key=lambda x: datetime.strptime(x[0].split('|')[1].strip(), "%Y-%m-%d"))
    sorted_fasta = '\n'.join(['\n'.join([header] + body) for header, body in sorted_data])
    return sorted_fasta

def process_fasta_file(input_path, output_path):
    with open(input_path, 'r') as file:
        file_content = file.read()

    parsed_data = parse_fasta_file(file_content)
    sorted_fasta_content = sort_and_reconstruct_fasta(parsed_data)

    with open(output_path, 'w') as sorted_file:
        sorted_file.write(sorted_fasta_content)

def is_file_sorted(file_path):
    with open(file_path, 'r') as file:
        file_content = file.read()

    parsed_data = parse_fasta_file(file_content)
    dates = [header.split('|')[1].strip() for header, _ in parsed_data]

    for i in range(len(dates) - 1):
        current_date = datetime.strptime(dates[i], "%Y-%m-%d")
        next_date = datetime.strptime(dates[i + 1], "%Y-%m-%d")
        if current_date > next_date:
            print(f"Sorting error between:\n{parsed_data[i][0]}\nand\n{parsed_data[i + 1][0]}")
            return False
    return True



# Unit Tests
class TestFastaSorting(unittest.TestCase):

    def test_date_standardization(self):
        sample_data = "> Country |2022\nATCG\n> Country |2021-06\nATCG\n> Country |2020-12-01\nATCG"
        parsed_data = parse_fasta_file(sample_data)
        for header, _ in parsed_data:
            _, date_str = header.split('|')
            self.assertEqual(len(date_str.strip()), 10)  # Check if date format is YYYY-MM-DD

    def test_sorting_order(self):
        sample_data = "> Country |2022\nATCG\n> Country |2021-06\nATCG\n> Country |2020-12-01\nATCG"
        parsed_data = parse_fasta_file(sample_data)
        sorted_data = sort_and_reconstruct_fasta(parsed_data)
        expected_order = ["> Country | 2020-12-01", "> Country | 2021-06-01", "> Country | 2022-01-01"]
        sorted_headers = [line for line in sorted_data.split('\n') if line.startswith('>')]
        self.assertEqual(sorted_headers, expected_order)
    
    def test_output_file_sorting(self):
        process_fasta_file('FIL-SARS-CoV-2.fa', 'sorted_output.fa')
        self.assertTrue(is_file_sorted('sorted_output.fa'))

if __name__ == '__main__':
    # Create an argument parser
    parser = argparse.ArgumentParser(description='Process a FASTA file.')
    parser.add_argument('input_path', type=str, help='Path to the input FASTA file')
    parser.add_argument('output_path', type=str, help='Path to the output FASTA file')
    
    # Parse the arguments
    args = parser.parse_args()

    # Process the FASTA file with the provided paths
    process_fasta_file(args.input_path, args.output_path)
    
    # Run unit tests
    unittest.main(argv=['first-arg-is-ignored'], exit=False)