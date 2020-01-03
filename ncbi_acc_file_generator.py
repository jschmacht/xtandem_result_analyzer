from pathlib import Path
from tqdm import tqdm
import argparse

# generate file consisting of all accession IDs to one fasta entry of NCBI-nr database
def main():
    parser = argparse.ArgumentParser(description='Read ncbi nr')
    parser.add_argument('-i', '--input', dest='input', default=None, help='path to nr database')

    options = parser.parse_args()

    out_path = Path(options.input).parents[0] / 'multispecies_acc'
    with open(str(out_path), 'w') as out:
        with open(options.input, 'r') as input:
            for line in tqdm(iter(input.readline, ''), total=1333693751):
                line = line.rstrip()
                if line.startswith('>'):
                    if '\x01' in line:
                        out.write('\t'.join([hdr.split(' ')[0] for hdr in line[1:].split('\x01')])+'\n')


if __name__ == '__main__':
    main()