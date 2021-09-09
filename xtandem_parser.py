import argparse

from pyteomics import tandem
import pandas as pd


class XtandemParser:

    def __init__(self, tandem_result):
        self.tandem_result = tandem_result

    def read_xtandem_xml_iterative(self):
        """
        id
        protein [list]
        -> peptide
        -> -> hyperscore
        :return:
        """
        reader = tandem.TandemXML(source=self.tandem_result,
                                  iterative=True,
                                  # huge_tree=True
                                  )
        dict = reader.read()
        # for protein in reader.read():
        #    print(protein)
        # Protein name
        #print(protein['label'].strip())


    def read_xtandem_xml_at_once(self):
        """
        id
        protein [list]
        -> peptide
        -> -> hyperscore
        :return:
        """
        reader = tandem.TandemXML(source=self.tandem_result,
                                 # iterative=True,
                                 # huge_tree=True
                                  )
        df = tandem.DataFrame(self.tandem_result)
        df = df [['scan', 'seq', 'protein_label', 'hyperscore']]
        df = df.rename(columns={'scan':'Title', 'seq':'Peptide', 'protein_label':'Protein', 'hyperscore':'Hyperscore'})
        df = df.explode('Protein')
        df["Title"].apply(lambda title: title.split("RTINSECONDS")[0].strip())
        pd.set_option('display.max_columns', None)
        print(df.head())
        return df

    def write_tsv(self, df):
        """
        old tsv format by java parser: #SpecFile, Title, Peptide, Protein, Hyperscore, Evalue
        :return: tsv Title | Peptide | Protein | Hyperscore
        """
        df.to_csv(str(self.tandem_result)+'.tsv', sep='\t')


def main():
    parser = argparse.ArgumentParser(description='Read xtandem xml to .tsv')
    parser.add_argument('-i', '--input', dest='input', default=None, help='xtandem output xml')
    options = parser.parse_args()

    xtandem_result = options.input
    obj = XtandemParser(xtandem_result)
    df = obj.read_xtandem_xml_at_once()
    obj.write_tsv(df)


if __name__ == '__main__':
    main()