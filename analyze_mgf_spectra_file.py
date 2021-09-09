import argparse

import numpy as np
from pyteomics import mgf
import pandas as pd

class mgf_parser:

    def __init__(self, spectra_file, output):
        self.spectra_file = spectra_file
        self.low_quality_spectra = {}
        self.output = output

    def read_all_information_from_mgf(self):
        """
        id
        protein [list]
        -> peptide
        -> -> hyperscore
        :return:
        """
        indexed_mgf_obj = mgf.read(source=self.spectra_file,
                                  iterative=True
                                  )
        return indexed_mgf_obj

    def remove_peaks(self, peaks, mz_min=200, dynamic_range=100):
        peaks = np.delete(peaks, np.where(peaks < mz_min))
        if len(peaks) > 4:
            max_peak = np.amax(peaks)
            peaks = peaks/max_peak*dynamic_range
            peaks = np.delete(peaks, np.where(peaks < 1))
        return peaks

    def read_mgf_iterative(self):
        # xtandem default spectrum parameter:
        # spectrum, minimum fragment mz Mfmin = 200 -> peaks caused by individual amino acid residues tend to be relatively small (m/z < 200)
        # spectrum, minimum peaks = 5 -> screen out spectra that contain too few fragment ions to be usefully interpreted.
        # spectrum, minimum parent m+h (mass + proton) = 850 ->  supress the analysis of spectra that were generated by low mass parent ions
        # -> This parameter is not used if the value of spectrum, use noise suppression = no = default
        # spectrum, dynamic range = 100 -> threshold peaks, highest normalized intensity values to 100, all others linear scaled down, if peak<1 rejected
        # {'params':
        # {title 'Run1_U1_2000ng.5651.5651.4 File:"Run1_U1_2000ng.raw", NativeID:"controllerType=0 controllerNumber=1 scan=5651"',
        # 'rtinseconds': int, 'pepmass': (423.71 (mass), 694374.1875(intensity)), 'charge': [int] },
        # 'm/z array': [np array], 'intensity array': [np array], 'charged array':[masked array]
        # searchgui xtandem param  "dynamicRange": 100.0,
        #           "nPeaks": 50,
        #           "minPrecursorMass": 500.0,
        #           "minFragmentMz": 200.0,
        #           "minPeaksPerSpectrum": 5,
        #           "useNoiseSuppression": false,
        # default settings: passing spectra = 163113
        # without charge one: passing spectra = 88666
        reader = mgf.MGF(source=self.spectra_file, use_header=True, convert_arrays=2, read_charges=True, dtype=None,
                                encoding='utf-8', read_ions=False)
        removed_spectra= set()
        number_of_spectra=0
        for spectrum in reader:
            number_of_spectra+=1
            charge = 1 if 'charge' not in spectrum['params'].keys() else int(spectrum['params']['charge'][0])
            if charge == 1:
                removed_spectra.add(spectrum['params']['title'])
            elif len(spectrum['m/z array']) < 5:
                removed_spectra.add(spectrum['params']['title'])
            elif charge > 4:
                removed_spectra.add(spectrum['params']['title'])
            # elif int(spectrum['params']['pepmass'][0]) < 500: # not used since "useNoiseSuppression": false,
            #     removed_spectra.append(spectrum['params']['title']) # passing spectra:
            else:
                intensity_peaks = self.remove_peaks(spectrum['intensity array'], mz_min=200, dynamic_range=100)
                if len(intensity_peaks) < 5:
                    removed_spectra.add(spectrum['params']['title'])
            # try:
            #     if int(spectrum['params']['pepmass'][0]) < 850:
            #         removed_spectra.append(spectrum['params']['title'])
            # except TypeError: # None
            #     pass
        print(f"Number of spectra: {number_of_spectra}")
        print(f"Number of spectra not passing quality control: {len(removed_spectra)}")
        print(f"Number of spectra passing quality control: {number_of_spectra-len(removed_spectra)}")

    def read_common_parameters_from_mgfr(self):
        """
        id
        protein [list]
        -> peptide
        -> -> hyperscore
        :return:
        """
        mgf_dict = mgf.read_header(source=self.spectra_file,
                                  # iterative=True,
                                  # huge_tree=True
                                  )
        return mgf_dict

    def write_tsv(self):
        """
        old tsv format by java parser: #SpecFile, Title, Peptide, Protein, Hyperscore, Evalue
        :return: tsv Title | Peptide | Protein | Hyperscore
        """
        df = self.low_quality_spectra
        df.to_csv(str(self.output)+'.tsv', sep='\t')


def main():
    parser = argparse.ArgumentParser(description='Read xtandem xml to .tsv')
    parser.add_argument('-i', '--input', dest='input', default=None, help='spectra file mgf')
    options = parser.parse_args()
    options.input = "/home/jules/Documents/Tax2Proteome/benchmarking/spectra/Run1_U1_2000ng.mgf"
    obj = mgf_parser(options.input , '')
    dict_o = obj.read_mgf_iterative()
    # print(dict_o)


if __name__ == '__main__':
    main()