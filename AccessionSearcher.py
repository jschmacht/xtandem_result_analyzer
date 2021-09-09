from pathlib import Path
import pickle


class AccessionSearcherNCBI:
    def __init__(self, path_to_file):
        """
        :param path_to_file: path to multiaccession file generated from NCBI-nr database
        one line one header with all accessions seperated by \t
        """
        self.path_to_file = Path(path_to_file)
        self.path_to_acc_dict = Path(str(self.path_to_file) +'_dict')
        self.acc_pos_dict = {}

    
    def search_ncbi_accession(self, accessions):
        """
        :param accessions: kist of lists of accessions
        """
        # build or load acc-pos dict: accession: (pos line start, pos line end)
        if not self.path_to_acc_dict.is_file():
            print('Start bulding accession dict (first accession:(position line start, position line end).')
            acc_pos_dict = {}
            f = open(str(self.path_to_file), 'r')
            while True:
                start = f.tell()
                line = f.readline()
                end = f.tell() - start
                if not line:
                    f.close()
                    break
                acc_pos_dict[line.split('\t')[0]] = (start, end)
            print('Save accession_position dict to %s' % open(str(self.path_to_acc_dict)))
            pickle.dump(self.acc_pos_dict, open(str(self.path_to_acc_dict), 'wb'))
        # load acc_pos_list
        else:
            try:
                print('Load accession_position list.')
                self.acc_pos_dict = pickle.load(open(str(self.path_to_acc_dict), 'rb'))
                print('Loaded')
            except UnicodeDecodeError or EOFError:
                print('Failed opening path to accessin_position list.')
                exit(1)

        # search accessions
        f = open(str(self.path_to_file))
        for acc in [item for sublist in accessions for item in sublist]:
            # search positions of start point line with accessions of interest
            try:
                chunk = self.acc_pos_dict[acc]
                f.seek(chunk[0])
                line = f.read(chunk[1]).decode("utf-8")
                self.acc_dict[acc] = line.split('\t')
            except ValueError:
                continue
        f.close()