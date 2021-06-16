import re

species_to_tax_dict = {'Enterobacteria phage ES18': 101570,
                       'Enterobacteria phage P22':  10754,
                       'Salmonella phage Felix O1 (isolate Felix O1-VT1)': 1283336,
                       'Enterobacteria phage M13': 1977402,
                       'Enterobacteria phage MS2': 12022,
                       'Nitrosospira multiformis (strain ATCC 25196 / NCIMB 11849)': 323848,
                       'Nitrososphaera viennensis EN76 NV_AIC17288.1': 926571,
                       'Nitrosomonas ureae': 44577,
                       'Nitrosomonas europaea (strain ATCC 19718 / NBRC 14298)': 228410,
                       'Escherichia coli (strain K12)': 83333,
                       'Agrobacterium fabrum (strain C58 / ATCC 33970)': 176299,
                       'Desulfovibrio vulgaris (strain Hildenborough / ATCC 29579 / NCIMB 8303)': 882,
                       'Thermus thermophilus (strain HB27 / ATCC BAA-163 / DSM 7039)': 262724,
                       'Bacillus subtilis (strain 168)': 224308,
                       'Rhizobium leguminosarum bv. viciae (strain 3841)': 216596,
                       'Alteromonas macleodii (strain English Channel 673)': 1004788,
                       'Burkholderia xenovorans (strain LB400)': 266265,
                       'Cupriavidus metallidurans (strain ATCC 43123 / DSM 2839 / NBRC 102507 / CH34)': 266264,
                       'Salmonella typhimurium (strain LT2 / SGSC1412 / ATCC 700720)': 99287,
                       'Pseudomonas denitrificans ATCC 13867': 1294143}
species_to_tax_dict_2 = {'Pseudomonas fluorescens ATCC 13525': 294,
                         'Pseudomonas pseudoalcaligenes KF707': 1149133,
                         'Chlamydomonas reinhardtii': 3055,
                         'Staphylococcus aureus subsp. aureus': 1280,
                         'Nitrososphaera viennensis EN76': 926571,
                         'Rhizobium leguminosarum': 384,
                         'Stenotrophomonas maltophilia': 40324,
                         'Stenotrophomonas maltophilia SeITE02': 1407502}
species_to_tax_dict_metagenome_seq = {'Chromobacterium violaceum CV026': 536, # parent: 536
                         'Staphylococcus aureus ATCC 13709': 1280, #Staphylococcus aureus subsp. aureus NCTC 8325: 93061, parent: 1280
                         'Paracoccus denitrificans JCM 21484': 1302247,
                        'Roseobacter sp. AK199': 1907202} # species: Roseobacter sp.
multispecies_to_species_dict ={'Stenotrophomonas' : 'Stenotrophomonas maltophilia SeITE02',
                               'Proteobacteria': 'Stenotrophomonas maltophilia SeITE02', # this might be uncorrect
                               'Xanthomonadaceae': 'Stenotrophomonas maltophilia SeITE02',
                               'Rhizobium': 'Rhizobium leguminosarum',
                               'Rhizobiaceae': 'Rhizobium leguminosarum',
                               'Rhizobiales': 'Rhizobium leguminosarum',
                               'Rhizobium/Agrobacterium group': 'Rhizobium leguminosarum'}
uniprot_set=set()
acc_set=set()
ncbi_acc_set=set()
with open('/home/jules/Documents/Tax2Proteome/benchmarking/Kleiner_ref_db/Mock_Comm_RefDB_V3.fasta', 'r') as input:
    with open('/home/jules/Documents/Tax2Proteome/benchmarking/Kleiner_ref_db/Kleiner_ref_crap.fasta', 'w') as output:
        with open('/home/jules/Documents/Tax2Proteome/benchmarking/Kleiner_ref_db/acc2tax_custom', 'w') as tax_output:
            for line in input:
                line=line.strip()
                if not line:
                    continue
                if line.startswith('>'):
                    acc0=line.split()[0][1:]

                    # CRAP
                    if line.split()[0].startswith('>CRAP'):
                        acc = line.split()[0][1:]
                        tax_output.write('CRAP'+ '\t' + acc.strip() + '\t' + acc + '\t' + 'crap' + '\t' + '0'+ '\n')
                    # metagenome
                    elif re.search(r'_peg.', line.split()[0]):
                        acc = line.split()[0][1:]
                        if acc.startswith('137'):
                            species='Staphylococcus aureus ATCC 13709'
                        elif acc.startswith('CV'):
                            species = 'Chromobacterium violaceum CV026'
                        elif acc.startswith('PaD'):
                            species = 'Paracoccus denitrificans JCM 21484'
                        elif acc.startswith('AK199'):
                            species = 'Paracoccus denitrificans JCM 21484'
                        tax_output.write('meta_seq'+ '\t' + acc.strip() + '\t' + acc + '\t' + species + '\t' + str(species_to_tax_dict_metagenome_seq[species])+ '\n')
                    # uniprot_acc
                    elif re.search(r'OS=.+ GN=', line):
                        species = re.search(r'OS=(.+) GN=', line).group(1)
                        acc=line.split()[1]
                        tax_output.write('uniprot'+ '\t' + acc0.strip()  + '\t' + acc.strip()  + '\t' + species + '\t' + str(species_to_tax_dict[species])+ '\n')
                    elif re.search(r'OS=.+ PE=', line):
                        species = re.search(r'OS=(.+) PE=', line).group(1)
                        acc=line.split()[1]
                        tax_output.write('uniprot'+ '\t' + acc0.strip()  + '\t' + acc.strip()  + '\t' + species + '\t' + str(species_to_tax_dict[species])+ '\n')
                    elif '_WP_' in line.split()[0]:
                        acc = re.search(r'>.*_(WP_.*)', line.split()[0]).group(1)
                        try:
                            species = re.search(r'\[(.*)\]', line).group(1)
                            if ']' in species:
                                species = re.search(r'\[.*\].*\[(.*)\]', line).group(1)
                            if 'MULTISPECIES' in line:
                                species = multispecies_to_species_dict[species]
                            tax_output.write('ncbi_multispecies' + '\t' + acc.strip()  + '\t' + acc.strip()  + '\t' + species + '\t' + str(species_to_tax_dict_2[species])+ '\n')
                        except KeyError:
                            # entries of type >SMS_WP_001878988.1 MULTISPECIES: hypothetical protein [Bacteria], not longer annotated to any genome
                            tax_output.write('outdated' + '\t' + acc.strip()  + '\t' + acc.strip()  + '\t' + 'no_species' + '\t' + str(0)+ '\n')
                            pass
                    elif '_NP_' in line.split()[0]:
                        acc = re.search(r'>.*_(NP_.*)', line.split()[0]).group(1)
                        tax_output.write('ncbi' + '\t' + acc0.strip()  + '\t' + acc.strip() + 'unknown_species' + '\t' + str(0)+ '\n')
                  #  elif 'SEED' in line:
                   #     acc = line.split()[0][1:]
                    #    tax_output.write('custom' + '\t' + acc0.strip() + '\t' + acc.strip() + 'unknown_species' + '\t' + str(0)+ '\n')
                    elif len(line.split()) == 1:
                        acc = line[1:]
                        tax_output.write('custom' + '\t' + acc0.strip() + '\t' + acc.strip() + '\t' + 'unknown_species' + '\t' + str(0)+ '\n')
                    elif re.search(r'\[(.+)\]', line):
                        species = re.search(r'\[(.+)\]', line).group(1)
                        if '[' in species:
                            #print(species)
                            species = species.split('[')[-1]
                        if '[' in species:
                            #print(species)
                            species = species.split('[')[1]
                        if ' : ' in species:
                            species = species.split(' : ')[0]
                        try:
                            acc = line.split()[0][1:]
                            if '_XP_' in line.split()[0]:
                                tax_output.write('ncbi' + '\t' + acc0.strip()  + '\t' + acc.strip()  + '\t' + species + '\t' + str(species_to_tax_dict_2[species])+ '\n')
                            else:
                                tax_output.write('unknown' + '\t' + acc0.strip()  + '\t' + acc.strip()  + '\t' + species + '\t' + str(species_to_tax_dict_2[species])+ '\n')
                        except KeyError:
                            print(line)
                    else:
                        print('line:', line)
                    output.write('>' + acc0 + '\n')
                else:
                    output.write(line + '\n')
