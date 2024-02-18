import subprocess
import pandas as pd
import re
import os
from Bio.Seq import Seq
from itertools import combinations, groupby

# Global variables:
num_of_nuc_around_mut = 20
mutation_size = 6
size_of_seq = 2*num_of_nuc_around_mut + 2*26 + 1
mut_num = 250 #max 250.  0 -> to get the common above "threshold"%
mer_lenght = 9
list_of_100_HLA = ['HLA-A0201', 'HLA-C0701', 'HLA-A0101', 'HLA-C0702', 'HLA-C0401', 'HLA-A0301', 'HLA-B0702', 'HLA-A2402', 'HLA-B0801', 'HLA-C0602', 'HLA-B4402', 'HLA-C0501', 'HLA-C0304', 'HLA-A1101', 'HLA-B3501', 'HLA-C1203', 'HLA-B5101', 'HLA-B1501', 'HLA-B4403', 'HLA-B1801', 'HLA-C0303', 'HLA-B4001', 'HLA-C0102', 'HLA-C0202', 'HLA-C1601', 'HLA-C0802', 'HLA-A2601', 'HLA-A3201', 'HLA-A6801', 'HLA-B2705', 'HLA-A2301', 'HLA-A2902', 'HLA-B5701', 'HLA-A3101', 'HLA-B1402', 'HLA-B1302', 'HLA-C1502', 'HLA-B3801', 'HLA-A2501', 'HLA-A3001', 'HLA-B4901', 'HLA-B3503', 'HLA-A3303', 'HLA-C1402',
      'HLA-C1701', 'HLA-B5201', 'HLA-C0704', 'HLA-B5501', 'HLA-B4002', 'HLA-A3002', 'HLA-C1202', 'HLA-B5801', 'HLA-B5301', 'HLA-B3901', 'HLA-A6802', 'HLA-B3701', 'HLA-B3502', 'HLA-B5001', 'HLA-A3301', 'HLA-B4501', 'HLA-A0205', 'HLA-B1401', 'HLA-C0801', 'HLA-C0302', 'HLA-C0210', 'HLA-B1503', 'HLA-C1505', 'HLA-B5601', 'HLA-A7401', 'HLA-B0705', 'HLA-A2901', 'HLA-A0206', 'HLA-B4601', 'HLA-B4101', 'HLA-B4201', 'HLA-B1502', 'HLA-B2702', 'HLA-B3906', 'HLA-B4102', 'HLA-B5802', 'HLA-A6601', 'HLA-A0202', 'HLA-B3508', 'HLA-C1602', 'HLA-B1510', 'HLA-A0207', 'HLA-B5703', 'HLA-B1517',
        'HLA-A3402', 'HLA-A0203', 'HLA-A3601', 'HLA-C1801', 'HLA-B4405', 'HLA-A3004', 'HLA-B3802', 'HLA-B1516', 'HLA-B1301', 'HLA-A0302', 'HLA-A2403', 'HLA-A0325']
# list_of_HLA = ['HLA-A0201', 'HLA-C0701', 'HLA-A0101', 'HLA-C0702', 'HLA-C0401', 'HLA-A0301', 'HLA-B0702', 'HLA-A2402', 'HLA-B0801', 'HLA-C0602']
list_of_HLA = ['HLA-A0201', 'HLA-C0701', 'HLA-A0101', 'HLA-C0702']

check_flag = True
check_OUT = []
double_editing = True
single_editing = False

def most_common_mut(mut_num, primary_site):
    f_path = f"/home/alu/netlandes/MHCpan/all_cancers/{primary_site}.tsv"  #frequent-mutations
    # Read the TSV file into a pandas DataFrame
    df = pd.read_csv(f_path, delimiter='\t')
    
    # Filter rows where the third column does not start with "Stop Gained"
    filtered_df = df[~df.iloc[:, 2].str.startswith("Stop Gained")]

    if mut_num == 0:
        # Extract the fourth column from the DataFrame
        filtered_df[['Value1', 'Value2']] = filtered_df.iloc[:, 3].str.split('/', expand=True)
        filtered_df['Value2'] = filtered_df['Value2'].str.rsplit(',', n=1, expand=True)[1].str.replace('%', '').astype(float)
        if primary_site == "ALL":
            threshold = 1
        else:
            threshold = 5
        filtered_df2 = filtered_df[filtered_df['Value2'] >= threshold]              
        selected_rows = filtered_df2.iloc[:,0:4].values
    else:
        if mut_num > 250:
            mut_num = 250
        selected_rows = filtered_df.iloc[:mut_num, 0:4].values

    # Create a list of rows
    original_list = selected_rows.tolist()
 

    new_list = []

    for item in original_list:
        original_string = item[0]
        gene_mutation = item[2].split()[-2:]  # Get the last two words
        gene_name, mutation = gene_mutation
        precent = item[3].split(',')[1]                                      
        if primary_site == "ALL":
            precent = item[3].split(',')[2]
        new_list.append([original_string, gene_name, mutation, precent])

    return new_list

def delete_files(folder_path):
    # iterate over all files in the folder
    for filename in os.listdir(folder_path):
        # construct the full file path
        file_path = os.path.join(folder_path, filename)
        
        # check if the file path is a regular file (not a folder or a symlink)
        if os.path.isfile(file_path):
            # delete the file
            os.remove(file_path)

def checking_input(mutation_input):
    flag_indel = False
    # We want deletion and insertion mutations only:
    if 'del' in mutation_input:
        flag_indel = True
        pattern = r'(\w+):g\.(\d+)([a-z]+)([A-Z]+)'
    elif 'ins' in mutation_input:
        flag_indel = True
        pattern = r'(\w+):g\.(\d+)_\d+([a-z]+)([A-Z]+)'
    elif '>' in mutation_input:
        pattern = r'(\w+):g\.(\d+)(\w)>(\w)'
    else:
        return None, flag_indel
    return re.search(pattern, mutation_input), flag_indel

def run_MHCpan(seq_list, patient_HLA, num_folder):
    count = 0
    a = ">Temp:"
    for seq in seq_list:
        count += 1
        mut_aa_seq = Seq(seq).translate()
        with open("/home/alu/netlandes/MHCpan/input.fsa", 'w') as f:
            f.write(a + '\n')
            f.write(str(mut_aa_seq))
            output_file = f"/home/alu/netlandes/MHCpan/results_common{num_folder}/result{count}.xls"  
        subprocess.run(f'"/private/common/Software/netMHC-4.0/Linux_x86_64/bin/netMHC" -f "/home/alu/netlandes/MHCpan/input.fsa" -tdir "/tmp" -version "/private/common/Software/netMHC-4.0/Linux_x86_64/data/version " -hlalist "/private/common/Software/netMHC-4.0/data/allelelist" -syn "/private/common/Software/netMHC-4.0/Linux_x86_64/data/synlists/%s.synlist" -thrfmt "/private/common/Software/netMHC-4.0/threshold/%s.thr" -rdir "/private/common/Software/netMHC-4.0/Linux_x86_64" -xlsfile {output_file} -xls -a {patient_HLA}', shell=True)

def create_results_file(tuples_of_SB, tuples_of_WB, hla):
    dir_files = os.listdir("/home/alu/netlandes/MHCpan/results_common1/")
    i = 0
    for file in dir_files:
        path = "/home/alu/netlandes/MHCpan/results_common1/" + file
        if os.path.getsize(path) == 0:
            i += 1
            continue
        results = pd.read_csv(path, sep='\t', skiprows=1)
        new_results = results.loc[results['N_binders'] == 1]
        if not new_results.empty:
            for index, row in new_results.iterrows():
                if check_flag:
                    check_results = pd.read_csv("/home/alu/netlandes/MHCpan/results_common2/" + file, sep='\t', skiprows=1)
                    check_new_results = check_results.loc[check_results['N_binders'] == 1]
                    if not check_new_results.empty and (index in check_new_results.index) and (check_new_results.at[index, 'Rank'] <= 0.5 or check_new_results.at[index, 'Rank'] <= row['Rank']):
                        check_OUT.append((check_new_results.at[index, 'Peptide'], check_new_results.at[index, 'Rank'] , row['Peptide'], row['Rank']))
                        continue
                if(row['Rank'] <= 0.5):
                    tuples_of_SB.append((i , row['Peptide'], row['Rank'], hla))
                else:
                    tuples_of_WB.append((i , row['Peptide'], row['Rank'], hla))
        i += 1


if __name__ == "__main__":
    if single_editing:
        file_name = "final_results_SINGLE"
    else:
        file_name = "final_results_DOUBLE"
    with open (f"/home/alu/netlandes/MHCpan/{file_name}.tsv", 'w') as final_file:
        final_file.write(f"Gene-Name:\tAA-change:\tMutation:\tStrand:\tPrevalence-percentage:\tTarget:\tGuide-RNA:\tRank:\tHLA:\n")
        file_names  = os.listdir("/home/alu/netlandes/MHCpan/all_cancers")  #frequent-mutations
        for primary_site in [os.path.splitext(file_name)[0] for file_name in file_names]:
            list_of_mutation = most_common_mut(mut_num, primary_site)
            for mut in list_of_mutation:
                # mut = ['chr10:g.87958013delA', 'PIK3CA', 'E545K', '9.80%']
                mutation_input = mut[0]
                match, flag_indel = checking_input(mutation_input)
                pos_list = []
                if match == None: # A mutation that is not a deletion insertion Or Substitution
                    continue
                mut_chr = match.group(1)
                mut_pos = match.group(2)
                if '>' in mutation_input:
                    mut_type = "Substitution"
                else:
                    mut_type = match.group(3)
                mut_seq = match.group(4)

                # We want the size of the insertion or deletion to be no more than 6bp:
                if (len(mut_seq) > mutation_size):
                    print(f"The mutation should be a maximum {mutation_size}bp\n")
                    continue

                # Run a program that accepts as input a chromosome and location and returns the 60 bases around it:
                p = subprocess.run(f'sh "/home/alu/netlandes/MHCpan/find_seq_of_mutation.sh" {mut_chr} {mut_pos} {int(size_of_seq / 2)}', capture_output=True,text=True, shell=True)
                seq = p.stdout.replace("\n", "")

                # We will run a program that will return the reading frame of the sequence:
                subprocess.run(f'sh "find.sh" {mut_chr} {mut_pos}', shell=True)
                output_file = "/home/alu/netlandes/MHCpan/intersects_output2.bed"
                if os.stat(output_file).st_size == 0:
                    print("The mutation is not in CDS\n")
                    continue
                temp = seq[::-1]
                isoforms_files = pd.read_csv(output_file, sep="\t", header=None)
                list_of_codons = []
                frame_init_pos = int(mut_pos) - int((size_of_seq/2) - 1)
                for index, row in isoforms_files.iterrows():
                    strand = row[4]
                    reading_frame = row[5]
                    if strand == "+":
                        cds_start_pos = row[1]
                        cds_end_pos = row[2]
                        count = 0
                        count2 =1
                    else:
                        count = 1
                        count2 =0
                        seq = temp
                        frame_init_pos = int(mut_pos) + int((size_of_seq/2) - 1)
                        cds_start_pos = row[2]
                        cds_end_pos = row[1]
                    end = len(seq)
                    flag = True
                    if abs(int(mut_pos) - cds_start_pos) < int((size_of_seq/2) - 1): 
                        start = int((size_of_seq/2) - 1) - abs(int(mut_pos) - cds_start_pos) + reading_frame + count
                        frame_init_pos = cds_start_pos + reading_frame 
                        flag = False
                    if abs(cds_end_pos - int(mut_pos)) < int(size_of_seq/2):
                        end = abs(cds_end_pos - frame_init_pos)
                        if not flag:
                            end += start
                    if flag:
                        if reading_frame == 0:
                            seq_frame = abs(frame_init_pos - cds_start_pos + count) % 3 
                        elif reading_frame == 1:
                            seq_frame = abs(frame_init_pos - 1 - cds_start_pos) % 3 
                        else:
                            seq_frame = abs(frame_init_pos - 2 - cds_start_pos - count) % 3
                        if seq_frame == 0:
                            start = seq_frame 
                        elif seq_frame == 1:
                            start = seq_frame + 1
                        elif seq_frame == 2:
                            start = seq_frame - 1
                    codons = [seq[i:i+3] for i in range(start, end, 3)]
                    
                    # Remove duplicates:
                    if codons not in list_of_codons:
                        list_of_codons.append(codons)
                # print(codons)

                # Creats mutatation in the seq for each isoform:
                mut_isoforms = []
                pos = int(size_of_seq/2) - start - count2 + count
                if strand == "-":
                    if mut_type == "ins":
                        mut_seq = mut_seq[::-1]
                    pos = pos - len(mut_seq)
 
                for codon in list_of_codons:
                    # Convert the list of codons to a DNA sequence
                    dna_seq = ''.join(codon)
                    if mut_type == "del":
                        # Remove the sequence from the DNA sequence:
                        mut_dna_seq = dna_seq[:pos] + dna_seq[pos + len(mut_seq):]
                    elif mut_type == "ins":
                        pos = int(size_of_seq/2) - start 
                        # Insert the sequence to the DNA sequence:
                        mut_dna_seq = dna_seq[:pos] + mut_seq + dna_seq[pos:]
                        for i in range(len(mut_seq)):
                            pos_list.append(pos + i)
                    else:
                        mut_dna_seq = dna_seq[:pos] + mut_seq + dna_seq[pos + 1:]
                        pos_list.append(pos)

                    # The netMHC tool accepts a peptide of at least 9mer in length
                    if len(mut_dna_seq) < mer_lenght * 3:
                        continue

                    # If the mutation is on the negative strand, we will look for editing on the negative strand:
                    if strand == "-":
                        mut_dna_seq = str(Seq(mut_dna_seq).complement())

                    mut_isoforms.append(mut_dna_seq)

                m = len(mut_seq)
                if mut_type == "del":
                    m = 0
        
                ADAR_mut_seq_list = []

                if single_editing:
                    for isoform in mut_isoforms:
                        for i in range(max(pos - num_of_nuc_around_mut + 1, 1), min(pos + num_of_nuc_around_mut + m - 1, len(isoform) - 1)):
                            sequence_list = list(isoform)
                            if sequence_list[i] == "A":                         
                                # We will not edit adenosines that are in the mutant sequence:
                                if i in pos_list:
                                    continue
                                new_string = sequence_list
                                new_string[i] = "G"
                                temp =isoform[:i] + "*" + isoform[i + 1:]
                                gRNA = temp[max(pos - num_of_nuc_around_mut, 0) : min(pos + num_of_nuc_around_mut + m, len(sequence_list))]
                                if i % 3 == 0:
                                    start_index = max(i - 24, 0)
                                    end_index = i + 27
                                elif i % 3 == 1:
                                    start_index = max(i - 25, 0)
                                    end_index = i + 26
                                else:
                                    start_index = max(i - 26, 0)
                                    end_index = i + 25
                                ADAR_mut_seq_list.append(("".join(new_string[start_index : end_index]), gRNA, isoform[start_index : end_index]))

                if double_editing:
                    for isoform in mut_isoforms:
                        gRNA_editing_sites = []
                        sequence_list = list(isoform)
                        for i in range(max(pos - num_of_nuc_around_mut + 1, 1), min(pos + num_of_nuc_around_mut + m - 1, len(isoform) - 1)):
                            if sequence_list[i] == "A":                         
                                # We will not edit adenosines that are in the mutant sequence:
                                if i in pos_list:
                                    continue
                                gRNA_editing_sites.append(i)

                        # print(gRNA_editing_sites)
                        unique_pairs = list(combinations(gRNA_editing_sites, 2))
                        for pair in unique_pairs:
                            sequence_list = list(isoform)
                            a = pair[1] - pair[0] + pair[0] % 3 + 1  # The distance in nucleotides
                            distance = -(-a // 3)  # The distance in amino acids
                            if distance <= mer_lenght:
                                new_string = sequence_list
                                new_string[pair[0]] = "G"
                                new_string[pair[1]] = "G"
                                temp = "".join(new_string)
                                # gRNA = temp[max(pos - n, 0) : min(pos + n + m, len(sequence_list))]
                                x = mer_lenght - distance
                                groups_of_three = []
                                for i in range(0, len(temp), 3):
                                    group = temp[i:i+3]  
                                    groups_of_three.append(group)
                                first_nuc = pair[0] + 1
                                second_nuc = pair[1] + 1
                                start_index = max(-(-first_nuc // 3) - 1 - x, 0)
                                end_index = -(-second_nuc // 3) + x
                                new_string[pair[0]] = "*"
                                new_string[pair[1]] = "*"
                                temp = "".join(new_string)
                                gRNA = temp[max(pos - num_of_nuc_around_mut, 0) : min(pos + num_of_nuc_around_mut + m, len(sequence_list))]
                                ADAR_mut_seq_list.append(("".join(groups_of_three[start_index : end_index]), gRNA, isoform[start_index*3 : end_index*3]))

                # for adar in ADAR_mut_seq_list:
                #     print(adar[0])
                # print("---------------------------------------")
                # for adar in ADAR_mut_seq_list:
                #     print(adar[2])
                # exit()

                # Remove duplicates from the adar_mut list:
                unique_ADAR_mut_list = list(set(ADAR_mut_seq_list))


                if unique_ADAR_mut_list:
                    sb_tuple = []
                    wb_tuple = []
                    for patient_HLA in list_of_100_HLA:
                        # Running the tool on the edited sequences:
                        run_MHCpan([t[0] for t in unique_ADAR_mut_list], patient_HLA, 1)
                        if check_flag:
                            # Running the tool on the non-edited sequences for comparison:
                            run_MHCpan([t[2] for t in unique_ADAR_mut_list], patient_HLA, 2)
                        create_results_file(sb_tuple, wb_tuple, patient_HLA) 
                
                    delete_files('/home/alu/netlandes/MHCpan/results_common1')

                    if check_flag:
                        delete_files('/home/alu/netlandes/MHCpan/results_common2')                    

                    uniq_sb = list(set(sb_tuple))
                    uniq_wb = list(set(wb_tuple))


                    # For each HLA we will find the result with the lowest rank:
                    sorted_data = sorted(uniq_sb + uniq_wb, key=lambda x: x[3])
                    grouped_data = {key: list(group) for key, group in groupby(sorted_data, key=lambda x: x[3])}
                    sb_and_wb = [min(group, key=lambda x: x[2]) for group in grouped_data.values()]

                    if not sb_and_wb:
                        continue              

                    for sb_wb in sb_and_wb:
                        sb_seq = Seq(unique_ADAR_mut_list[sb_wb[0]][0])
                        protein = sb_seq.translate()
                        if "*" in str(protein):
                            protein = str(protein)
                            protein = protein.replace("*", "X")
                            protein = Seq(protein)
                        pos = protein.find(sb_wb[1])
                        ASO_target = sb_seq[(pos * 3):((pos + len(sb_wb[1])) * 3)]
                        final_file.write(f"{mut[1]}\t{mut[2]}\t{mutation_input}\t{strand}\t{mut[3]}\t{str(ASO_target)}\t{str(unique_ADAR_mut_list[sb_wb[0]][1])}\t{sb_wb[2]}\t{sb_wb[3]}\n")
        
        if check_flag:
            uniq_check_OUT = list(set(check_OUT))
            final_file.write("Original:\tEditing:\n")
            for out in uniq_check_OUT:
                final_file.write(f"{out[0]} {out[1]}\t{out[2]} {out[3]}\n")               
