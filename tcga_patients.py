import subprocess
import pandas as pd
import re
import os
from Bio.Seq import Seq
import gzip
from itertools import combinations


# Global variables:
num_of_nuc_around_mut = 20
mutation_size = 6
size_of_seq = 2*num_of_nuc_around_mut + 2*26 + 1
mer_lenght = 9

edits_list = [0,1,2]
# check_flag = False
# without_editing = False
# single_editing = False
# double_editing = False
# check_OUT = []

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

def checking_HLA(input_HLA):
    with open('/home/alu/netlandes/MHCpan/HLA-list.txt', 'r') as file:
        lines = [line.strip() for line in file.readlines()]
    if input_HLA not in lines:
        print("Look at HLA-list.txt and insert HLA according to the required format")
        exit(1)

def run_MHCpan(seq_list, patient_HLA, num_folder):
    count = 0
    a = ">Temp:"
    for seq in seq_list:
        count += 1
        mut_aa_seq = Seq(seq).translate()
        with open("/home/alu/netlandes/MHCpan/input_seq.fsa", 'w') as f:
            f.write(a + '\n')
            f.write(str(mut_aa_seq))
            output_file = f"/home/alu/netlandes/MHCpan/results{num_folder}/result{count}.xls"  
        subprocess.run(f'"/private/common/Software/netMHC-4.0/Linux_x86_64/bin/netMHC" -f "/home/alu/netlandes/MHCpan/input_seq.fsa" -tdir "/tmp" -version "/private/common/Software/netMHC-4.0/Linux_x86_64/data/version " -hlalist "/private/common/Software/netMHC-4.0/data/allelelist" -syn "/private/common/Software/netMHC-4.0/Linux_x86_64/data/synlists/%s.synlist" -thrfmt "/private/common/Software/netMHC-4.0/threshold/%s.thr" -rdir "/private/common/Software/netMHC-4.0/Linux_x86_64" -xlsfile {output_file} -xls -a {patient_HLA}', shell=True)

def create_results_file(tuples_of_SB, tuples_of_WB, hla):
    dir_files = os.listdir("/home/alu/netlandes/MHCpan/results1/")
    i = 0
    for file in dir_files:
        path = "/home/alu/netlandes/MHCpan/results1/" + file
        if os.path.getsize(path) == 0:
            i += 1
            continue
        results = pd.read_csv(path, sep='\t', skiprows=1)
        new_results = results.loc[results['N_binders'] == 1]
        if not new_results.empty:
            for index, row in new_results.iterrows():
                # if check_flag:
                #     check_results = pd.read_csv("/home/alu/netlandes/MHCpan/results2/" + file, sep='\t', skiprows=1)
                #     check_new_results = check_results.loc[check_results['N_binders'] == 1]
                #     if not check_new_results.empty and (index in check_new_results.index) and (check_new_results.at[index, 'Rank'] <= 0.5 or check_new_results.at[index, 'Rank'] <= row['Rank']):
                #         check_OUT.append((check_new_results.at[index, 'Peptide'], check_new_results.at[index, 'Rank'] , row['Peptide'], row['Rank']))
                #         continue                
                if(row['Rank'] <= 0.5):
                    tuples_of_SB.append((i , row['Peptide'], row['Rank'], hla))
                else:
                    tuples_of_WB.append((i , row['Peptide'], row['Rank'], hla))
        i += 1

def creats_list_of_mutation(mut_path):
    gz_file_path = None
    for file_name in os.listdir(mut_path):
        if file_name.endswith('.gz'):
            gz_file_path = os.path.join(mut_path, file_name)
            break
    list_mutations = []
    line_count = 0
    with gzip.open(gz_file_path , 'rt') as file:
        for line in file:
            line_count += 1
            if line_count > 7:  # Start extracting mutations from the sixth line
                columns = line.strip().split('\t')
                if not columns[35]: # Non-CDS mutations will be skipped
                    continue
                mutation = (columns[4], columns[5], columns[9], columns[10], columns[12])  # Assuming the mutation column is the first column
                if columns[9] == "SNP":
                    new_mutation= f"{mutation[0]}:g.{mutation[1]}{mutation[3]}>{mutation[4]}"
                elif columns[9] == "DEL":
                    new_mutation= f"{mutation[0]}:g.{mutation[1]}del{mutation[3]}"
                elif columns[9] == "INS":
                    new_mutation= f"{mutation[0]}:g.{mutation[1]}_{columns[6]}ins{mutation[4]}"
                else:
                    continue
                list_mutations.append((new_mutation, columns[0], columns[37]))
    return list_mutations, columns[15], line_count - 8

def find_patientID(file_id, cancer_type):
    df = pd.read_csv(f"/private7/projects/Netanel/dowmoloding_MHC/HLA-types/{cancer_type}.txt", delimiter='\t', header=None)
    all_files_id = df.iloc[:, 1]
    index = all_files_id[all_files_id == file_id].index
    row_index = index[0]
    return df.iloc[row_index, 0]

def find_HLA(patient_id, cancer_type):
    df = pd.read_csv(f"/private7/projects/Netanel/dowmoloding_MHC/HLA-types/{cancer_type}.tsv",  sep='\t')
    patient_row = df[df.iloc[:, 0] == patient_id]
    HLA_list = []
    if not patient_row.empty:
        hlas = patient_row.iloc[:, 2:].values.flatten().tolist()
        for hla in hlas:
            parts = hla.split('*')
            if len(parts) == 2:
                converted_hla = f'{parts[0]}{parts[1][:2]}{parts[1][3:]}'
                HLA_list.append(converted_hla)
    return HLA_list

def main():
    for primary_site in os.listdir("/home/alu/netlandes/MHCpan/projects2"):
        # First we will check if there is a folder for this type of cancer:
        fp = f"/home/alu/netlandes/MHCpan/final_results/{primary_site}"
        if not os.path.isdir(fp):
            os.mkdir(fp)
        else: # If there is already a folder from a previous run - we will delete the previous results
            dir_list = os.listdir(fp)
            for dir_name in dir_list:
                dir_path = os.path.join(f"/home/alu/netlandes/MHCpan/final_results/{primary_site}/", dir_name)
                os.system('rm -rf {}'.format(dir_path))
        for file_id in os.listdir(f"/home/alu/netlandes/MHCpan/projects2/{primary_site}"):
            list_of_mutation, p_id, num_mutation = creats_list_of_mutation(f"/home/alu/netlandes/MHCpan/projects2/{primary_site}/{file_id}")
            t = p_id.split('-')[0:3]
            patient_ID = '-'.join(t)
            list_HLA = find_HLA(patient_ID, primary_site)
            list_of_HLA = list(set(list_HLA))
            
            # final_results = []
            
            # If no HLA is found for the current patient, it is skipped:
            if not list_of_HLA:
                continue
            
            with open (f"/home/alu/netlandes/MHCpan/final_results/{primary_site}/{patient_ID}.tsv", 'w') as final_file:
                final_file.write("#NEO-ADAR-TIGEN 1.0\n")
                final_file.write(f"#Patient ID: {patient_ID}\n")
                final_file.write(f"#File ID: {file_id}\n")
                final_file.write(f"#Num of mutation: {num_mutation}\n")
                final_file.write(f"#Num of mutation in CDS: {len(list_of_mutation)}\n")
                final_file.write(f"Gene_Name:\tTranscript_ID:\tMutation:\tStrand:\tBest-target:\tGuide-RNA:\tRank:\tHLA:\tEDITS:\n")
                for mutation in list_of_mutation:
                    for edit_option in edits_list:
                        if edit_option == 0:
                            without_editing = True
                            single_editing = False
                            double_editing = False
                        elif edit_option == 1:
                            without_editing = False
                            single_editing = True
                            double_editing = False
                        else:
                            without_editing = False
                            single_editing = False
                            double_editing = True
                            
                        mutation_input = mutation[0]
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
                        subprocess.run(f'sh "find_position.sh" {mut_chr} {mut_pos}', shell=True)
                        output_file = "/home/alu/netlandes/MHCpan/intersects_output.bed"
                        if os.stat(output_file).st_size == 0:
                            print("The mutation is not in CDS\n")
                            continue
                        
                        isoforms_files = pd.read_csv(output_file, sep="\t", header=None)
                
                        # We will take the patient's specific transcript:
                        for index, row in isoforms_files.iterrows():
                            if mutation[2] in row[3]:
                                break

                        temp = seq[::-1]
                        # list_of_codons = []
                        frame_init_pos = int(mut_pos) - int((size_of_seq/2) - 1)
                        strand = row[4]
                        reading_frame = row[5]
                        if strand == "+":
                            cds_start_pos = row[1]
                            cds_end_pos = row[2]
                            count = 0
                            count2 = 1
                        else:
                            count = 1
                            count2 = 0
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
                            if (reading_frame == 0):
                                seq_frame = abs(frame_init_pos - cds_start_pos + count) % 3 
                            elif(reading_frame == 1):
                                seq_frame = abs(frame_init_pos - 1 - cds_start_pos) % 3 
                            else:
                                seq_frame = abs(frame_init_pos - 2 - cds_start_pos - count) % 3

                            if seq_frame == 0:
                                start = seq_frame 
                            elif seq_frame == 1:
                                start = seq_frame + 1
                            elif seq_frame == 2:
                                start = seq_frame - 1
                                
                        codon = [seq[i:i+3] for i in range(start, end, 3)]
                        # # Remove duplicates:
                        # if codon not in list_of_codons:
                        #     list_of_codons.append(codon)

                        # Creats mutatation in the seq for each isoform:
                        mut_isoforms = []
                        pos = int(size_of_seq/2) - start - count2 + count
                        if strand == "-":
                            if mut_type == "ins":
                                mut_seq = mut_seq[::-1]
                            pos = pos - len(mut_seq)

                        # for codon in list_of_codons:
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

                        # print(dna_seq)
                        # print(mut_dna_seq)
                        # print(Seq(mut_dna_seq).translate())

                        m = len(mut_seq)
                        if mut_type == "del":
                            m = 0

                        sb_tuple = []
                        wb_tuple = []

                        if without_editing:
                            for patient_HLA in list_of_HLA:
                                # Running the tool on the edited sequences:
                                run_MHCpan(mut_isoforms, patient_HLA, 1)
                                create_results_file(sb_tuple, wb_tuple, patient_HLA) 
                        
                            delete_files('/home/alu/netlandes/MHCpan/results1')
                            sequence_list = list(mut_isoforms[0]) 
                            gRNA = "".join(sequence_list[max(pos - num_of_nuc_around_mut, 0) : min(pos + num_of_nuc_around_mut + m, len(sequence_list))])
                        
                        if single_editing or double_editing:
                            
                            # Creating the edit in the mutant sequence:
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
                            
                            # Remove duplicates from the adar_mut list:
                            unique_ADAR_mut_list = list(set(ADAR_mut_seq_list))

                            if unique_ADAR_mut_list:
                                for patient_HLA in list_of_HLA:
                                    # Running the tool on the edited sequences:
                                    run_MHCpan([t[0] for t in unique_ADAR_mut_list], patient_HLA, 1)
                                    # if check_flag:
                                    #     # Running the tool on the non-edited sequences for comparison:
                                    #     run_MHCpan([t[2] for t in unique_ADAR_mut_list], patient_HLA, 2)
                                    create_results_file(sb_tuple, wb_tuple, patient_HLA) 
                            
                                delete_files('/home/alu/netlandes/MHCpan/results1')

                            # if check_flag:
                            #     delete_files('/home/alu/netlandes/MHCpan/results2')                    

                        uniq_sb = list(set(sb_tuple))
                        uniq_wb = list(set(wb_tuple))

                        if not uniq_sb and not uniq_wb:
                            sequence_list = list(mut_isoforms[0]) 
                            gRNA_NA = "".join(sequence_list[max(pos - num_of_nuc_around_mut, 0) : min(pos + num_of_nuc_around_mut + m, len(sequence_list))])
                            final_file.write(f"{mutation[1]}\t{mutation[2]}\t{mutation[0]}\t{strand}\t{'NA'}\t{gRNA_NA}\t{'NA'}\t{'NA'}\t{edit_option}\n")
                            continue

                        # Finding the strongest neoantigen that can be produced:
                        current_best_aso = sorted(uniq_sb + uniq_wb, key=lambda x: x[2])[0]
                        if single_editing or double_editing:
                            sb_seq = Seq(unique_ADAR_mut_list[current_best_aso[0]][0])
                        else:
                            sb_seq = Seq(mut_isoforms[current_best_aso[0]])

                        protein = sb_seq.translate()
                        if "*" in str(protein):
                            protein = str(protein)
                            protein = protein.replace("*", "X")
                            protein = Seq(protein)
                        pos = protein.find(current_best_aso[1])
                        ASO_target = sb_seq[(pos * 3):((pos + len(current_best_aso[1])) * 3)]
                        if single_editing or double_editing:
                            best_ASO = (unique_ADAR_mut_list[current_best_aso[0]][1], ASO_target, current_best_aso[2], mutation, current_best_aso[3], strand)
                        else:
                            best_ASO = (gRNA, ASO_target, current_best_aso[2], mutation, current_best_aso[3], strand)
                        final_file.write(f"{best_ASO[3][1]}\t{best_ASO[3][2]}\t{best_ASO[3][0]}\t{best_ASO[5]}\t{str(best_ASO[1])}\t{str(best_ASO[0])}\t{best_ASO[2]}\t{best_ASO[4]}\t{edit_option}\n")

                        # final_results.append((unique_ADAR_mut_list[current_best_aso[0]][1], ASO_target, current_best_aso[2], mutation, current_best_aso[3], strand))
                
                # sorted_final_results = sorted(final_results, key=lambda x: x[2])
                
                # # For each patient we will include in the report the three mutations with the lowest rank:
                # for best_ASO in sorted_final_results[0:3]:    
                #     final_file.write(f"{best_ASO[3][1]}\t{best_ASO[3][2]}\t{best_ASO[3][0]}\t{best_ASO[5]}\t{str(best_ASO[1])}\t{str(best_ASO[0])}\t{best_ASO[2]}\t{best_ASO[4]}\n")
                
if __name__ == "__main__":
    main()


