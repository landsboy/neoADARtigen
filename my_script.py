import subprocess
import pandas as pd
import re
import os
from Bio.Seq import Seq

# Global variables:
mutation_size = 6
size_of_seq = 52
mut_num = 50 #max 100
mer_lenght = 9

def most_common_mut(mut_num, primary_site):
    f_path = f"/home/alu/netlandes/MHCpan/frequent-mutations/{primary_site}.tsv"
    # Read the TSV file into a pandas DataFrame
    df = pd.read_csv(f_path, delimiter='\t')

    # Filter rows where the third column does not start with "Stop Gained"
    filtered_df = df[~df.iloc[:, 2].str.startswith("Stop Gained")]

    # Extract the first n rows
    if mut_num > 100:
        mut_num = 100

    selected_rows = filtered_df.iloc[:mut_num, 0].values

    # Create a list of rows
    return selected_rows.tolist()

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
    # We want deletion and insertion mutations only:
    if 'del' in mutation_input:
        pattern = r'(\w+):g\.(\d+)([a-z]+)([A-Z]+)'
    elif 'ins' in mutation_input:
        pattern = r'(\w+):g\.(\d+)_\d+([a-z]+)([A-Z]+)'
    elif '>' in mutation_input:
        pattern = r'(\w+):g\.(\d+)(\w)>(\w)'
    else:
        return None
    return re.search(pattern, mutation_input)

def checking_HLA(input_HLA):
    with open('/home/alu/netlandes/MHCpan/HLA-list.txt', 'r') as file:
        lines = [line.strip() for line in file.readlines()]
    if input_HLA not in lines:
        print("Look at HLA-list.txt and insert HLA according to the required format")
        exit(1)

def run_MHCpan(seq_list, patient_HLA, strand):
    count = 0
    a = ">Temp:"
    for seq in seq_list:
        count += 1
        if strand == "+":
            mut_aa_seq = Seq(seq).translate()
        else:
            mut_aa_seq = Seq(seq).complement().translate()
        with open("/home/alu/netlandes/MHCpan/input_seq.fsa", 'w') as f:
            f.write(a + '\n')
            f.write(str(mut_aa_seq))
            output_file = f"/home/alu/netlandes/MHCpan/results/result{count}.xls"  
        subprocess.run(f'"/private/common/Software/netMHC-4.0/Linux_x86_64/bin/netMHC" -f "/home/alu/netlandes/MHCpan/input_seq.fsa" -tdir "/tmp" -version "/private/common/Software/netMHC-4.0/Linux_x86_64/data/version " -hlalist "/private/common/Software/netMHC-4.0/data/allelelist" -syn "/private/common/Software/netMHC-4.0/Linux_x86_64/data/synlists/%s.synlist" -thrfmt "/private/common/Software/netMHC-4.0/threshold/%s.thr" -rdir "/private/common/Software/netMHC-4.0/Linux_x86_64" -xlsfile {output_file} -xls -a {patient_HLA}', shell=True)

def pre_proccesing(f):
    with open(f"/home/alu/netlandes/MHCpan/results/{f}", 'r') as file:
        first_line = file.readline()
        # Open a new file in write mode
        with open(f"/home/alu/netlandes/MHCpan/new_results/{f}", 'w') as f_new:
            # Iterate over the remaining lines and write them to the new file
            for line in file:
                f_new.write(line)

def create_results_file():
    dir_files = os.listdir("/home/alu/netlandes/MHCpan/results/")
    i = 0
    tuples_of_SB = []
    tuples_of_WB = []
    for file in dir_files:
        pre_proccesing(file)
        path = "/home/alu/netlandes/MHCpan/new_results/" + file
        results = pd.read_csv(path, sep='\t')
        new_results = results.loc[results['N_binders'] == 1]
        if not new_results.empty:
            for _, row in new_results.iterrows():
                if(row['Rank'] <= 0.5):
                    tuples_of_SB.append((i , row['Peptide'], row['Rank']))
                else:
                    tuples_of_WB.append((i , row['Peptide'], row['Rank']))
        i += 1
    return tuples_of_SB, tuples_of_WB


if __name__ == "__main__":
    best_ASO = ()
    primary_site = "brain"
    list_of_mutation = ["chr3:g.179199150delGAA"]
    # list_of_mutation = most_common_mut(mut_num, primary_site)
    list_of_HLA = ["HLA-A0101"]
    # checking_HLA(patient_HLA)
    first = True

    for mutation_input in list_of_mutation:
        match = checking_input(mutation_input)
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

        # If there is a deletion or insertion mutation, we dont want to have FS:
        if (mut_type != "Substitution" and len(mut_seq) % 3 != 0):
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
            if flag:
                if (reading_frame == 0):
                    seq_frame = abs(frame_init_pos - cds_start_pos + count) % 3 
                elif(reading_frame == 1):
                    seq_frame = abs(frame_init_pos - 1 - cds_start_pos) % 3 
                else:
                    seq_frame = abs(frame_init_pos - 2 - cds_start_pos - count) % 3
                if seq_frame == 0:
                    start = seq_frame 
                if seq_frame == 1:
                    start = seq_frame + 1
                if seq_frame == 2:
                    start = seq_frame - 1
            codons = [seq[i:i+3] for i in range(start, end, 3)]
            # print(seq_frame)
            # print(frame_init_pos - cds_start_pos + count)
            # Remove duplicates:
            if codons not in list_of_codons:
                list_of_codons.append(codons)
        print(list_of_codons)

        # Creats mutatation in the seq for each isoform:
        mut_isoforms = []
        pos1 = abs(int(mut_pos) - frame_init_pos) + 1
        pos = int(size_of_seq/2) - start - count2 + count
        if strand == "-":
            pos = pos -len(mut_seq)
        print (pos1)
        print (pos)
        for codon in list_of_codons:
            # Convert the list of codons to a DNA sequence
            dna_seq = ''.join(codon)
            if mut_type == "del":
                # Remove the sequence from the DNA sequence:
                mut_dna_seq = dna_seq[:pos] + dna_seq[pos + len(mut_seq):]
            elif mut_type == "ins":
                # Insert the sequence to the DNA sequence:
                mut_dna_seq = dna_seq[:pos1] + mut_seq + dna_seq[pos1:]
            else:
                mut_dna_seq = dna_seq[:pos] + mut_seq + dna_seq[pos + 1:]
            mut_isoforms.append(mut_dna_seq)
        if strand == "+":
            print(Seq(dna_seq))
            print(Seq(mut_isoforms[0]))
        else:
            print(Seq(dna_seq))
            print(Seq(mut_isoforms[0]))


        ADAR_mut_seq_list = []
        for isoform in mut_isoforms:
            if len(isoform) < mer_lenght*3:
                continue
            sequence_list = list(isoform)
            # Iterate over the sequence and replace "A" with "G" when necessary
            for i in range(1, len(sequence_list)):
                sequence_list = list(isoform)
                if sequence_list[i] == "A" and sequence_list[i - 1] != "G":
                    new_string = sequence_list
                    new_string[i] = "G"
                    ADAR_mut_seq_list.append("".join(new_string))

        # Remove duplicates from the adar_mut list:
        unique_ADAR_mut_list = list(set(ADAR_mut_seq_list))
        print(unique_ADAR_mut_list)

        if unique_ADAR_mut_list:
            for patient_HLA in list_of_HLA:
                run_MHCpan(unique_ADAR_mut_list, patient_HLA, strand)
                SB, WB = create_results_file() 
        
            delete_files('/home/alu/netlandes/MHCpan/results')
            delete_files('/home/alu/netlandes/MHCpan/new_results')
            # Get a list of all files in the folder
            if first:
                file_list = os.listdir("/home/alu/netlandes/MHCpan/final_results")
                # Iterate over the files and delete each one
                for file_name in file_list:
                    file_path = os.path.join("/home/alu/netlandes/MHCpan/final_results", file_name)
                    if os.path.isfile(file_path):
                        os.remove(file_path)
            
            uniq_sb = list(set(SB))
            uniq_wb = list(set(WB))
            # print(uniq_sb)
            # print(uniq_wb)

            if not uniq_sb and not uniq_wb:
                continue
            current_best_aso = ()
            current_first = True
            if uniq_sb:
                for sb in uniq_sb:
                    sb_seq = Seq(unique_ADAR_mut_list[sb[0]])
                    if strand == "+":
                        protein = sb_seq.translate()
                    else:
                        protein = sb_seq.complement().translate()
                    if "*" in str(protein):
                        protein = str(protein)
                        protein = protein.replace("*", "X")
                        protein = Seq(protein)
                    pos = protein.find(sb[1])
                    ASO_target = sb_seq[(pos * 3):((pos + len(sb[1])) * 3)]
                    if current_first:
                        current_best_aso = (unique_ADAR_mut_list[sb[0]], ASO_target, sb[2], mutation_input)
                        current_first = False
                    else:
                        if sb[2] < current_best_aso[2]:
                            current_best_aso = (unique_ADAR_mut_list[sb[0]], ASO_target, sb[2], mutation_input)
                    if first:
                        first = False
                        best_ASO = (unique_ADAR_mut_list[sb[0]], ASO_target, sb[2], mutation_input)
                        continue
                    if sb[2] < best_ASO[2]:
                        best_ASO = (unique_ADAR_mut_list[sb[0]], ASO_target, sb[2], mutation_input)
            if uniq_wb:
                for wb in uniq_wb:
                    wb_seq = Seq(unique_ADAR_mut_list[wb[0]])
                    if strand == "+":
                        protein = wb_seq.translate()
                    else:
                        protein = wb_seq.complement().translate()
                    if "*" in str(protein):
                        protein = str(protein)
                        protein = protein.replace("*", "X")
                        protein = Seq(protein)
                    pos = protein.find(wb[1])
                    ASO_target = wb_seq[(pos * 3):((pos + len(wb[1])) * 3)]
                    if current_first:
                        current_best_aso = (unique_ADAR_mut_list[wb[0]], ASO_target, wb[2], mutation_input)
                        current_first = False
                    else:
                        if wb[2] < current_best_aso[2]:
                            current_best_aso = (unique_ADAR_mut_list[wb[0]], ASO_target, wb[2], mutation_input)
                    if first:
                        first = False
                        best_ASO = (unique_ADAR_mut_list[wb[0]], ASO_target, wb[2] ,mutation_input)
                        continue
                    if wb[2] < best_ASO[2]:
                        best_ASO = (unique_ADAR_mut_list[wb[0]], ASO_target, wb[2], mutation_input)
            if current_best_aso:
                with open(f"/home/alu/netlandes/MHCpan/final_results/{mutation_input}", 'w') as f:
                    f.write(f"The best ASO to this mut is: {str(current_best_aso[1])} with rank of {current_best_aso[2]}")
    if best_ASO:
        print(f"\nThe best ASO is: {str(best_ASO[1])} with {best_ASO[2]} rank in the {best_ASO[3]} mutation\n")
    else:
        print("\nNo results!!")