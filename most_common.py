import subprocess
import pandas as pd
import re
import os
from Bio.Seq import Seq
from itertools import combinations, groupby


# Global variables:
k=10

num_of_nuc_around_mut = 20
mutation_size = 6
size_of_seq = 2*num_of_nuc_around_mut + 2*26 + 1
mer_lenght = 9

mut_num = 100 #max 250. 

list_of_weizmann_HLA = ['HLA-A0101', 'HLA-A0201', 'HLA-A0301', 'HLA-A1101', 'HLA-A2402', 'HLA-A2501', 'HLA-A2902 ', 'HLA-A3001', 'HLA-A3101', 'HLA-A3201 ', 'HLA-A3301 ', 'HLA-A6801', 'HLA-B0702', 'HLA-B0801', 'HLA-B1501', 'HLA-B1801', 'HLA-B2705', 'HLA-B3501', 'HLA-B4001', 'HLA-B4002', 'HLA-B4032', 'HLA-B4402', 'HLA-B4403', 'HLA-B5101', 'HLA-B5701', 'HLA-C0102', 'HLA-C0205', 'HLA-C0303', 'HLA-C0304', 'HLA-C0401', 'HLA-C0501', 'HLA-C0602', 'HLA-C0701', 'HLA-C0702', 'HLA-C0802', 'HLA-C1203', 'HLA-C1601']
# list_of_100_HLA = ['HLA-A0201', 'HLA-C0701', 'HLA-A0101', 'HLA-C0702', 'HLA-C0401', 'HLA-A0301', 'HLA-B0702', 'HLA-A2402', 'HLA-B0801', 'HLA-C0602', 'HLA-B4402', 'HLA-C0501', 'HLA-C0304', 'HLA-A1101', 'HLA-B3501', 'HLA-C1203', 'HLA-B5101', 'HLA-B1501', 'HLA-B4403', 'HLA-B1801', 'HLA-C0303', 'HLA-B4001', 'HLA-C0102', 'HLA-C0202', 'HLA-C1601', 'HLA-C0802', 'HLA-A2601', 'HLA-A3201', 'HLA-A6801', 'HLA-B2705', 'HLA-A2301', 'HLA-A2902', 'HLA-B5701', 'HLA-A3101', 'HLA-B1402', 'HLA-B1302', 'HLA-C1502', 'HLA-B3801', 'HLA-A2501', 'HLA-A3001', 'HLA-B4901', 'HLA-B3503', 'HLA-A3303', 'HLA-C1402', 'HLA-C1701', 'HLA-B5201', 'HLA-C0704', 'HLA-B5501', 'HLA-B4002', 'HLA-A3002', 'HLA-C1202', 'HLA-B5801', 'HLA-B5301', 'HLA-B3901', 'HLA-A6802', 'HLA-B3701', 'HLA-B3502', 'HLA-B5001', 'HLA-A3301', 'HLA-B4501', 'HLA-A0205', 'HLA-B1401', 'HLA-C0801', 'HLA-C0302', 'HLA-C0210', 'HLA-B1503', 'HLA-C1505', 'HLA-B5601', 'HLA-A7401', 'HLA-B0705', 'HLA-A2901', 'HLA-A0206', 'HLA-B4601', 'HLA-B4101', 'HLA-B4201', 'HLA-B1502', 'HLA-B2702', 'HLA-B3906', 'HLA-B4102', 'HLA-B5802', 'HLA-A6601', 'HLA-A0202', 'HLA-B3508', 'HLA-C1602', 'HLA-B1510', 'HLA-A0207', 'HLA-B5703', 'HLA-B1517', 'HLA-A3402', 'HLA-A0203', 'HLA-A3601', 'HLA-C1801', 'HLA-B4405', 'HLA-A3004', 'HLA-B3802', 'HLA-B1516', 'HLA-B1301', 'HLA-A0302', 'HLA-A2403', 'HLA-A0325']
edits_list = [0,1,2]

def setting_Environment_Variables():
    # Define the new values for the environment variables
    os.environ["PLATFORM"] = "Linux_x86_64"
    os.environ["NMHOME"] = "/private7/projects/Netanel/netMHCpan-4.1"
    os.environ["NETMHCpan"] = os.environ["NMHOME"] + "/" + os.environ["PLATFORM"]
    os.environ["NetMHCpanWWWPATH"] = "/services/NetMHCpan/tmp/"
    os.environ["NetMHCpanWWWDIR"] = "/usr/opt/www/pub/CBS/services/NetMHCpan/tmp"
    os.environ["TMPDIR"] = "/tmp"
    os.environ["DTUIBSWWW"] = "www"

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

def creat_FASTA_file(list_of_seq, k):
    with open(f"/home/alu/netlandes/MHCpan/temp_folder/input_seq{k}.fsa", "w") as fasta_file:
        for i, sequence in enumerate(list_of_seq):
            fasta_file.write(f">{i}\n{str(Seq(sequence).translate())}\n")

def run_MHCpan(patient_HLA, k):
    output = subprocess.run(f'"/private7/projects/Netanel/netMHCpan-4.1/Linux_x86_64/bin/netMHCpan" -f "/home/alu/netlandes/MHCpan/temp_folder/input_seq{k}.fsa" -a {patient_HLA} -l {mer_lenght} -BA -t 2.5', shell=True, capture_output=True, text=True)
    return output.stdout

def create_results_file(tuples_of_SB, tuples_of_WB, output):
    # If we didn't get any result:
    if not output:
        return
    output_list = output.split("\n")
    for result in output_list:   
        if result[-2:] == "WB" and result[0] != "#":
            final_result = [value for value in result.split(" ") if value.strip()]
            tuples_of_WB.append((int(final_result[10]), final_result[2], final_result[12], final_result[1].replace('*', ''), final_result[15]))
        if result[-2:] == "SB" and result[0] != "#":
            final_result = [value for value in result.split(" ") if value.strip()]
            tuples_of_SB.append((int(final_result[10]), final_result[2], final_result[12], final_result[1].replace('*', ''), final_result[15])) 

def most_common_mut(mut_num, primary_site):
    mut_list = []
    f_path = f"/private7/projects/Netanel/MHC_COMMON_project/all_cancers/{primary_site}.tsv"  #frequent-mutations

    with open (f_path, 'r') as mut_file:
        next(mut_file)
        if mut_num > 250:
            mut_num = 250
        for line in mut_file:
            if len(mut_list) == mut_num:
                break
            columns = line.strip().split('\t')
            if columns[2].startswith("Intron") or columns[2].startswith("Stop Gained"):
                continue
            mut_val = columns[2].split(" ")
            if "*" in mut_val[-1]: # Mutations that cause FS within our editing framework - will be removed
                if "fs" in mut_val[-1]:
                    num_AA = mut_val[-1].split("*")[1]
                    if  num_AA == "?" or int(num_AA) < 7:
                        continue
                else:
                    continue

            mut_list.append([columns[0], mut_val[-2], mut_val[-1], columns[3].split(",")[2]])
            # mut_list.append([columns[0], mut_val[1], mut_val[2], columns[3].split(",")[2]])

    return mut_list
            
def preprocessing(mutation):
    mutation_input = mutation[0]
    match = checking_input(mutation_input)
    pos_list = []
    if match == None: # A mutation that is not a deletion insertion Or Substitution
        return
    mut_chr = match.group(1)
    mut_pos = match.group(2)

    if '>' in mutation_input:
        mut_type = "Substitution"
    else:
        mut_type = match.group(3)
    mut_seq = match.group(4)

    # Run a program that accepts as input a chromosome and location and returns the bases around it:
    p = subprocess.run(f'sh "/home/alu/netlandes/MHCpan/find_seq_of_mutation.sh" {mut_chr} {mut_pos} {int(size_of_seq / 2)} {k}', capture_output=True,text=True, shell=True)
    seq = p.stdout.replace("\n", "") 

    # We will run a program that will return the reading frame of the sequence:
    q = subprocess.run(f'sh "/home/alu/netlandes/MHCpan/find_position.sh" {mut_chr} {mut_pos} {k}', shell=True, capture_output=True,text=True)
    res_list = [res for res in q.stdout.split("\n") if res]

    if not res_list:
        return
                
    # We will take the patient's specific transcript:
    for res in res_list:
        row = [value for value in res.split("\t") if value.strip()]
        if mutation[2] in row[3]:
            break

    temp = seq[::-1]
    frame_init_pos = int(mut_pos) - int((size_of_seq/2) - 1)
    strand = row[4]
    reading_frame = int(row[5])
    if strand == "+":
        cds_start_pos = int(row[1])
        cds_end_pos = int(row[2])
        count = 0
        count2 = 1
    else:
        count = 1
        count2 = 0
        seq = temp
        frame_init_pos = int(mut_pos) + int((size_of_seq/2) - 1)
        cds_start_pos = int(row[2])
        cds_end_pos = int(row[1])
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

    # Creats mutatation in the seq for each isoform:
    mut_isoforms = []
    pos = int(size_of_seq/2) - start - count2 + count
    if strand == "-":
        if mut_type == "ins":
            mut_seq = mut_seq[::-1]
        pos = pos - len(mut_seq)

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
        return 

    # If the mutation is on the negative strand, we will look for editing on the negative strand:
    if strand == "-":
        mut_dna_seq = str(Seq(mut_dna_seq).complement())

    mut_isoforms.append(mut_dna_seq)

    m = len(mut_seq)
    if mut_type == "del":
        m = 0

    return mut_isoforms, strand, pos_list, m, pos
        
def edit_option_initialization(edit_option):
    if edit_option == 0:
        return True, False, False
    elif edit_option == 1:
        return False, True, False
    else:
        return False, False, True
    

def main():
    setting_Environment_Variables()
    with open (f"/home/alu/netlandes/MHCpan/most_common_results_weizmann.tsv", 'w') as final_file:
        final_file.write(f"Gene-Name\tAA-change\tMutation\tStrand\tPrevalence-percentage\tTarget\tGuide-RNA\tRank\tAffinity(nM)\tHLA\tEdits\n")
        file_names  = os.listdir("/private7/projects/Netanel/MHC_COMMON_project/all_cancers")  
        for primary_site in [os.path.splitext(file_name)[0] for file_name in file_names]:
            list_of_mutation = most_common_mut(mut_num, primary_site)
            for mut in list_of_mutation:
                fun_result = preprocessing(mut)
                if fun_result is not None:
                    mut_isoforms, strand, pos_list, m, pos = fun_result
                else:
                    continue

                for edit_option in edits_list:
                    without_editing, single_editing, double_editing = edit_option_initialization(edit_option)

                    sb_tuple = []
                    wb_tuple = []

                    if without_editing:
                        creat_FASTA_file(mut_isoforms, k)
                        for hla in list_of_weizmann_HLA:
                            output = run_MHCpan(hla, k)
                            create_results_file(sb_tuple, wb_tuple, output) 
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
                            creat_FASTA_file([t[0] for t in unique_ADAR_mut_list], k)
                            # Running the tool on the edited sequences:
                            for hla in list_of_weizmann_HLA:
                                output = run_MHCpan(hla, k)
                                create_results_file(sb_tuple, wb_tuple, output)
                
                    uniq_sb = list(set(sb_tuple))
                    uniq_wb = list(set(wb_tuple))


                    # For each HLA we will find the result with the lowest rank:
                    sorted_data = sorted(uniq_sb + uniq_wb, key=lambda x: x[3])
                    grouped_data = {key: list(group) for key, group in groupby(sorted_data, key=lambda x: x[3])}
                    sb_and_wb = [min(group, key=lambda x: x[2]) for group in grouped_data.values()]

                    if not sb_and_wb:
                        for hla in sorted(list_of_weizmann_HLA):
                            final_file.write(f"{mut[1]}\t{mut[2]}\t{mut[0]}\t{strand}\t{mut[3]}\tNA\tNA\tNA\tNA\t{hla}\t{edit_option}\n")
                        continue    

                    result_dict = {tup[3]: tup for tup in sb_and_wb}

                    for hla in sorted(list_of_weizmann_HLA):
                        if hla in result_dict:
                            sb_wb = result_dict[hla]
                            if single_editing or double_editing:
                                sb_seq = Seq(unique_ADAR_mut_list[sb_wb[0]][0])
                                guide_RNA = unique_ADAR_mut_list[sb_wb[0]][1]
                            else:
                                sb_seq = Seq(mut_isoforms[sb_wb[0]])
                                guide_RNA = gRNA
                            protein = sb_seq.translate()
                            if "*" in str(protein):
                                protein = str(protein)
                                protein = protein.replace("*", "X")
                                protein = Seq(protein)
                            pos_p = protein.find(sb_wb[1])
                            ASO_target = sb_seq[(pos_p * 3):((pos_p + len(sb_wb[1])) * 3)]
                            final_file.write(f"{mut[1]}\t{mut[2]}\t{mut[0]}\t{strand}\t{mut[3]}\t{str(ASO_target)}\t{str(guide_RNA)}\t{sb_wb[2]}\t{sb_wb[4]}\t{sb_wb[3]}\t{edit_option}\n")
                        else:
                            final_file.write(f"{mut[1]}\t{mut[2]}\t{mut[0]}\t{strand}\t{mut[3]}\tNA\tNA\tNA\tNA\t{hla}\t{edit_option}\n")
          
if __name__ == "__main__":
    main()
