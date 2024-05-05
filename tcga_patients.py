import argparse
import subprocess
import pandas as pd
import re
import os
from Bio.Seq import Seq
import gzip
from itertools import combinations
import csv

# Global variables:
k=5

num_of_nuc_around_mut = 20
mutation_size = 6
size_of_seq = 2*num_of_nuc_around_mut + 2*26 + 1
mer_lenght = 9
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

def delete_previouse_results(primary_site):
    fp = f"/home/alu/netlandes/MHCpan/final_results/{primary_site}"
    dir_list = os.listdir(fp)
    for dir_name in dir_list:
        dir_path = os.path.join(f"/home/alu/netlandes/MHCpan/final_results/{primary_site}/", dir_name)
        os.system('rm -rf {}'.format(dir_path))

def find_duplicate_patients_in_data(primary_site):
    all_patient = []
    for file_id in os.listdir(f"/private7/projects/Netanel/MHC_TCGA_projects/projects{k}/{primary_site}"):
        mut_path = f"/private7/projects/Netanel/MHC_TCGA_projects/projects{k}/{primary_site}/{file_id}"
        gz_file_path = None
        for file_name in os.listdir(mut_path):
            if file_name.endswith('.gz'):
                gz_file_path = os.path.join(mut_path, file_name)
                break
        line_count = 0
        with gzip.open(gz_file_path , 'rt') as file:
            for line in file:
                line_count += 1
                if line_count > 8:  # Start extracting mutations from the 7th line
                    columns = line.strip().split('\t')
                    temp_val = columns[15].split('-')[0:3]
                    all_patient.append(('-').join(temp_val))
                    break
    # Create an empty dictionary to store the count of occurrences of each patient name
    name_count = {}
    # Count the occurrences of each patient name in the list
    for name in all_patient:
        if name in name_count:
            name_count[name] += 1
        else:
            name_count[name] = 1
    return [name for name, count in name_count.items() if count > 1]

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
            if line_count > 8:  # Start extracting mutations from the 7th line
                columns = line.strip().split('\t')
                if not columns[35]: # Non-CDS mutations will be skipped
                    continue
                if columns[4] == "chrM":
                    continue
                HGVSp_Short = columns[36].split('.')[-1]
                if "*" in HGVSp_Short: # Mutations that cause FS within our editing framework - will be removed
                    if "fs" in HGVSp_Short:
                        num_AA = HGVSp_Short.split("*")[-1]
                        if num_AA == "?" or int(num_AA) < 7:
                            continue
                    else:
                        continue
                mutation = (columns[4], columns[5], columns[9], columns[10], columns[12])  
                if columns[9] == "SNP":
                    new_mutation= f"{mutation[0]}:g.{mutation[1]}{mutation[3]}>{mutation[4]}"
                elif columns[9] == "DEL":
                    new_mutation= f"{mutation[0]}:g.{mutation[1]}del{mutation[3]}"
                elif columns[9] == "INS":
                    new_mutation= f"{mutation[0]}:g.{mutation[1]}_{columns[6]}ins{mutation[4]}"
                else:
                    continue
                list_mutations.append((new_mutation, columns[0], columns[37], HGVSp_Short))
    return list_mutations, columns[15], line_count - 8

def find_HLA(patient_id, df):
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

    # # We want the size of the insertion to be no more than 6bp:
    # if not (len(mut_seq) > mutation_size):
    #     return

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
        
def adding_TPM(primary_site, duplicats_patients):
    with open(f'/private7/projects/Netanel/MHC_TCGA_projects/projects{k}/UUIDtoCASE.tsv', 'r') as tsvfile:
        reader = csv.reader(tsvfile, delimiter='\t')
        data = [(row[1], row[2], row[0]) for row in reader]

    if not os.path.isdir(f"/home/alu/netlandes/MHCpan/final_results/final_{primary_site}"):
            os.mkdir(f"/home/alu/netlandes/MHCpan/final_results/final_{primary_site}")


    for filename in os.listdir(f"/home/alu/netlandes/MHCpan/final_results/{primary_site}"):
        mapping = {}
        if filename.split(".")[0].count("-") > 3:
            parts = (filename.split(".")[0]).split("-")
            patient_id = "-".join(parts[:4])
            temp = "-".join(parts[:3])
        else:
            patient_id = filename.split(".")[0]
            temp = filename.split(".")[0]
        for item in data:
            if (patient_id not in duplicats_patients and item[0] == patient_id) or (temp in duplicats_patients and item[2] == patient_id):
                folder_path = f'/private7/projects/Netanel/MHC_TCGA_projects/projects{k}/TPM/{item[1]}/'
                if not os.path.exists(folder_path):
                    break
                files_in_folder = os.listdir(folder_path)
                tsv_files = [file for file in files_in_folder if file.endswith('.tsv')]
                with open(folder_path + tsv_files[0] , 'r') as TPM_file:
                    for i, line in enumerate(TPM_file):
                        if i >= 6:
                            columns = line.strip().split('\t')
                            key = columns[6]  
                            value = columns[1]  
                            mapping[value] = key

                with open(f"/home/alu/netlandes/MHCpan/final_results/{primary_site}/" + filename , 'r') as patient_file:
                    with open(f"/home/alu/netlandes/MHCpan/final_results/final_{primary_site}/" + filename, 'w') as new_patient_file:
                        for i, line in enumerate(patient_file):
                            if i >= 6:
                                gene = line.strip().split('\t')[0]
                                TPM_value = mapping.get(gene, '')
                                if TPM_value:
                                    new_patient_file.write(line.replace("\n", "") + "\t" + TPM_value + "\n")
                                else:
                                    new_patient_file.write(line) 
                            elif i == 5:
                                new_patient_file.write(line.replace("\n", "") + "\t" + "TPM" + "\n")
                            else:
                                new_patient_file.write(line) 
                break

def edit_option_initialization(edit_option):
    if edit_option == 0:
        return True, False, False
    elif edit_option == 1:
        return False, True, False
    else:
        return False, False, True

def main():
    setting_Environment_Variables()
    for primary_site in os.listdir(f"/private7/projects/Netanel/MHC_TCGA_projects/projects{k}"):
        if not os.path.isdir(os.path.join(f"/private7/projects/Netanel/MHC_TCGA_projects/projects{k}/", primary_site)) or primary_site == "TPM":
            continue
        if not os.path.isdir(f"/home/alu/netlandes/MHCpan/final_results/{primary_site}"):
            os.mkdir(f"/home/alu/netlandes/MHCpan/final_results/{primary_site}")
        # delete_previouse_results(primary_site)
        df = pd.read_csv(f"/private7/projects/Netanel/MHC_TCGA_HLA_types/{primary_site}.tsv",  sep='\t')
        duplicats_patients = find_duplicate_patients_in_data(primary_site)
        for file_id in os.listdir(f"/private7/projects/Netanel/MHC_TCGA_projects/projects{k}/{primary_site}"):
            list_of_mutation, p_id, num_mutation = creats_list_of_mutation(f"/private7/projects/Netanel/MHC_TCGA_projects/projects{k}/{primary_site}/{file_id}")
            t = p_id.split('-')[0:3]
            patient_ID = '-'.join(t)

            if patient_ID in duplicats_patients:
                res_path = f"/home/alu/netlandes/MHCpan/final_results/{primary_site}/{p_id}.tsv"
            else:
                res_path = f"/home/alu/netlandes/MHCpan/final_results/{primary_site}/{patient_ID}.tsv"

            list_HLA = find_HLA(patient_ID, df)
            list_of_HLA = list(set(list_HLA))

            # If no HLA is found for the current patient, it is skipped:
            if not list_of_HLA:
                continue

            with open (res_path, 'w') as final_file:
                final_file.write("#NEO-ADAR-TIGEN 1.0\n")
                final_file.write(f"#Patient ID: {patient_ID}\n")
                final_file.write(f"#File ID: {file_id}\n")
                final_file.write(f"#Num of mutation: {num_mutation}\n")
                final_file.write(f"#Num of mutation in CDS: {len(list_of_mutation)}\n")
                final_file.write("Gene_Name\tTranscript_ID\tMutation\tStrand\tHGVSp_Short\tBest-target\tGuide-RNA\tRank\tAffinity(nM)\tHLA\tEDITS\n")

                if len(list_of_mutation) >= 6000:
                    final_file.write("TODO\n")
                    continue

                for mutation in list_of_mutation:
                    fun_result = preprocessing(mutation)
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
                            for hla in list_of_HLA:
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
                                for hla in list_of_HLA:
                                    output = run_MHCpan(hla, k)
                                    create_results_file(sb_tuple, wb_tuple, output)
                 
                        uniq_sb = list(set(sb_tuple))
                        uniq_wb = list(set(wb_tuple))

                        if not uniq_sb and not uniq_wb:
                            sequence_list = list(mut_isoforms[0]) 
                            gRNA_NA = "".join(sequence_list[max(pos - num_of_nuc_around_mut, 0) : min(pos + num_of_nuc_around_mut + m, len(sequence_list))])
                            final_file.write(f"{mutation[1]}\t{mutation[2]}\t{mutation[0]}\t{strand}\t{mutation[3]}\t{'NA'}\t{gRNA_NA}\t{'NA'}\t{'NA'}\t{'NA'}\t{edit_option}\n")
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
                        pos_p = protein.find(current_best_aso[1])
                        ASO_target = sb_seq[(pos_p * 3):((pos_p + len(current_best_aso[1])) * 3)]
                        if single_editing or double_editing:
                            best_ASO = (unique_ADAR_mut_list[current_best_aso[0]][1], ASO_target, current_best_aso[2], mutation, current_best_aso[3], strand, current_best_aso[4])
                        else:
                            best_ASO = (gRNA, ASO_target, current_best_aso[2], mutation, current_best_aso[3], strand, current_best_aso[4])
                        final_file.write(f"{best_ASO[3][1]}\t{best_ASO[3][2]}\t{best_ASO[3][0]}\t{best_ASO[5]}\t{best_ASO[3][3]}\t{str(best_ASO[1])}\t{str(best_ASO[0])}\t{best_ASO[2]}\t{best_ASO[6]}\t{best_ASO[4]}\t{edit_option}\n")
        
        adding_TPM(primary_site, duplicats_patients)


if __name__ == "__main__":
    main()









































# NumEditingOptionsDict = {"All":[0,1,2], "NoEditing":[0], "SingleEditing":[1],"DoubleEditing":[2]}

    # # Create parser
    # class MyFormatter(argparse.ArgumentDefaultsHelpFormatter, argparse.RawTextHelpFormatter):
    #     pass

    # parser = argparse.ArgumentParser(formatter_class=MyFormatter,
    #                                  description='TODO')
    
    # parser.add_argument('-i', '--site_dir', dest='site_dir', action='store', metavar='root_dir',
    #                     nargs='?', default="/home/alu/netlandes/MHCpan/projects2", help='The input directory')
    # parser.add_argument('-o', '--output_dir', dest='output_dir', action='store', metavar='output_file', nargs='?',
    #                     default="/home/alu/netlandes/MHCpan/final_results", help='Output dir')
    # parser.add_argument('-l', '--log_path', dest='log_path', action='store', metavar='log_path', nargs='?',
    #                     default=None, help='Log file, default is to create in input dir')
    # parser.add_argument('--num_editing', dest='num_editing', choices=[key for key in NumEditingOptionsDict], nargs='+',
    #                     default=NumEditingOptionsDict['All'], help='TODO')

    # options = parser.parse_args()





            # with open (res_path, 'r+') as final_file:                
                # # Move the file pointer to the end of the file
                # final_file.seek(0, os.SEEK_END)
                
                # # Get the position of the last non-empty line
                # pos = final_file.tell()
                # while pos > 0:
                #     pos -= 1
                #     final_file.seek(pos, os.SEEK_SET)
                #     if final_file.read(1) == '\n':
                #         break

                # # Save the position to start appending new data
                # start_pos = pos + 1

                # # Write the new data
                # final_file.seek(start_pos, os.SEEK_SET)
