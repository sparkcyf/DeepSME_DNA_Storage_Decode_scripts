import pickle, os, tqdm
from Bio import SeqIO

def get_all_anchor(anchor):
    all_anchor = [anchor]
    # delete one
    subanchors = []
    for i in range(1,len(anchor)):
        subanchor = anchor[:i]+anchor[i+1:]
        subanchors.append(subanchor)
    all_anchor += list(set(subanchors))
    # insert one
    subanchors = []
    for i in range(1,len(anchor)):
        for insert in ['A','C','G','T']:
            subanchor = anchor[:i]+insert+anchor[i:]
            subanchors.append(subanchor)
    all_anchor += list(set(subanchors))
    return all_anchor

def search_anchor(dna, anchor, ref_position):
    # all_anchor = get_all_anchor(anchor)
    all_anchor = [anchor]
    
    true_position = [0,0,0]
    anchor_len = [0,0,0]
    
    ps_offset = 0 
    
    for ps in ref_position:
        ps += ps_offset
        for anchor in all_anchor:
            is_found = False
            search_path = [0] + [num if num == 0 else (-1) ** num * (num // 2) for num in range(2,(10+1)*2)] # 10+1 改 5+1
            for walk in search_path:
                if dna[ps+walk:ps+walk+len(anchor)] == anchor:
                    true_position[ref_position.index(ps-ps_offset)] = ps+walk
                    anchor_len[ref_position.index(ps-ps_offset)] = len(anchor)
                    is_found = True
                    ps_offset += walk
                    break
            if is_found:
                break
    
    return sum(1 for num in true_position if num != 0), abs(ps_offset), anchor_len

def reverse_complement(dna_sequence):
    complement_dict = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}
    reverse_sequence = dna_sequence[::-1]
    complement_sequence = ''.join([complement_dict[base] for base in reverse_sequence])
    return complement_sequence

if __name__ == "__main__":
    basecalled_dir = r'./'

    pickle_path = r'./readisc.P'

    bcd_path = r'./reference/barcodeljy_384.fasta'

    all_bcd_id = [(str(item.id)) for item in SeqIO.parse(bcd_path, 'fasta')]

    print('Loading pickle data...')
    with open(pickle_path, 'rb') as f:
        pickle_data = pickle.load(f)
    print('Pickle data loaded...')

    filelist = os.listdir(basecalled_dir)
    fastq_list = list(
                filter(lambda filename: os.path.split(filename)[-1].split('.')[-1] == 'fastq', filelist))

    import math, tqdm

    ref_num = 487

    st_flag = 'TTTCT'
    ed_flag = 'TAGAGC'

    grouped_seqs = [[] for _ in range(ref_num)]
    print('Reading all seqs...')
    fastq_folder = os.getenv('FASTQ_FOLDER')
    all_seqs = [(str(item.id), str(item.seq)) for item in SeqIO.parse(f'{fastq_folder}/basecalling.fastq', 'fastq')]
    print('All seqs read...')
    
    with tqdm.trange(len(all_seqs), desc='genome grouping...') as tbar:
        for seq in all_seqs:
            try:
                target_name = pickle_data[seq[0]]['T_n']
            except:
                tbar.update(1)
                continue
            query_start  = pickle_data[seq[0]]['Q_s']
            query_end  = pickle_data[seq[0]]['Q_e']

            # target_index = int(re.findall(r'\d+', target_name.split('_')[0])[0])

            if pickle_data[seq[0]]['R_s'] == '+':
                sequence = seq[1]
            if pickle_data[seq[0]]['R_s'] == '-':
                sequence = reverse_complement(seq[1])

            
            bcd_index = all_bcd_id.index(target_name)
            # print(target_name, bcd_index)
            st_index = sequence.find(st_flag)
            ed_index = sequence.find(ed_flag)

            if (st_index != -1) and (ed_index != -1) and 223<len(sequence[st_index:ed_index+len(ed_flag)])<263:
                
                temp_seq = (seq[0], sequence[st_index:ed_index+len(ed_flag)])

                dna = sequence[st_index:ed_index+len(ed_flag)]
            # if 223<len(sequence)<263:
                
                # temp_seq = (seq[0], sequence[st_index:ed_index+len(ed_flag)])

                # dna = sequence

                ref_position = [97, 142, 187]
                anchors = ['GGCCT', 'AAGGC', 'GCGCT']
                
                GGCCT_counts, GGCCT_offset, GGCCT_len = search_anchor(dna, 'GGCCT', ref_position)
                AAGGC_counts, AAGGC_offset, AAGGC_len = search_anchor(dna, 'AAGGC', ref_position)
                GCGCT_counts, GCGCT_offset, GCGCT_len = search_anchor(dna, 'GCGCT', ref_position)

                # GGCCT_counts, GGCCT_offset = search_anchor(dna, 'GGCCT', ref_position)
                # AAGGC_counts, AAGGC_offset = search_anchor(dna, 'AAGGC', ref_position)
                # GCGCT_counts, GCGCT_offset = search_anchor(dna, 'GCGCT', ref_position)

                
                anchor_counts = [GGCCT_counts, AAGGC_counts, GCGCT_counts]
                anchor_offsets = [GGCCT_offset, AAGGC_offset, GCGCT_offset]
                anchor_lens = [GGCCT_len, AAGGC_len, GCGCT_len]
                
                if all(item <= 1 for item in anchor_counts):
                # 宽松一点
                # if all(item < 1 for item in anchor_counts):
                    tbar.update(1)
                    continue

                max_counts_value = max(anchor_counts)
                max_counts_indices = [index for index, value in enumerate(anchor_counts) if value == max_counts_value]
                
                for i, data in enumerate(anchor_counts):
                    if data == 0:
                        anchor_offsets[i] = math.inf

                min_offsets_value = min(anchor_offsets)
                min_offsets_indices = [index for index, value in enumerate(anchor_offsets) if value == min_offsets_value]                    
                # print(min_offsets_value, min_offsets_indices)
                if len(max_counts_indices) == 0:
                    tbar.update(1)
                    continue
                elif len(max_counts_indices) == 1:
                    if max_counts_indices[0] == 0:
                        target_index = bcd_index
                    elif max_counts_indices[0] == 1:
                        target_index = bcd_index + 55
                    elif max_counts_indices[0] == 2:
                        target_index = bcd_index + 55 + 384
                elif len(max_counts_indices) > 1:
                    if min_offsets_indices[0] == 0:
                        target_index = bcd_index
                    elif min_offsets_indices[0] == 1:
                        target_index = bcd_index + 55
                    elif min_offsets_indices[0] == 2:
                        target_index = bcd_index + 55 + 384
                    else:
                        tbar.update(1)
                        continue

                # if all_bcd_id.index(target_name) > 7:
                #     target_index = all_bcd_id.index(target_name) + 7
                # else:
                #     GGCCT_counts, GGCCT_offset, GGCCT_len = search_anchor(dna, 'GGCCT', ref_position)
                #     AAGGC_counts, AAGGC_offset, AAGGC_len = search_anchor(dna, 'AAGGC', ref_position)


                #     # AAGGC_counts, AAGGC_offset = search_anchor(dna, 'AAGGC', ref_position)
                #     # GCGCT_counts, GCGCT_offset = search_anchor(dna, 'GCGCT', ref_position)
                    
                #     anchor_counts = [GGCCT_counts, AAGGC_counts]
                #     anchor_offsets = [GGCCT_offset, AAGGC_offset]
                #     anchor_lens = [GGCCT_len, AAGGC_len]
                    
                #     if all(item <= 1 for item in anchor_counts):
                #         continue 

                #     max_counts_value = max(anchor_counts)
                #     max_counts_indices = [index for index, value in enumerate(anchor_counts) if value == max_counts_value]
                    
                #     for i, data in enumerate(anchor_counts):
                #         if data == 0:
                #             anchor_offsets[i] = math.inf

                #     min_offsets_value = min(anchor_offsets)
                #     min_offsets_indices = [index for index, value in enumerate(anchor_offsets) if value == min_offsets_value]

                if target_index > 486:
                    tbar.update(1)
                    continue
                
                grouped_seqs[target_index].append((seq[0], dna))
                
            tbar.update(1)
    print('Seqs Grouping Done...')

    print('\n')

    print('Grouping result write working...')

    grouping_res_dir = r'./grouping_res/'

    if not os.path.exists(grouping_res_dir):
        os.makedirs(grouping_res_dir)

    for idx in tqdm.tqdm(range(ref_num)):
        with open(grouping_res_dir+f'grouping_res_{idx}.fasta', 'w') as f:
            for item in grouped_seqs[idx]:
                f.write(f'>{item[0]}\n')
                f.write(f'{item[1]}\n')

    print('Grouping result write Done...')