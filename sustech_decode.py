from base64 import decode
import numpy as np
from copy import deepcopy
from typing import List
from utils import composite_letter, Hypothesis, code_pattern, HypothesisTree
from scipy.stats import wasserstein_distance
from scipy.spatial import distance
from Bio import SeqIO
import tqdm
from reedsolo import RSCodec
import math
import time
import random
import itertools
import os
# import nupyck

DEBUG0 = False
DEBUG1 = False
INPUT_FASTA_PATH = os.environ['INPUT_FASTA_PATH']
"""
This version added GC content and max runs of homopolymer constrant.
Date : 
"""
def ran_hash(u):
    """
    哈希编码
    :param u: 编码值
    :return: 哈希值
    """
    v = u * 3935559000370003845 + 2691343689449507681
    v ^= v >> 21
    v ^= v << 37
    v ^= v >> 4
    v *= 4768777513237032717
    v ^= v << 20
    v ^= v >> 41
    v ^= v << 5
    return v


# def unpack_vbits(msg, pattern):
#     """
#     将msg转化成二进制数据
#     :param msg: 输入数据
#     :param pattern: 进行哈希变换时的数据位
#     :return: 二进制数据
#     """
#     vbits = []
#     for word in msg:
#         mask = 1 << 7
#         for i in range(8):
#             vbits.append((word & mask) >> (7 - i))
#             mask = mask >> 1
#     for i in range(pattern - 1):
#         vbits.append(0)
#     return vbits


def unpack_vbits(msg, pattern):
    """
    将msg转化成二进制数据, little-endian
    :param msg: 输入数据
    :param pattern: 进行哈希变换时的数据位
    :return: 二进制数据
    """
    vbits = []
    for word in msg:
        mask = 1
        for _ in range(pattern):
            vbits.append((word & mask))
            word = word >> 1
    # for i in range(pattern - 1):
    #     vbits.append(0)
    return vbits


def pack_vbits(vbits, pattern):
    """
    将二进制数据转化成msg, little-endian
    :param vbits: 输入数据
    :param pattern: 进行哈希变换时的数据位
    :return: 二进制数据
    """
    msg = []
    for i in range(0, len(vbits), pattern):
        temp_word = vbits[i:i + pattern]
        temp_msg = 0
        for _ in range(pattern):
            try:
                temp_msg += temp_word[_] * (1 << _)
            except IndexError:
                # print('error:%d'%_)
                # print(len(vbits))
                pass
        msg.append(temp_msg)
    return msg


def KL_dist(ratio1, ratio2):
    """
    计算两个正态分布之间的距离
    :param ratio1: 正态分布，即ATCG的比例关系
    :param ratio2: 同ratio1
    :return: KL散度，即两个正态分布之间的距离
    """
    eps = 0.001
    tweaked_ratio2 = [r + eps for r in ratio2]
    tweaked_ratio2 = [r / sum(tweaked_ratio2) for r in tweaked_ratio2]
    return sum(r1 * np.log(r1 / r2) if r1 > 0.0 else 0
               for (r1, r2) in zip(ratio1, tweaked_ratio2))


def wasserstein_dist(refer: List[float], input: List[float]):
    """
    计算两个正态分布之间的距离
    :param refer: 正态分布，即ATCG的比例关系
    :param input: 同ratio1
    :return: wasserstein_distance
    """
    p = np.array(refer)
    q = np.array(input)
    emd = wasserstein_distance(p, q)
    return emd


def js_dist(refer: List[float], input: List[float]):
    res = distance.jensenshannon(refer, input)
    return res


def calgc(ratio: dict):
    """
    计算GC含量
    :param ratio: ATCG比例
    :return: GC含量
    """
    gccontent = (ratio['G'] + ratio['C']) / sum(
        [ratio[l] for l in ratio.keys()])
    return gccontent


def calindex(ratio, resolution: int):
    """
    返回ratio对应的哈希值，用于统计是否访问过该ratio
    :param ratio: ATCG比例
    :param resolution: ATCG加起来的值
    :return: ratio对应的哈希值
    """
    return ratio['A'] * (resolution**3) + ratio['C'] * (
        resolution**2) + ratio['G'] * resolution + ratio['T']


class HEDGES:

    def __init__(self, resolution, sigma, step):
        self.pattern = 8  # 8个二进制位对应一个letter
        self.step = step  # 每次滑窗步进的二进制位数，--> code rate
        self.pattern_mask = (1 << 7 * self.pattern) - 1  # 对应pattern的mask
        self.MAX_SEQ = 250000  # 允许的最大Hypothesis数量 250000
        self.PENALTY_THRESHOLD = 10  # 10 控制解码退出，惩罚太高直接退出
        self.NODE_THRESHOLD = 50 # 50  # 解码树每层最多节点数
        self.BUCKET_NUM = 13 # 10

        self.mod = 0
        self.duplicate = 1
        self.resolution = resolution
        self.sigma = sigma
        self.alphabet = []
        self.alphabet_set(self.resolution, self.sigma)

    def alphabet_set(self, resolution: int, sigma: int):
        """
        计算合法的全部letter, sigma mod
        ACGTMKRY
        :param resolution: ATCG加起来的值
        :return: None
        """
        base_A = composite_letter(code_pattern(resolution)['A'])
        base_C = composite_letter(code_pattern(resolution)['C'])
        base_G = composite_letter(code_pattern(resolution)['G'])
        base_T = composite_letter(code_pattern(resolution)['T'])
        base_M = composite_letter(code_pattern(resolution)['M'])
        base_K = composite_letter(code_pattern(resolution)['K'])
        base_R = composite_letter(code_pattern(resolution)['R'])
        base_Y = composite_letter(code_pattern(resolution)['Y'])

        if sigma == 4:
            self.alphabet = [base_A, base_C, base_G, base_T]
        elif sigma == 6:
            self.alphabet = [base_A, base_C, base_G, base_T, base_M, base_K]
        elif sigma == 8:
            self.alphabet = [
                base_A, base_C, base_G, base_T, base_M, base_K, base_R, base_Y
            ]
        self.mod = sigma
        

    def alphabet_show(self, n=10):
        """
        show alphabet
        :param n: 每行显示的letter个数
        :return: None
        """
        count = 0
        for i in range(self.mod):
            print(self.alphabet[i], end=' ')
            count += 1
            if count == n:
                print()
                count = 0
        print()
        print(f'resolution: {self.resolution}\t\tsigma: {self.sigma}')

    def dnacallowed(self, codetext):
        """
        Limit homopolymer length and GC content
        """
        all_base = ['A','C','G','T','M','K','R','Y']
        related_alphabet = [[[2,0,0,0],[1,1,0,0],[1,0,1,0]], # A,M,R
                            [[0,2,0,0],[1,1,0,0],[0,1,0,1]], # C,M,Y
                            [[0,0,2,0],[0,0,1,1],[1,0,1,0]], # G,K,R
                            [[0,0,0,2],[0,0,1,1],[0,1,0,1]]] # T,K,Y
        maxgcratio = 0.55 # 0.55
        mingcratio = 0.45 # 0.45
        max_runs = 4 # 4
        max_window = 160
        is_homopolymer = False
 
        alphabetletter = [[2, 0, 0, 0], [0, 2, 0, 0], [0, 0, 2, 0], [0, 0, 0, 2],
                          [1, 1, 0, 0], [0, 0, 1, 1], [1, 0, 1, 0], [0, 1, 0, 1]]
        

        if len(codetext) > max_runs-1:
            alphabet_to_delete = []
            for sub_related_alphabet in related_alphabet:
                tmp = []
                for item in codetext[-max_runs:]:
                    apb = [item.ratio_dict[x] for x in all_base[:4]]
                    tmp.append(apb in sub_related_alphabet)
                is_homopolymer = all(tmp)
                if is_homopolymer:
                    alphabet_to_delete = []
                    alphabet_to_delete = [alphabetletter.index(i) for i in sub_related_alphabet]
                    break
        alphabet_to_delete = alphabet_to_delete if is_homopolymer else []
        window = len(codetext) if len(codetext)<= max_window else max_window
        if len(codetext) > max_runs-1:
            alphabet_count_in_window = [0] * 8
            for code_idx in range(-window , 0):
                current_code_base_ratio = [codetext[code_idx].ratio_dict[x] for x in all_base[:4]]
                for _ in range(len(alphabetletter)):
                    alphabet_count_in_window[_] += 1 if current_code_base_ratio == alphabetletter[_] else 0
            
            gc_count = (alphabet_count_in_window[1]+alphabet_count_in_window[2])
            my_count = (alphabet_count_in_window[4]+alphabet_count_in_window[7])
            kr_count = (alphabet_count_in_window[5]+alphabet_count_in_window[6])

            if mingcratio *len(codetext) > gc_count and self.sigma == 4:
                alphabet_to_delete += [0,3]
            if maxgcratio *len(codetext) < gc_count and self.sigma == 4:
                alphabet_to_delete += [1,2]
            
            if my_count < kr_count:
                alphabet_to_delete += [5,6]
            elif my_count > kr_count:
                alphabet_to_delete += [4,7]
            
            alphabet = self.alphabet
            alphabet = [item for index, item in enumerate(alphabet) if index not in list(set(alphabet_to_delete))]
            mod = len(alphabet)
        else:
            alphabet = self.alphabet[:self.sigma]
            mod = len(alphabet)
        
        return mod, alphabet


    def encode(self, vbits, change_alphabet=False):
        """
        编码
        if DEBUG0: msg为全0数据
        if DEBUG1: msg为全1数据
        :param msg: input 数据
        :return: composite hedges编码后letter
        """
        flag_high = True
        if DEBUG0:
            vbits = [0] * 256
        if DEBUG1:
            vbits = [1] * 256
        # vbits = [1, 0] * 128
        # print(vbits)
        nm = len(vbits)
        # print(nm)
        if nm > self.MAX_SEQ:
            raise ValueError('encode: MAXSEQ too small')
        codetext = []
        rawtext = []
        prev_bits = 0
        k = 0
        while k < nm:
            bit = 0
            for _ in range(self.step):
                try:
                    bit = (bit << 1) + vbits[k]
                    k += 1
                except:
                    bit = bit << 1
                    k += 1
            prev_bits = ((prev_bits << self.step) + bit) & self.pattern_mask
            if self.sigma == 4 or self.sigma == 8:
                mod, alphabet = self.dnacallowed(codetext)
            else:
                mod, alphabet = self.mod, self.alphabet
            # mod, alphabet = self.mod,self.alphabet
            # print(mod , alphabet)
            if change_alphabet:
                alphabet = alphabet[1:] + [alphabet[0]]
            digest_code = self.digest(prev_bits, k - 1, mod)
            # print('digest: %d, prev_bits:%d, index:%d' % (digest_code, prev_bits,  k - 1))
            rawtext.append(digest_code)
            letter = alphabet[digest_code]
            # if gc_content = 50%
            temp = deepcopy(letter)
            codetext.append(temp)
            # if len(codetext) % 40 == 0 and len(codetext) > 0 and len(codetext)!= 160:
            #     codetext += [self.alphabet[0],self.alphabet[0],self.alphabet[2],self.alphabet[2],self.alphabet[1]]

        return codetext

    def decode(self, data: List[composite_letter]):
        """
        通过遍历每一步的全部可能性，计算可能生成的当前letter可能性与data比较，并返回概率最大的数据
        :param data: letter list
        :return: decode后的二进制数据及其的得分
        """
        hypotree = HypothesisTree(step=self.step)
        hypotree_plist = [hypotree.root]
        hypotree_qlist = []
        index = 0
        
        for _ in range(len(data)):
            hypotree_qlist = []
            index += self.step
            for hypotree_p in hypotree_plist:
                penalty_list = []
                hypo_prev_bits_list = []
                
                for bit in range(1 << self.step):                   
                    hypo_prev_bits = ((hypotree_p.prev_bits << self.step) +
                                      bit) & self.pattern_mask

                    mod, alphabet = self.dnacallowed(data[:_])
                    hypo_code = self.digest(hypo_prev_bits, index - 1,
                                            mod)
                    
                    hypo_letter = alphabet[hypo_code]

                    penalty = js_dist(hypo_letter.ratio_cal(),
                                      data[_].ratio_cal())
                    
                    penalty_list.append(penalty)
                    hypo_prev_bits_list.append(hypo_prev_bits)

                hypotree.add_children(node=hypotree_p,
                                    penalty_list=penalty_list,
                                    prev_bits_list=hypo_prev_bits_list)
                # delete high penalty node
                for child in hypotree_p.children:
                    if child.penalty < self.PENALTY_THRESHOLD:
                        hypotree_qlist.append(child)

            if len(hypotree_qlist) > self.NODE_THRESHOLD:
                bucket_list = [[] for _ in range(self.BUCKET_NUM)]
                bucket_size = self.PENALTY_THRESHOLD / self.BUCKET_NUM
                for i in range(len(hypotree_qlist)):
                    bucket_list[int(hypotree_qlist[i].penalty //
                                    bucket_size)].append(hypotree_qlist[i])
                cnt = 0
                templist = []
                min_bucket_size = len(bucket_list[0])
                for i in range(self.BUCKET_NUM):
                    templist.append(bucket_list[i])
                    cnt += len(bucket_list[i])
                    if cnt >= self.NODE_THRESHOLD:
                        break
                templist = [
                    element for sublist in templist for element in sublist
                ]
                hypotree_qlist = templist[:self.NODE_THRESHOLD
                                          if min_bucket_size < self.
                                          NODE_THRESHOLD else min_bucket_size]

            hypotree_plist = hypotree_qlist

        del hypotree_qlist
        return hypotree_plist


    def hypo_backtrack(self, hypotree_list):
        min_penalty = np.inf
        chosen = []
        res = []

        for hypo in hypotree_list:
            if hypo.penalty < min_penalty:
                min_penalty = hypo.penalty
        for i in range(len(hypotree_list)):
            if hypotree_list[i].penalty == min_penalty:
                chosen.append(i)

        for i in range(len(chosen)):
            stack = []
            stack_n = []
            p = hypotree_list[chosen[i]]
            while p.parent != None:
                stack.append(p.bits)
                p = p.parent
            stack.reverse()

            for num in stack:
                if self.step != 1:
                    stack_n.append((num // (self.step)) % self.step)
                    stack_n.append(num % self.step)
                else:
                    stack_n.append(num)
            # print(stack_n)
            # stack_n.reverse()
            # print(stack_n)
            res.append(stack_n)

        return res, min_penalty

    def digest(self, bits, k: int, mod: int):
        """
        计算并返回对应的哈希值
        :param bits: 前k位二进制数据
        :param k: 当前二进制位的index
        :param mod: 模
        :return:
        """
        return int(ran_hash(bits + (k << self.pattern)) % mod)

    def to_combine(self, encode_res: List[dict]) -> str:
        # res = []
        tempstr = ''
        cnt = 0
        for item in encode_res:
            for key in code_pattern(self.resolution):
                if item.ratio_dict == code_pattern(self.resolution)[key]:
                    tempstr += key
                    cnt += 1
                    # if cnt % fragment == 0:
                    #     res.append(tempstr)
                    #     tempstr = ''
        # if tempstr != '':
        #     res.append(tempstr)
        return tempstr

    def combine_to_letters(self, string: str) -> List[composite_letter]:
        data = []
        for i in range(len(string)):
            data.append(
                composite_letter(code_pattern(self.resolution)[string[i]]))
        return data

def eight2four(seq: str, num_copies):
    alphabet_dict = {'A':['A'], 'C':['C'], 'T':['T'], 'G':['G'],'M':['A','C'], 'K':['G','T'], 'R':['A','G'], 'Y':['C','T']}
    gc_window = 10
    gc_max = 0.525
    gc_min = 0.475
    copies = [''] * num_copies
    for idx in range(len(seq)):
        if idx + 1 >= gc_window:
            current_base = seq[idx]
            current_alphabet_list = alphabet_dict[current_base]
            if len(current_alphabet_list) == 2:
                initial_list = current_alphabet_list * (num_copies//2)
                all_permutations = list(itertools.permutations(initial_list))
                unique_permutations = set(all_permutations)
                result = [list(perm) for perm in unique_permutations]
                all_alphabet_list = result
                dist = 9999
                for alphabet_mode in all_alphabet_list:
                    temp = ['']*num_copies
                    for _ in range(num_copies):
                        temp[_] = copies[_] + str(alphabet_mode[_])
                    gc_ratio = [(i[:].count('C')+i[:].count('G'))/len(i[:]) for i in temp]
                    dd = distance.jensenshannon([0.5]*num_copies,gc_ratio)
                    dd = math.sqrt(sum([(i-0.5)**2 for i in gc_ratio]))
                    if dd < dist:
                        gc_ok_mode = alphabet_mode
                        dist = dd
                for _ in range(num_copies):
                    copies[_] = copies[_] + str(gc_ok_mode[_])
            elif len(current_alphabet_list) == 1:
                all_alphabet_list = current_alphabet_list * num_copies
                for _ in range(num_copies):
                    copies[_] += str(all_alphabet_list[_])
        else:
            current_base = seq[idx]
            current_alphabet_list = alphabet_dict[current_base]
            if len(current_alphabet_list) == 2:
                initial_list = current_alphabet_list * (num_copies//2)
                all_permutations = list(itertools.permutations(initial_list))
                unique_permutations = set(all_permutations)
                result = [list(perm) for perm in unique_permutations]
                all_alphabet_list = result
                for alphabet_mode in all_alphabet_list:
                    temp = ['']*num_copies
                    for _ in range(num_copies):
                        temp[_] = copies[_] + str(alphabet_mode[_])
                    gc_ratio = [(i[:].count('C')+i[:].count('G'))/len(i[:]) for i in temp]
                    gc_large = [int(gc_ratio_i >= gc_max) for gc_ratio_i in gc_ratio].count(1) > 0
                    gc_small = [int(gc_ratio_i <= gc_min) for gc_ratio_i in gc_ratio].count(1) > 0
                    gc_ok = not(gc_large) and not(gc_small)
                    if gc_ok:
                        copies = temp
                        break
                copies = temp
            elif len(current_alphabet_list) == 1:
                all_alphabet_list = current_alphabet_list * num_copies
                for _ in range(num_copies):
                    copies[_] += str(all_alphabet_list[_])
    return copies

def make_gc_inrange(seq: str):
    def check_homopolymer(seq: str):
        alphabet = ['A','C','T','G']
        homopolymer_list = [seq.count(_*5)== 0 for _ in alphabet]
        return all(homopolymer_list)
    
    ori_seq = seq
    gc_min = 0.45
    gc_max = 0.6
    allgc_index = []
    gc_ratio = (ori_seq[:].count('C')+ori_seq[:].count('G'))/len(ori_seq[:])
    choice_list = [['A','T'],['C','G']]
    
    if gc_ratio > gc_max:
        gc = 1 
    elif gc_ratio < gc_min: 
        gc = 0
    else:
        gc = 'ok'
        
    if gc == 'ok':
        return seq
    for i, c in enumerate(seq):
        if c in choice_list[gc]:
            allgc_index.append(i)
    while True:
        random_index = random.choice(allgc_index)
        new_seq = ori_seq[:random_index] + random.choice(choice_list[abs(gc-1)]) + ori_seq[random_index+1:]
        gc_ratio = (new_seq[:].count('C')+new_seq[:].count('G'))/len(new_seq[:])
        if not(gc_ratio > gc_min and gc_ratio < gc_max):
            ori_seq = new_seq
            continue
        if check_homopolymer(new_seq):
            break
        else:
            ori_seq = seq
            continue
    seq = new_seq
    return seq

def four2eight(seqs):
    ratio_list = []
    ratio_div = []
    seqin8 = []
    for m in range(min([len(item) for item in seqs])):
        temp = []
        for n in range(len(seqs)):
                temp.append(seqs[n][m])
        lis = str(temp)
        ratio_list.append([lis.count('A'),lis.count('C'),lis.count('G'),lis.count('T')])
        ratio_div.append(lis.count('-'))
    for ratio in ratio_list:
        if ratio_div[ratio_list.index(ratio)] > int(len(seqs)*0.8):
            continue
        res_M = distance.jensenshannon([1,1,0,0],ratio)
        res_K = distance.jensenshannon([0,0,1,1],ratio)
        res_R = distance.jensenshannon([1,0,1,0],ratio)
        res_Y = distance.jensenshannon([0,1,0,1],ratio)
        res_A = distance.jensenshannon([1,0,0,0],ratio)
        res_C = distance.jensenshannon([0,1,0,0],ratio)
        res_G = distance.jensenshannon([0,0,1,0],ratio)
        res_T = distance.jensenshannon([0,0,0,1],ratio)
        res = [res_M,res_K,res_R,res_Y,res_A,res_C,res_G,res_T]
        min_index = res.index(min(res))
        if min_index == 0:
            seqin8.append('M')
        if min_index == 1:
            seqin8.append('K')
        if min_index == 2:
            seqin8.append('R')
        if min_index == 3:
            seqin8.append('Y')
        if min_index == 4:
            seqin8.append('A')
        if min_index == 5:
            seqin8.append('C')
        if min_index == 6:
            seqin8.append('G')
        if min_index == 7:
            seqin8.append('T')
    return ''.join(seqin8)
        

if __name__ == '__main__':
    resolution = 2
    refs = [str(i.seq) for i in SeqIO.parse('./reference/zxy_sustech_seqs.fasta','fasta')]
    
    sustech_pic_ref = r'./reference/sustech_logo.jpg'
    sustech_txt_ref = r'./reference/sustech_introduction.txt'
    txt_ref = open(sustech_txt_ref, 'rb').read()
    pic_ref = open(sustech_pic_ref, 'rb').read()
    txt_ref = [txt_ref[i:i+18] for i in range(0, len(txt_ref), 18)]
    pic_ref = [pic_ref[i:i+18] for i in range(0, len(pic_ref), 18)]
    source_file_ref = txt_ref + pic_ref

    sigma = 4
    step = 1
    N, K = 20, 18

    # Initialize counters
    decode_fail_num = 0
    total_mismatched_bytes = 0
    total_bytes = 0
    total_segments = 0
    total_failed_segments = 0

    sustech_pic_output = f'./sustech_logo_decode.jpg'
    sustech_txt_output = f'./sustech_text_decode.txt'
    rsc = RSCodec(N - K, nsize=N)

    hedges = HEDGES(resolution=resolution, sigma=sigma, step=step)
    hedges.alphabet_show()

    consensuses = [str(i.seq) for i in SeqIO.parse(INPUT_FASTA_PATH,'fasta')]
    consensuses_noprimers = [i[57:-11] for i in consensuses]
    payloads = [i[:40] + i[45:85] + i[90:130] + i[135:] for i in consensuses_noprimers]
    decode_st_ticks = time.time()
    pic_decode_res = []
    txt_decode_res = []
    fail_index = []
    
    with tqdm.trange(len(payloads), desc='decoding...') as tbar:
        for i in range(len(payloads)):
            total_segments += 1  # Count total number of segments processed
            data = hedges.combine_to_letters(payloads[i])
            ref = refs[i][57:-11]
            ref = ref[:40]+ref[45:85]+ref[90:130]+ref[135:]

            if len(payloads[i]) != len(ref):
                print(f'index {i} LEN not match ❌, ref length is {len(ref)}, payload length is {len(payloads[i])}'),
                decode_fail_num += 1
                total_failed_segments += 1
                fail_index.append(i)
                continue

            error_count = sum([1 for j in range(len(payloads[i])) if payloads[i][j] != ref[j]])

            res = hedges.decode(data)
            decode_vbits, min_penalty = hedges.hypo_backtrack(res)

            if decode_vbits == []:
                rs_decode_res = bytearray(18)
            
            decode_success = False
            for idx, item in enumerate(decode_vbits):
                try:
                    msg = pack_vbits(item, hedges.pattern)
                    rs_decode_res = rsc.decode(bytearray(msg))[0]
                    decode_success = True
                    break
                except:
                    rs_decode_res = bytearray(18)
                    continue

            # Compare decoded result with reference byte by byte
            ref_bytes = source_file_ref[i]
            total_bytes += len(ref_bytes)
            
            if decode_success:
                # Optional, only when we have a reference to compare with.
                mismatched_bytes = sum(1 for a, b in zip(rs_decode_res, ref_bytes) if a != b)
                if mismatched_bytes > 0:
                    total_mismatched_bytes += mismatched_bytes
                    mismatch_percentage = (mismatched_bytes / len(ref_bytes)) * 100
                    print(f'index {i} not match with reference! {mismatched_bytes}/{len(ref_bytes)} bytes mismatched ({mismatch_percentage:.2f}%) ❌')
                    decode_fail_num += 1
                    total_failed_segments += 1
                    fail_index.append(i)
            else:
                # If decode failed completely
                total_mismatched_bytes += len(ref_bytes)
                print(f'index {i} decode failed completely! {len(ref_bytes)}/{len(ref_bytes)} bytes mismatched (100.00%) ❌')
                decode_fail_num += 1
                total_failed_segments += 1
                fail_index.append(i)

            if i <= 54:
                txt_decode_res.append(rs_decode_res)
            if i > 54:
                pic_decode_res.append(rs_decode_res)
            
            tbar.update()

        decode_end_ticks = time.time()
        print('decode time: %.3f' % (decode_end_ticks - decode_st_ticks))

        # Write decoded results to files
        with open(sustech_pic_output, 'wb') as f:
            for item in pic_decode_res:
                f.write(item)
        with open(sustech_txt_output, 'wb') as f:
            for item in txt_decode_res:
                f.write(item)

        # Compare final decoded files with reference files
        with open(sustech_pic_ref, 'rb') as f:
            original_pic = f.read()
        with open(sustech_txt_ref, 'rb') as f:
            original_txt = f.read()
        with open(sustech_pic_output, 'rb') as f:
            decoded_pic = f.read()
        with open(sustech_txt_output, 'rb') as f:
            decoded_txt = f.read()

        # Calculate accuracy for picture
        pic_matched_bytes = sum(1 for a, b in zip(original_pic, decoded_pic) if a == b)
        pic_accuracy = (pic_matched_bytes / len(original_pic)) * 100

        # Calculate accuracy for text
        txt_matched_bytes = sum(1 for a, b in zip(original_txt, decoded_txt) if a == b)
        txt_accuracy = (txt_matched_bytes / len(original_txt)) * 100

        # Calculate weighted average accuracy
        total_bytes = len(original_pic) + len(original_txt)
        weighted_accuracy = (pic_accuracy * len(original_pic) + txt_accuracy * len(original_txt)) / total_bytes

        print('\nDecoding Summary:')
        print(f'Total segments processed: {total_segments}')
        print(f'Failed segments: {total_failed_segments} ({(total_failed_segments/total_segments)*100:.2f}%)')
        print(f'Total bytes processed during decoding: {total_bytes}')
        print(f'Total mismatched bytes during decoding: {total_mismatched_bytes}/{total_bytes} ({(total_mismatched_bytes/total_bytes)*100:.2f}%)')
        print('\nFile Comparison Results:')
        print(f'Picture accuracy: {pic_matched_bytes}/{len(original_pic)} ({pic_accuracy:.2f}%)')
        print(f'Text accuracy: {txt_matched_bytes}/{len(original_txt)} ({txt_accuracy:.2f}%)')
        print(f'Weighted average accuracy: {weighted_accuracy:.2f}%')
        print(f'Failed indices: {fail_index}')
        
        # Save failed indices
        with open(os.environ['FAILED_INDEX_TXT'], 'w') as f:
            f.write(str(fail_index))
        print('decode finished!')