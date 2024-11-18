import copy


class composite_letter:

    def __repr__(self):
        return repr(self.ratio_dict)

    def __str__(self):
        return '{}'.format(''.join(
            (str(self.ratio_dict[i]) + i for i in self.basic_alphabet)))

    def __init__(self, ratio_dict, basic_alphabet=None):
        if basic_alphabet is None:
            basic_alphabet = ['A', 'C', 'G', 'T']
        self.basic_alphabet = basic_alphabet
        self.ratio_dict = ratio_dict

    def ratio_cal(self):
        """
        计算ACGT比例
        :return: ACGT比例
        """
        return [
            self.ratio_dict['A'], self.ratio_dict['C'], self.ratio_dict['G'],
            self.ratio_dict['T']
        ]

    def trans2base(self):
        """
        将ACGT比例转化为ACGTMKRY
        :return: 字母
        """
        alphabet_dict = {'A':[2,0,0,0],'C':[0,2,0,0],'G':[0,0,2,0],'T':[0,0,0,2],
                         'M':[1,1,0,0],'K':[0,0,1,1],'R':[1,0,1,0],'Y':[0,1,0,1]}
        ratio = self.ratio_cal()
        for key, value in alphabet_dict.items():
            if ratio == value:
                return key
        # return alphabet_dict[ratio]

class Hypothesis:

    def __init__(self,
                 step,
                 pattern_mask,
                 MAX_SEQ,
                 index=-1,
                 pred=None,
                 msg_bit=None):
        self.pred = pred
        self.msg_bit = msg_bit
        self.index = index  # bit index
        self.score = -1
        self.MAX_SEQ = MAX_SEQ
        if self.pred:
            self.prev_bits = (
                (self.pred.prev_bits << step) + self.msg_bit) & pattern_mask
        else:
            self.prev_bits = self.msg_bit
        # self.letter = None

    def __str__(self):
        return f'bit:{self.msg_bit} index:{self.index} score:{self.score}'

    def init_from_predecessor(self, pred, skew):
        self.pred = pred
        self.index = pred.index + 1
        if self.index > self.MAX_SEQ:
            raise ValueError('init_from_predecessor: MAX_SEQ too small')
        


class HypothesisNode:

    def __init__(self, step, prev_bits, bits=0, penalty=0, parent=None, children=None, letters=None):
        self.step = step
        self.penalty = penalty
        self.parent = parent
        self.prev_bits = prev_bits
        self.bits = bits
        # self.letters = letters
        if children:
            self.children = copy.deepcopy(children)
        else:
            self.children = [None] * (1 << self.step)


class HypothesisTree:

    def __init__(self, step):
        self.step = step
        # self.root = HypothesisNode(step=step,prev_bits=0, penalty=0, parent=None, children=None, letters=[])
        self.root = HypothesisNode(step=step,prev_bits=0, penalty=0, parent=None, children=None)
        self.cur_node = self.root

    # def add_children(self, node: HypothesisNode, penalty_list, prev_bits_list, letters_list):
    def add_children(self, node: HypothesisNode, penalty_list, prev_bits_list):
        for i in range(1 << self.step):
            child = HypothesisNode(step=self.step,
                                   prev_bits=prev_bits_list[i],
                                   bits=i,
                                   penalty=(node.penalty + penalty_list[i]),
                                   parent=node,
                                   children=None)
            # child = HypothesisNode(step=self.step,
            #                        prev_bits=prev_bits_list[i],
            #                        bits=i,
            #                        penalty=(node.penalty + penalty_list[i]),
            #                        parent=node,
            #                        children=None,
            #                        letters= letters_list[i])
            node.children[i] = child

    def delete_children(self, node: HypothesisNode):
        parent = node.parent
        for idx, child in enumerate(parent.children):
            if child == node:
                parent.children[idx] = None


def code_pattern(resolution: int) -> dict:
    CODE_PATTERN = {
        'A': {
            'A': resolution,
            'C': 0,
            'G': 0,
            'T': 0
        },
        'C': {
            'A': 0,
            'C': resolution,
            'G': 0,
            'T': 0
        },
        'G': {
            'A': 0,
            'C': 0,
            'G': resolution,
            'T': 0
        },
        'T': {
            'A': 0,
            'C': 0,
            'G': 0,
            'T': resolution
        },
        'M': {
            'A': resolution // 2,
            'C': resolution // 2,
            'G': 0,
            'T': 0
        },
        'K': {
            'A': 0,
            'C': 0,
            'G': resolution // 2,
            'T': resolution // 2
        },
        'R': {
            'A': resolution // 2,
            'C': 0,
            'G': resolution // 2,
            'T': 0
        },
        'Y': {
            'A': 0,
            'C': resolution // 2,
            'G': 0,
            'T': resolution // 2
        }
    }
    return CODE_PATTERN
