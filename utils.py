import yaml
import re
import selfies as sf


class CharVocab:
    def __init__(self, vocab_path):
        self.name = "char"

        # load the pre-computed vocabulary
        with open(vocab_path, 'r') as f:
            self.vocab = yaml.full_load(f)

        # a dictionary to map integer back to SMILES
        # tokens for sampling
        self.int2tocken = {}
        for token, num in self.vocab.items():
            self.int2tocken[num] = token

        # a hashset of tokens for O(1) lookup
        self.tokens = self.vocab.keys()

    def tokenize_smiles(self, smiles):
        """
        Takes a SMILES string and returns a list of tokens.
        Atoms with 2 characters are treated as one token. The 
        logic references this code piece:
        https://github.com/topazape/LSTM_Chem/blob/master/lstm_chem/utils/smiles_tokenizer2.py
        """
        n = len(smiles)
        tokenized = ['<sos>']
        i = 0

        # process all characters except the last one
        while (i < n - 1):
            # procoss tokens with length 2 first
            c2 = smiles[i:i + 2]
            if c2 in self.tokens:
                tokenized.append(c2)
                i += 2
                continue

            # tokens with length 2
            c1 = smiles[i]
            if c1 in self.tokens:
                tokenized.append(c1)
                i += 1
                continue

            raise ValueError(
                "Unrecognized charater in SMILES: {}, {}".format(c1, c2))

        # process last character if there is any
        if i == n:
            pass
        elif i == n - 1 and smiles[i] in self.tokens:
            tokenized.append(smiles[i])
        else:
            raise ValueError(
                "Unrecognized charater in SMILES: {}".format(smiles[i]))

        tokenized.append('<eos>')

        tokenized = [self.vocab[token] for token in tokenized]
        return tokenized

    def combine_list(self, smiles):
        return "".join(smiles)


class RegExVocab:
    def __init__(self, vocab_path):
        self.name = "regex"

        # load the pre-computed vocabulary
        with open(vocab_path, 'r') as f:
            self.vocab = yaml.full_load(f)

        # a dictionary to map integer back to SMILES
        # tokens for sampling
        self.int2tocken = {}
        for token, num in self.vocab.items():
            if token == "R":
                self.int2tocken[num] = "Br"
            elif token == "L":
                self.int2tocken[num] = "Cl"
            else:
                self.int2tocken[num] = token

    def tokenize_smiles(self, smiles):
        """Takes a SMILES string and returns a list of tokens.
        This will swap 'Cl' and 'Br' to 'L' and 'R' and treat
        '[xx]' as one token."""
        regex = '(\[[^\[\]]{1,6}\])'
        smiles = self.replace_halogen(smiles)
        char_list = re.split(regex, smiles)

        tokenized = ['<sos>']

        for char in char_list:
            if char.startswith('['):
                tokenized.append(char)
            else:
                chars = [unit for unit in char]
                [tokenized.append(unit) for unit in chars]
        tokenized.append('<eos>')

        # convert tokens to integer tokens
        tokenized = [self.vocab[token] for token in tokenized]

        return tokenized

    def replace_halogen(self, string):
        """Regex to replace Br and Cl with single letters"""
        br = re.compile('Br')
        cl = re.compile('Cl')
        string = br.sub('R', string)
        string = cl.sub('L', string)

        return string

    def combine_list(self, smiles):
        return "".join(smiles)


class SELFIESVocab:
    def __init__(self, vocab_path):
        self.name = "selfies"

        # load the pre-computed vocabulary
        with open(vocab_path, 'r') as f:
            self.vocab = yaml.full_load(f)

        self.int2tocken = {value: key for key, value in self.vocab.items()}

    def tokenize_smiles(self, smiles):
        """convert the smiles to selfies, then return 
        integer tokens."""
        ints = [self.vocab['<sos>']]

        encoded_selfies = sf.encoder(smiles)
        
        selfies_list = list(sf.split_selfies(encoded_selfies))
        for token in selfies_list:
            ints.append(self.vocab[token])
        ints.append(self.vocab['<eos>'])

        return ints

    def combine_list(self, selfies):
        return "".join(selfies)
