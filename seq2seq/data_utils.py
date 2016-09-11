# Copyright 2015 The TensorFlow Authors. All Rights Reserved.
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
# ==============================================================================

"""Utilities for downloading data from WMT, tokenizing, vocabularies."""
from __future__ import absolute_import
from __future__ import division
from __future__ import print_function


from tensorflow.python.platform import gfile

# Special vocabulary symbols - we always put them at the start.
_PAD = "_PAD"
_GO = "_GO"
_EOS = "_EOS"
_START_VOCAB = [_PAD, _GO, _EOS]

PAD_ID = 0
GO_ID = 1
EOS_ID = 2


def write_vocabulary(vocabulary_path, data_path, tokenizer=list):
  """Create vocabulary file (if it does not exist yet) from data file.

  Data file is assumed to contain one sequence per line. Each sentence is
  tokenized.
  We write it to vocabulary_path in a one-token-per-line format, so that later
  token in the first line gets id=0, second line gets id=1, and so on.

  Args:
    vocabulary_path: path where the vocabulary will be created.
    data_path: data file that will be used to create vocabulary.
    tokenizer: a function to use to tokenize each seq
  """
  if not gfile.Exists(vocabulary_path):
    print("Creating vocabulary %s from data %s" % (vocabulary_path, data_path))
    vocab = {}
    with gfile.GFile(data_path, mode="r") as f:
      counter = 0
      for line in f:
        counter += 1
        if counter % 5000 == 0:
          print("  processing line %d" % counter)
        tokens = line.strip().split()
        assert len(tokens) == 1
        tokens = tokenizer(tokens[0]) #split into single characters
        for word in tokens:
          if word in vocab:
            vocab[word] += 1
          else:
            vocab[word] = 1

      vocab_list = _START_VOCAB + sorted(vocab)
      with gfile.GFile(vocabulary_path, mode="w") as vocab_file:
        for w in vocab_list:
          vocab_file.write(w + "\n")


def read_vocabulary(vocabulary_path):
  """Initialize vocabulary from file.

  We assume the vocabulary is stored one-item-per-line, so a file:
    A
    C
  will result in a vocabulary {"A": 0, "C": 1}, and this function will
  also return the reversed-vocabulary ["A", "C"].

  Args:
    vocabulary_path: path to the file containing the vocabulary.

  Returns:
    a pair: the vocabulary (a dictionary mapping string to integers), and
    the reversed vocabulary (a list, which reverses the vocabulary mapping).

  Raises:
    ValueError: if the provided vocabulary_path does not exist.
  """
  if gfile.Exists(vocabulary_path):
    rev_vocab = []
    with gfile.GFile(vocabulary_path, mode="r") as f:
      rev_vocab.extend(f.readlines())
    rev_vocab = [line.strip() for line in rev_vocab]
    vocab = dict([(x, y) for (y, x) in enumerate(rev_vocab)])
    return vocab, rev_vocab
  else:
    raise ValueError("Vocabulary file %s not found.", vocabulary_path)


def seq_to_token_ids(seq, vocabulary, tokenizer=list):
  """Convert a string to list of integers representing token-ids.

  For example, a sequence "ACGCT" may become tokenized into
  ["A", "C", "G", "C", "T"] and with vocabulary {"A": 1, "C": 2,
  "G": 3, "T": 4} this function will return [1, 2, 3, 2, 4].

  Args:
    seq: the seq in bytes format to convert to token-ids.
    vocabulary: a dictionary mapping tokens to integers.
    tokenizer: a function to use to tokenize each sentence;

  Returns:
    a list of integers, the token-ids for the sentence.
  """

  words = tokenizer(seq)
  return [vocabulary[w] for w in words]


def data_to_token_ids(data_path, target_path, vocabulary_path,
                      tokenizer=list):
  """Tokenize data file and turn into token-ids using given vocabulary file.

  This function loads data line-by-line from data_path, calls the above
  seq_to_token_ids, and saves the result to target_path. See comment
  for seq_to_token_ids on the details of token-ids format.

  Args:
    data_path: path to the data file in one-sentence-per-line format.
    target_path: path where the file with token-ids will be created.
    vocabulary_path: path to the vocabulary file.
    tokenizer: a function to use to tokenize each sentence;
  """
  if not gfile.Exists(target_path):
    print("Tokenizing data in %s" % data_path)
    vocab, _ = read_vocabulary(vocabulary_path)
    with gfile.GFile(data_path, mode="r") as data_file:
      with gfile.GFile(target_path, mode="w") as tokens_file:
        counter = 0
        for line in data_file:
          counter += 1
          if counter % 5000 == 0:
            print("  tokenizing line %d" % counter)
          token_ids = seq_to_token_ids(line.strip(), vocab, tokenizer)
          tokens_file.write(" ".join([str(tok) for tok in token_ids]) + "\n")


def prepare_sequence_data(training_file, validation_file, tokenizer=list):
  """Create vocabularies and tokenize data.

  Args:
    training_file: training sequence file.
    validation_file: validation sequence file.
    tokenizer: a function to use to tokenize each data sentence;

  Returns:
    A tuple of 3 elements:
      (1) path to the token-ids for the training data-set,
      (2) path to the token-ids for the validation data-set,
      (3) path to the vocabulary file,
  """
  vocab_path = training_file + ".sequence.vocab"
  training_ids_path = training_file + ".sequence.ids" # trainin file encoded with new IDs
  validation_ids_path = validation_file + ".sequence.ids"  # validation file encoded with new IDs

  # Create vocabularies of the appropriate sizes.
  write_vocabulary(vocab_path, training_file, tokenizer)

  # Create token ids for the training data.
  data_to_token_ids(training_file, training_ids_path, vocab_path, tokenizer)
  data_to_token_ids(validation_file, validation_ids_path, vocab_path, tokenizer)

  return training_ids_path, validation_ids_path, vocab_path