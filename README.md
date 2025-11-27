# SnapGene File Format Converter

This is a reversed engineered parser for SnapGene file formats (SGFF in short) - for dna, rna and protein.

Currently parser partially do it's job, the result is a JSON dictionary.

Ideally this project will create a set of writer and reader scripts with the following scheme:

```mermaid
graph LR;
    A[SGFF] ---> B[Reader];
    B ---> C[JSON];
    C ---> D[Writer];
    D ---> E[SGFF];
```

## ImHex pattern parser

The project also has a pattern parser for ImHex: `snapgene.hexpat`, you may use it for examining binary.

## File structure

The file is a Type-Length-Value format, which also partially has compressed blocks (e.g. history node/tree) with LZMA. (see documentation in the acknowledgment section for some insights)

Currently file structure looks like this:

```bash
# file structure table
# 0: sequence_dna:ascii
# 1: compressed sequence:complex structure (2bpb format)
# 2:
# 3: enzyme_library:mixed (enzyme sites + id?)
# 4:
# 5: primers:xml
# 6: notes:xml
# 7: history tree:xml
# 8: additional sequence properties:xml
# 9: file Description:?
# 10: features:xml
# 11: history_node:tlv container for ((0/1/32)+30)/29
# 12:
# 13: enzyme_info:mixed
# 14: enzyme_custom:xml
# 15:
# 16: sequence trace(legacy): 4 empty bytes
# 17: alignable sequences:xml
# 18: sequence trace:zrt-trace format
# 19: uracil_positions:?
# 20: custom_colors:xml
# 21: sequence_protein:utf-8
# 22:
# 23:
# 24:
# 25:
# 26:
# 27: unknown:binary
# 28: enzyme_vizualisation:xml
# 29: history_modifier:lzma
# 30: history_content:lzma (content from file it was taken, except for sequence)
# 31:
# 32: sequence_rna:ascii
```

## Install

This project uses `uv`, which automatically handles `venv` and packages.

```bash

# to run scripts do
uv run main.py # or history_analysis.py

```

## Roadmap

- [ ] Improve SGFF parsing, unify TLV strategy
- [ ] Understand whole file structure
- [ ] Correctly parse into readable form every block
- [ ] Parse XML into pure JSON format
- [ ] Create writer (if possible)
- [ ] Refine, refactor reader/writer
- [ ] Proper documentation and README cleanup

# Acknowledgments

This project would not have been possible without previous work done by
- Damien Goutte-Gattat, see his PDF on SGFF structure: https://incenp.org/dvlpt/docs/binary-sequence-formats/binary-sequence-formats.pdf
- Isaac Luo, for his version of SnapGene reader: https://github.com/IsaacLuo/SnapGeneFileReader

# License

Distributed under MIT licence, see `LICENSE` for more.
