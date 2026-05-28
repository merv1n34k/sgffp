from pathlib import Path
import sgffp
import glob

test_folder = Path(__file__).parent / "data"


class TestExamples:
    def parse_fasta(self, file_path: Path) -> dict[str, str]:
        key = ""
        out_dict = dict()
        for line in file_path.read_text().splitlines():
            if line.startswith(">"):
                key = line[1:].strip()
            else:
                out_dict[key] = line.strip()
        return out_dict

    def test_example_gibson_assembly(self):
        sgff = sgffp.SgffReader.from_file(test_folder / "gibson_assembly.dna")
        fasta_dict = self.parse_fasta(test_folder / "gibson_assembly.fasta")
        assert sgff.block(0)["sequence"] == fasta_dict["product"]
        root_node = sgff.history.tree.root

        linearized_vector_tree_node, fragment_tree_node = root_node.children

        linear_vector_node = sgff.history.get_node(linearized_vector_tree_node.id)
        assert (
            linear_vector_node.sequence.upper()
            == fasta_dict["linearized_vector"].upper()
        )

        fragment_node = sgff.history.get_node(fragment_tree_node.id)
        assert fragment_node.sequence.upper() == fasta_dict["fragment"].upper()

    def test_example_degenerate_examples(self):
        """Test a few files that have degenerate sequences, and the history is
        just a circularization, so the sequence at the root node should be the
        same as the sequence at the first node.
        """

        for file in glob.glob(str(test_folder / "degenerate_examples" / "*.dna")):
            sgff = sgffp.SgffReader.from_file(Path(file))
            assert sgff.blocks[0][0]["sequence"] == sgff.history.get_sequence_at(
                0
            ), f"Sequence mismatch for {file}"
