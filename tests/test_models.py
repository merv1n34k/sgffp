"""
Tests for SGFF data models
"""

import pytest
from sgffp.models import (
    SgffSequence,
    SgffFeature,
    SgffFeatureList,
    SgffSegment,
    SgffHistory,
    SgffHistoryNode,
    SgffHistoryNodeContent,
    SgffHistoryTree,
    SgffHistoryTreeNode,
    SgffHistoryOligo,
    SgffInputSummary,
    HistoryOperation,
    SgffPrimer,
    SgffPrimerList,
    SgffNotes,
    SgffProperties,
    SgffAlignment,
    SgffAlignmentList,
    SgffTrace,
    SgffTraceList,
    SgffTraceClip,
    SgffTraceSamples,
)


class TestSgffSequence:
    def test_empty_blocks(self):
        """Empty blocks return empty sequence"""
        seq = SgffSequence({})
        assert seq.value == ""
        assert seq.length == 0

    def test_load_from_block_0(self):
        """Load DNA sequence from block 0"""
        blocks = {
            0: [{"sequence": "ATCG", "topology": "circular", "strandedness": "double"}]
        }
        seq = SgffSequence(blocks)
        assert seq.value == "ATCG"
        assert seq.length == 4
        assert seq.topology == "circular"
        assert seq.is_circular
        assert seq.is_double_stranded

    def test_modify_sequence(self):
        """Modify sequence updates blocks"""
        blocks = {0: [{"sequence": "ATCG"}]}
        seq = SgffSequence(blocks)
        seq.value = "GGGG"
        assert blocks[0][0]["sequence"] == "GGGG"

    def test_modify_topology(self):
        """Modify topology updates blocks"""
        blocks = {0: [{"sequence": "ATCG", "topology": "linear"}]}
        seq = SgffSequence(blocks)
        seq.topology = "circular"
        assert blocks[0][0]["topology"] == "circular"


class TestSgffFeature:
    def test_from_dict(self):
        """Create feature from dict"""
        data = {
            "name": "GFP",
            "type": "CDS",
            "strand": "+",
            "segments": [{"range": "1-100", "color": "#00FF00"}],
            "qualifiers": {"note": "Green fluorescent protein"},
        }
        feature = SgffFeature.from_dict(data)
        assert feature.name == "GFP"
        assert feature.type == "CDS"
        assert feature.strand == "+"
        assert len(feature.segments) == 1
        assert feature.start == 0
        assert feature.end == 100

    def test_to_dict(self):
        """Convert feature back to dict"""
        feature = SgffFeature(
            name="Test",
            type="gene",
            strand="-",
            segments=[SgffSegment(start=0, end=50)],
        )
        data = feature.to_dict()
        assert data["name"] == "Test"
        assert data["type"] == "gene"
        assert data["strand"] == "-"


class TestSgffFeatureList:
    def test_empty_blocks(self):
        """Empty blocks return empty list"""
        fl = SgffFeatureList({})
        assert len(fl) == 0

    def test_load_features(self):
        """Load features from block 10"""
        blocks = {
            10: [
                {
                    "features": [
                        {"name": "A", "type": "gene", "segments": []},
                        {"name": "B", "type": "CDS", "segments": []},
                    ]
                }
            ]
        }
        fl = SgffFeatureList(blocks)
        assert len(fl) == 2
        assert fl[0].name == "A"
        assert fl[1].name == "B"

    def test_add_feature(self):
        """Add feature updates blocks"""
        blocks = {10: [{"features": []}]}
        fl = SgffFeatureList(blocks)
        fl.add(SgffFeature(name="New", type="gene"))
        assert len(fl) == 1
        assert len(blocks[10][0]["features"]) == 1

    def test_remove_feature(self):
        """Remove feature updates blocks"""
        blocks = {
            10: [{"features": [{"name": "A", "type": "gene", "segments": []}]}]
        }
        fl = SgffFeatureList(blocks)
        fl.remove(0)
        assert len(fl) == 0

    def test_find_by_name(self):
        """Find feature by name"""
        blocks = {
            10: [{"features": [{"name": "Target", "type": "gene", "segments": []}]}]
        }
        fl = SgffFeatureList(blocks)
        f = fl.find_by_name("Target")
        assert f is not None
        assert f.name == "Target"

    def test_find_by_type(self):
        """Find features by type"""
        blocks = {
            10: [
                {
                    "features": [
                        {"name": "A", "type": "CDS", "segments": []},
                        {"name": "B", "type": "gene", "segments": []},
                        {"name": "C", "type": "CDS", "segments": []},
                    ]
                }
            ]
        }
        fl = SgffFeatureList(blocks)
        cds = fl.find_by_type("CDS")
        assert len(cds) == 2


class TestSgffHistoryNode:
    def test_from_dict(self):
        """Create node from dict"""
        data = {
            "node_index": 5,
            "sequence": "ATCG",
            "sequence_type": 0,
            "length": 4,
            "node_info": {30: [{8: [{"test": "value"}]}]},
        }
        node = SgffHistoryNode.from_dict(data)
        assert node.index == 5
        assert node.sequence == "ATCG"
        assert node.sequence_type == 0
        assert node.length == 4
        assert node.content is not None
        assert node.content.has_properties
        assert node.properties.exists

    def test_to_dict(self):
        """Convert node back to dict"""
        node = SgffHistoryNode(index=3, sequence="GGG", sequence_type=1)
        data = node.to_dict()
        assert data["node_index"] == 3
        assert data["sequence"] == "GGG"
        assert data["sequence_type"] == 1


class TestSgffHistory:
    def test_empty_blocks(self):
        """Empty blocks return empty history"""
        h = SgffHistory({})
        assert len(h) == 0
        assert not h.exists

    def test_load_nodes(self):
        """Load history nodes from block 11"""
        blocks = {
            11: [
                {"node_index": 0, "sequence": "AAA", "sequence_type": 0},
                {"node_index": 1, "sequence": "BBB", "sequence_type": 0},
            ]
        }
        h = SgffHistory(blocks)
        assert len(h) == 2
        assert h.get_node(0).sequence == "AAA"
        assert h.get_node(1).sequence == "BBB"

    def test_get_sequence_at(self):
        """Get sequence at node index"""
        blocks = {11: [{"node_index": 2, "sequence": "TCGA", "sequence_type": 0}]}
        h = SgffHistory(blocks)
        assert h.get_sequence_at(2) == "TCGA"
        assert h.get_sequence_at(99) is None

    def test_add_node(self):
        """Add node updates blocks"""
        blocks = {}
        h = SgffHistory(blocks)
        h.add_node(SgffHistoryNode(index=0, sequence="NEW"))
        assert 11 in blocks
        assert len(h) == 1

    def test_remove_node(self):
        """Remove node updates blocks"""
        blocks = {11: [{"node_index": 0, "sequence": "DEL", "sequence_type": 0}]}
        h = SgffHistory(blocks)
        assert h.remove_node(0)
        assert len(h) == 0

    def test_update_node(self):
        """Update node attributes"""
        blocks = {11: [{"node_index": 0, "sequence": "OLD", "sequence_type": 0}]}
        h = SgffHistory(blocks)
        h.update_node(0, sequence="NEW")
        assert h.get_node(0).sequence == "NEW"

    def test_clear(self):
        """Clear removes all history blocks"""
        blocks = {7: [{}], 11: [{}], 29: [{}], 30: [{}]}
        h = SgffHistory(blocks)
        h.clear()
        assert 7 not in blocks
        assert 11 not in blocks
        assert 29 not in blocks
        assert 30 not in blocks

    def test_tree_node_linking(self):
        """Tree nodes are linked to sequence nodes"""
        blocks = {
            7: [{"HistoryTree": {"Node": {
                "ID": "1", "name": "current.dna", "type": "DNA", "seqLen": "100",
                "strandedness": "double", "circular": "0", "operation": "makeDna",
                "upstreamModification": "Unmodified", "downstreamModification": "Unmodified",
                "Node": {
                    "ID": "0", "name": "original.rna", "type": "RNA", "seqLen": "100",
                    "strandedness": "single", "circular": "0", "operation": "invalid",
                    "upstreamModification": "Unmodified", "downstreamModification": "Unmodified",
                    "resurrectable": "1"
                }
            }}}],
            11: [
                {"node_index": 0, "sequence": "RNA_SEQ", "sequence_type": 32},
                {"node_index": 1, "sequence": "DNA_SEQ", "sequence_type": 0},
            ]
        }
        h = SgffHistory(blocks)

        # Tree should be parsed
        assert h.tree is not None
        assert h.tree.root.name == "current.dna"

        # Nodes should be linked to tree
        node_0 = h.get_node(0)
        assert node_0.tree_node is not None
        assert node_0.tree_node.name == "original.rna"
        assert node_0.tree_node.type == "RNA"


class TestSgffHistoryOligo:
    def test_from_dict(self):
        """Create oligo from dict"""
        data = {"name": "FWD_Primer", "sequence": "ATGCATGC", "phosphorylated": "1"}
        oligo = SgffHistoryOligo.from_dict(data)
        assert oligo.name == "FWD_Primer"
        assert oligo.sequence == "ATGCATGC"
        assert oligo.phosphorylated is True

    def test_to_dict(self):
        """Convert oligo to dict"""
        oligo = SgffHistoryOligo(name="REV", sequence="GCTAGCTA", phosphorylated=False)
        data = oligo.to_dict()
        assert data["name"] == "REV"
        assert data["sequence"] == "GCTAGCTA"
        assert "phosphorylated" not in data  # Only included if True


class TestSgffInputSummary:
    def test_from_dict(self):
        """Create input summary from dict"""
        data = {"manipulation": "amplify", "val1": "10", "val2": "500"}
        summary = SgffInputSummary.from_dict(data)
        assert summary.manipulation == "amplify"
        assert summary.val1 == 10
        assert summary.val2 == 500

    def test_to_dict(self):
        """Convert input summary to dict"""
        summary = SgffInputSummary(manipulation="select", val1=0, val2=100)
        data = summary.to_dict()
        assert data["manipulation"] == "select"
        assert data["val1"] == "0"
        assert data["val2"] == "100"


class TestSgffHistoryTreeNode:
    def test_from_dict_simple(self):
        """Create tree node from simple dict"""
        data = {
            "ID": "5",
            "name": "test.dna",
            "type": "DNA",
            "seqLen": "1000",
            "strandedness": "double",
            "circular": "1",
            "operation": "amplifyFragment",
            "upstreamModification": "FivePrimePhosphorylated",
            "downstreamModification": "Unmodified",
            "resurrectable": "1",
        }
        node = SgffHistoryTreeNode.from_dict(data)
        assert node.id == 5
        assert node.name == "test.dna"
        assert node.type == "DNA"
        assert node.seq_len == 1000
        assert node.circular is True
        assert node.operation == "amplifyFragment"
        assert node.resurrectable is True

    def test_from_dict_with_oligos(self):
        """Parse oligos from tree node"""
        data = {
            "ID": "1", "name": "amp.dna", "type": "DNA", "seqLen": "100",
            "strandedness": "double", "circular": "0", "operation": "amplifyFragment",
            "upstreamModification": "Unmodified", "downstreamModification": "Unmodified",
            "Oligo": [
                {"name": "FWD", "sequence": "ATGC", "phosphorylated": "1"},
                {"name": "REV", "sequence": "GCTA"},
            ],
        }
        node = SgffHistoryTreeNode.from_dict(data)
        assert len(node.oligos) == 2
        assert node.oligos[0].name == "FWD"
        assert node.oligos[0].phosphorylated is True
        assert node.oligos[1].phosphorylated is False

    def test_from_dict_with_children(self):
        """Parse nested child nodes"""
        data = {
            "ID": "2", "name": "child2.dna", "type": "DNA", "seqLen": "200",
            "strandedness": "double", "circular": "0", "operation": "makeDna",
            "upstreamModification": "Unmodified", "downstreamModification": "Unmodified",
            "Node": {
                "ID": "1", "name": "child1.dna", "type": "DNA", "seqLen": "200",
                "strandedness": "double", "circular": "0", "operation": "invalid",
                "upstreamModification": "Unmodified", "downstreamModification": "Unmodified",
            }
        }
        node = SgffHistoryTreeNode.from_dict(data)
        assert len(node.children) == 1
        assert node.children[0].id == 1
        assert node.children[0].parent is node

    def test_to_dict(self):
        """Convert tree node back to dict"""
        node = SgffHistoryTreeNode(
            id=1, name="test.dna", type="DNA", seq_len=500,
            strandedness="double", circular=False, operation="makeDna",
            upstream_modification="Unmodified", downstream_modification="Unmodified",
            resurrectable=True,
            oligos=[SgffHistoryOligo(name="P1", sequence="AAAA")],
            input_summaries=[SgffInputSummary(manipulation="select", val1=0, val2=499)],
        )
        data = node.to_dict()
        assert data["ID"] == "1"
        assert data["name"] == "test.dna"
        assert data["resurrectable"] == "1"
        assert data["Oligo"]["name"] == "P1"
        assert data["InputSummary"]["manipulation"] == "select"


class TestSgffHistoryTree:
    def test_empty_tree(self):
        """Empty data creates empty tree"""
        tree = SgffHistoryTree(None)
        assert tree.root is None
        assert len(tree) == 0

    def test_parse_tree(self):
        """Parse tree from block 7 data"""
        data = {"HistoryTree": {"Node": {
            "ID": "2", "name": "current.dna", "type": "DNA", "seqLen": "100",
            "strandedness": "double", "circular": "0", "operation": "makeDna",
            "upstreamModification": "Unmodified", "downstreamModification": "Unmodified",
            "Node": {
                "ID": "1", "name": "previous.rna", "type": "RNA", "seqLen": "100",
                "strandedness": "single", "circular": "0", "operation": "invalid",
                "upstreamModification": "Unmodified", "downstreamModification": "Unmodified",
            }
        }}}
        tree = SgffHistoryTree(data)
        assert tree.root is not None
        assert tree.root.id == 2
        assert len(tree) == 2

    def test_get_node(self):
        """Get node by ID"""
        data = {"HistoryTree": {"Node": {
            "ID": "1", "name": "test.dna", "type": "DNA", "seqLen": "100",
            "strandedness": "double", "circular": "0", "operation": "invalid",
            "upstreamModification": "Unmodified", "downstreamModification": "Unmodified",
        }}}
        tree = SgffHistoryTree(data)
        assert tree.get(1).name == "test.dna"
        assert tree.get(999) is None

    def test_walk(self):
        """Walk tree depth-first"""
        data = {"HistoryTree": {"Node": {
            "ID": "2", "name": "node2", "type": "DNA", "seqLen": "100",
            "strandedness": "double", "circular": "0", "operation": "makeDna",
            "upstreamModification": "Unmodified", "downstreamModification": "Unmodified",
            "Node": {
                "ID": "1", "name": "node1", "type": "DNA", "seqLen": "100",
                "strandedness": "double", "circular": "0", "operation": "invalid",
                "upstreamModification": "Unmodified", "downstreamModification": "Unmodified",
            }
        }}}
        tree = SgffHistoryTree(data)
        names = [n.name for n in tree.walk()]
        assert names == ["node2", "node1"]

    def test_ancestors(self):
        """Get ancestor chain"""
        data = {"HistoryTree": {"Node": {
            "ID": "2", "name": "grandchild", "type": "DNA", "seqLen": "100",
            "strandedness": "double", "circular": "0", "operation": "makeDna",
            "upstreamModification": "Unmodified", "downstreamModification": "Unmodified",
            "Node": {
                "ID": "1", "name": "child", "type": "DNA", "seqLen": "100",
                "strandedness": "double", "circular": "0", "operation": "makeRna",
                "upstreamModification": "Unmodified", "downstreamModification": "Unmodified",
                "Node": {
                    "ID": "0", "name": "root", "type": "DNA", "seqLen": "100",
                    "strandedness": "double", "circular": "0", "operation": "invalid",
                    "upstreamModification": "Unmodified", "downstreamModification": "Unmodified",
                }
            }
        }}}
        tree = SgffHistoryTree(data)
        ancestors = tree.ancestors(0)
        assert len(ancestors) == 3
        assert ancestors[0].name == "root"
        assert ancestors[2].name == "grandchild"

    def test_to_dict(self):
        """Serialize tree back to dict"""
        data = {"HistoryTree": {"Node": {
            "ID": "1", "name": "test", "type": "DNA", "seqLen": "100",
            "strandedness": "double", "circular": "0", "operation": "invalid",
            "upstreamModification": "Unmodified", "downstreamModification": "Unmodified",
        }}}
        tree = SgffHistoryTree(data)
        result = tree.to_dict()
        assert "HistoryTree" in result
        assert result["HistoryTree"]["Node"]["name"] == "test"


class TestSgffHistoryNodeContent:
    def test_empty_content(self):
        """Empty data creates empty content"""
        content = SgffHistoryNodeContent.from_dict({})
        assert not content.exists
        assert not content.has_features
        assert not content.has_primers
        assert not content.has_notes
        assert len(content.features) == 0

    def test_parse_content(self):
        """Parse content from node_info"""
        data = {30: [{
            8: [{"AdditionalSequenceProperties": {"UpstreamStickiness": "0"}}],
            5: [{"Primers": {"nextValidID": "5"}}],
            6: [{"Notes": {"Type": "Synthetic"}}],
        }]}
        content = SgffHistoryNodeContent.from_dict(data)
        assert content.exists
        assert content.has_properties
        assert content.has_primers
        assert content.has_notes

    def test_to_dict(self):
        """Serialize content back to dict"""
        blocks = {
            8: [{"test": "value"}],
            5: [{"Primers": {}}],
        }
        content = SgffHistoryNodeContent(blocks)
        data = content.to_dict()
        assert 30 in data
        assert 8 in data[30][0]
        assert 5 in data[30][0]

    def test_block_access(self):
        """Access raw blocks"""
        blocks = {10: [{"features": []}], 18: [{"bases": "ATCG"}]}
        content = SgffHistoryNodeContent(blocks)
        assert content.block_types == [10, 18]
        assert 10 in content
        assert 18 in content
        assert 99 not in content

    def test_model_accessors(self):
        """Model accessors return appropriate types"""
        blocks = {10: [{"features": [{"name": "test", "type": "gene"}]}]}
        content = SgffHistoryNodeContent(blocks)
        assert len(content.features) == 1
        assert content.features[0].name == "test"

    def test_traces_accessor(self):
        """Multiple traces accessible via traces property"""
        blocks = {16: [
            {"flags": 0, "blocks": {18: [{"bases": "ATCG"}]}},
            {"flags": 0, "blocks": {18: [{"bases": "GGGG"}]}},
        ]}
        content = SgffHistoryNodeContent(blocks)
        assert content.has_traces
        assert len(content.traces) == 2
        assert content.traces[0].bases == "ATCG"
        assert content.traces[1].bases == "GGGG"


class TestHistoryOperation:
    def test_constants(self):
        """Operation constants are defined"""
        assert HistoryOperation.INVALID == "invalid"
        assert HistoryOperation.MAKE_DNA == "makeDna"
        assert HistoryOperation.AMPLIFY == "amplifyFragment"


class TestSgffPrimerList:
    def test_empty_blocks(self):
        """Empty blocks return empty list"""
        pl = SgffPrimerList({})
        assert len(pl) == 0

    def test_load_primers(self):
        """Load primers from block 5"""
        blocks = {
            5: [{"Primers": {"Primer": [{"name": "FWD", "sequence": "ATCG"}]}}]
        }
        pl = SgffPrimerList(blocks)
        assert len(pl) == 1
        assert pl[0].name == "FWD"
        assert pl[0].sequence == "ATCG"


class TestSgffNotes:
    def test_empty_blocks(self):
        """Empty blocks return empty notes"""
        n = SgffNotes({})
        assert n.description == ""

    def test_load_notes(self):
        """Load notes from block 6"""
        blocks = {6: [{"Notes": {"Description": "Test plasmid"}}]}
        n = SgffNotes(blocks)
        assert n.description == "Test plasmid"

    def test_set_note(self):
        """Set note updates blocks"""
        blocks = {6: [{"Notes": {}}]}
        n = SgffNotes(blocks)
        n.description = "Updated"
        assert blocks[6][0]["Notes"]["Description"] == "Updated"


class TestSgffProperties:
    def test_empty_blocks(self):
        """Empty blocks return empty properties"""
        p = SgffProperties({})
        assert not p.exists

    def test_load_properties(self):
        """Load properties from block 8"""
        blocks = {8: [{"key": "value"}]}
        p = SgffProperties(blocks)
        assert p.get("key") == "value"


class TestSgffAlignmentList:
    def test_empty_blocks(self):
        """Empty blocks return empty list"""
        al = SgffAlignmentList({})
        assert len(al) == 0

    def test_load_alignments(self):
        """Load alignments from block 17"""
        blocks = {
            17: [
                {"AlignableSequences": {"Sequence": [{"name": "Ref", "sequence": "ATCG"}]}}
            ]
        }
        al = SgffAlignmentList(blocks)
        assert len(al) == 1
        assert al[0].name == "Ref"


class TestSgffTraceClip:
    def test_from_dict(self):
        """Create clip from dict"""
        data = {"left": 10, "right": 500}
        clip = SgffTraceClip.from_dict(data)
        assert clip.left == 10
        assert clip.right == 500

    def test_to_dict(self):
        """Convert clip to dict"""
        clip = SgffTraceClip(left=5, right=100)
        data = clip.to_dict()
        assert data["left"] == 5
        assert data["right"] == 100

    def test_defaults(self):
        """Missing values default to 0"""
        clip = SgffTraceClip.from_dict({})
        assert clip.left == 0
        assert clip.right == 0


class TestSgffTraceSamples:
    def test_from_dict(self):
        """Create samples from dict"""
        data = {"A": [1, 2, 3], "C": [4, 5, 6], "G": [7, 8, 9], "T": [10, 11, 12]}
        samples = SgffTraceSamples.from_dict(data)
        assert samples.a == [1, 2, 3]
        assert samples.c == [4, 5, 6]
        assert samples.g == [7, 8, 9]
        assert samples.t == [10, 11, 12]

    def test_to_dict(self):
        """Convert samples to dict"""
        samples = SgffTraceSamples(a=[1, 2], c=[3, 4], g=[5, 6], t=[7, 8])
        data = samples.to_dict()
        assert data["A"] == [1, 2]
        assert data["C"] == [3, 4]
        assert data["G"] == [5, 6]
        assert data["T"] == [7, 8]

    def test_to_dict_empty_channels(self):
        """Empty channels are omitted from dict"""
        samples = SgffTraceSamples(a=[1, 2], c=[], g=[], t=[])
        data = samples.to_dict()
        assert data == {"A": [1, 2]}

    def test_length(self):
        """Length returns max channel length"""
        samples = SgffTraceSamples(a=[1, 2, 3], c=[1, 2], g=[], t=[1])
        assert samples.length == 3
        assert len(samples) == 3

    def test_empty_length(self):
        """Empty samples have length 0"""
        samples = SgffTraceSamples()
        assert len(samples) == 0


class TestSgffTrace:
    def test_from_dict(self):
        """Create trace from dict"""
        data = {
            "bases": "ATCGATCG",
            "positions": [10, 20, 30, 40, 50, 60, 70, 80],
            "confidence": [40, 45, 50, 55, 40, 45, 50, 55],
            "samples": {"A": [1, 2], "C": [3, 4], "G": [5, 6], "T": [7, 8]},
            "clip": {"left": 5, "right": 100},
            "text": {"MACH": "ABI3730"},
            "comments": ["Test trace"],
        }
        trace = SgffTrace.from_dict(data)
        assert trace.bases == "ATCGATCG"
        assert trace.sequence == "ATCGATCG"
        assert trace.length == 8
        assert len(trace) == 8
        assert trace.positions == [10, 20, 30, 40, 50, 60, 70, 80]
        assert trace.confidence == [40, 45, 50, 55, 40, 45, 50, 55]
        assert trace.samples.a == [1, 2]
        assert trace.sample_count == 2
        assert trace.clip.left == 5
        assert trace.clip.right == 100
        assert trace.text == {"MACH": "ABI3730"}
        assert trace.comments == ["Test trace"]

    def test_to_dict(self):
        """Convert trace to dict"""
        trace = SgffTrace(
            bases="ATCG",
            positions=[10, 20, 30, 40],
            confidence=[40, 45, 50, 55],
            samples=SgffTraceSamples(a=[1, 2], c=[3, 4], g=[5, 6], t=[7, 8]),
            clip=SgffTraceClip(left=5, right=100),
            text={"MACH": "ABI3730"},
            comments=["Test"],
        )
        data = trace.to_dict()
        assert data["bases"] == "ATCG"
        assert data["positions"] == [10, 20, 30, 40]
        assert data["samples"]["A"] == [1, 2]
        assert data["clip"]["left"] == 5

    def test_empty_trace(self):
        """Empty trace has defaults"""
        trace = SgffTrace()
        assert trace.bases == ""
        assert trace.positions == []
        assert trace.confidence == []
        assert trace.samples is None
        assert trace.clip is None
        assert trace.length == 0

    def test_get_confidence_at(self):
        """Get confidence at specific base"""
        trace = SgffTrace(bases="ATCG", confidence=[10, 20, 30, 40])
        assert trace.get_confidence_at(0) == 10
        assert trace.get_confidence_at(2) == 30
        assert trace.get_confidence_at(10) is None

    def test_get_position_at(self):
        """Get sample position at specific base"""
        trace = SgffTrace(bases="ATCG", positions=[100, 200, 300, 400])
        assert trace.get_position_at(0) == 100
        assert trace.get_position_at(3) == 400
        assert trace.get_position_at(10) is None

    def test_get_metadata(self):
        """Get metadata by key"""
        trace = SgffTrace(text={"MACH": "ABI3730", "LANE": "5"})
        assert trace.get_metadata("MACH") == "ABI3730"
        assert trace.get_metadata("LANE") == "5"
        assert trace.get_metadata("MISSING", "default") == "default"

    def test_repr(self):
        """String representation"""
        trace = SgffTrace(bases="ATCG", samples=SgffTraceSamples(a=[1, 2, 3]))
        assert "SgffTrace" in repr(trace)
        assert "bases=4" in repr(trace)
        assert "samples=3" in repr(trace)


class TestSgffTraceList:
    def test_empty_blocks(self):
        """Empty blocks return empty list"""
        traces = SgffTraceList({})
        assert len(traces) == 0

    def test_load_traces_from_containers(self):
        """Load traces from block 16 containers"""
        blocks = {16: [
            {"flags": 0, "blocks": {18: [{"bases": "ATCG"}]}},
            {"flags": 1, "blocks": {18: [{"bases": "GGGG"}]}},
        ]}
        traces = SgffTraceList(blocks)
        assert len(traces) == 2
        assert traces[0].bases == "ATCG"
        assert traces[1].bases == "GGGG"

    def test_load_from_single_container(self):
        """Load single trace from block 16 container"""
        blocks = {16: [{"flags": 0, "blocks": {18: [{"bases": "ATCG"}]}}]}
        traces = SgffTraceList(blocks)
        assert len(traces) == 1
        assert traces[0].bases == "ATCG"

    def test_add_trace(self):
        """Add trace syncs to block 16 containers"""
        blocks = {16: [{"flags": 0, "blocks": {18: [{"bases": "ATCG"}]}}]}
        traces = SgffTraceList(blocks)
        traces.add(SgffTrace(bases="GGGG"))
        assert len(traces) == 2
        assert len(blocks[16]) == 2
        assert blocks[16][1]["blocks"][18][0]["bases"] == "GGGG"

    def test_remove_trace(self):
        """Remove trace from list"""
        blocks = {16: [
            {"flags": 0, "blocks": {18: [{"bases": "ATCG"}]}},
            {"flags": 0, "blocks": {18: [{"bases": "GGGG"}]}},
        ]}
        traces = SgffTraceList(blocks)
        traces.remove(0)
        assert len(traces) == 1
        assert traces[0].bases == "GGGG"

    def test_clear_traces(self):
        """Clear all traces"""
        blocks = {16: [{"flags": 0, "blocks": {18: [{"bases": "ATCG"}]}}]}
        traces = SgffTraceList(blocks)
        traces.clear()
        assert len(traces) == 0
        assert 16 not in blocks

    def test_iterate(self):
        """Iterate over traces"""
        blocks = {16: [
            {"flags": 0, "blocks": {18: [{"bases": "A"}]}},
            {"flags": 0, "blocks": {18: [{"bases": "T"}]}},
            {"flags": 0, "blocks": {18: [{"bases": "C"}]}},
        ]}
        traces = SgffTraceList(blocks)
        bases = [t.bases for t in traces]
        assert bases == ["A", "T", "C"]

    def test_repr(self):
        """String representation"""
        blocks = {16: [
            {"flags": 0, "blocks": {18: [{"bases": "ATCG"}]}},
            {"flags": 0, "blocks": {18: [{"bases": "GGGG"}]}},
        ]}
        traces = SgffTraceList(blocks)
        assert "SgffTraceList" in repr(traces)
        assert "count=2" in repr(traces)
