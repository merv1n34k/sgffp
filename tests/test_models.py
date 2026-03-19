"""
Tests for SGFF data models
"""

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
    SgffBindingSite,
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
        blocks = {10: [{"features": [{"name": "A", "type": "gene", "segments": []}]}]}
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
            7: [
                {
                    "HistoryTree": {
                        "Node": {
                            "ID": "1",
                            "name": "current.dna",
                            "type": "DNA",
                            "seqLen": "100",
                            "strandedness": "double",
                            "circular": "0",
                            "operation": "makeDna",
                            "upstreamModification": "Unmodified",
                            "downstreamModification": "Unmodified",
                            "Node": {
                                "ID": "0",
                                "name": "original.rna",
                                "type": "RNA",
                                "seqLen": "100",
                                "strandedness": "single",
                                "circular": "0",
                                "operation": "invalid",
                                "upstreamModification": "Unmodified",
                                "downstreamModification": "Unmodified",
                                "resurrectable": "1",
                            },
                        }
                    }
                }
            ],
            11: [
                {"node_index": 0, "sequence": "RNA_SEQ", "sequence_type": 32},
                {"node_index": 1, "sequence": "DNA_SEQ", "sequence_type": 0},
            ],
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
            "ID": "1",
            "name": "amp.dna",
            "type": "DNA",
            "seqLen": "100",
            "strandedness": "double",
            "circular": "0",
            "operation": "amplifyFragment",
            "upstreamModification": "Unmodified",
            "downstreamModification": "Unmodified",
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
            "ID": "2",
            "name": "child2.dna",
            "type": "DNA",
            "seqLen": "200",
            "strandedness": "double",
            "circular": "0",
            "operation": "makeDna",
            "upstreamModification": "Unmodified",
            "downstreamModification": "Unmodified",
            "Node": {
                "ID": "1",
                "name": "child1.dna",
                "type": "DNA",
                "seqLen": "200",
                "strandedness": "double",
                "circular": "0",
                "operation": "invalid",
                "upstreamModification": "Unmodified",
                "downstreamModification": "Unmodified",
            },
        }
        node = SgffHistoryTreeNode.from_dict(data)
        assert len(node.children) == 1
        assert node.children[0].id == 1
        assert node.children[0].parent is node

    def test_to_dict(self):
        """Convert tree node back to dict"""
        node = SgffHistoryTreeNode(
            id=1,
            name="test.dna",
            type="DNA",
            seq_len=500,
            strandedness="double",
            circular=False,
            operation="makeDna",
            upstream_modification="Unmodified",
            downstream_modification="Unmodified",
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
        data = {
            "HistoryTree": {
                "Node": {
                    "ID": "2",
                    "name": "current.dna",
                    "type": "DNA",
                    "seqLen": "100",
                    "strandedness": "double",
                    "circular": "0",
                    "operation": "makeDna",
                    "upstreamModification": "Unmodified",
                    "downstreamModification": "Unmodified",
                    "Node": {
                        "ID": "1",
                        "name": "previous.rna",
                        "type": "RNA",
                        "seqLen": "100",
                        "strandedness": "single",
                        "circular": "0",
                        "operation": "invalid",
                        "upstreamModification": "Unmodified",
                        "downstreamModification": "Unmodified",
                    },
                }
            }
        }
        tree = SgffHistoryTree(data)
        assert tree.root is not None
        assert tree.root.id == 2
        assert len(tree) == 2

    def test_get_node(self):
        """Get node by ID"""
        data = {
            "HistoryTree": {
                "Node": {
                    "ID": "1",
                    "name": "test.dna",
                    "type": "DNA",
                    "seqLen": "100",
                    "strandedness": "double",
                    "circular": "0",
                    "operation": "invalid",
                    "upstreamModification": "Unmodified",
                    "downstreamModification": "Unmodified",
                }
            }
        }
        tree = SgffHistoryTree(data)
        assert tree.get(1).name == "test.dna"
        assert tree.get(999) is None

    def test_walk(self):
        """Walk tree depth-first"""
        data = {
            "HistoryTree": {
                "Node": {
                    "ID": "2",
                    "name": "node2",
                    "type": "DNA",
                    "seqLen": "100",
                    "strandedness": "double",
                    "circular": "0",
                    "operation": "makeDna",
                    "upstreamModification": "Unmodified",
                    "downstreamModification": "Unmodified",
                    "Node": {
                        "ID": "1",
                        "name": "node1",
                        "type": "DNA",
                        "seqLen": "100",
                        "strandedness": "double",
                        "circular": "0",
                        "operation": "invalid",
                        "upstreamModification": "Unmodified",
                        "downstreamModification": "Unmodified",
                    },
                }
            }
        }
        tree = SgffHistoryTree(data)
        names = [n.name for n in tree.walk()]
        assert names == ["node2", "node1"]

    def test_ancestors(self):
        """Get ancestor chain"""
        data = {
            "HistoryTree": {
                "Node": {
                    "ID": "2",
                    "name": "grandchild",
                    "type": "DNA",
                    "seqLen": "100",
                    "strandedness": "double",
                    "circular": "0",
                    "operation": "makeDna",
                    "upstreamModification": "Unmodified",
                    "downstreamModification": "Unmodified",
                    "Node": {
                        "ID": "1",
                        "name": "child",
                        "type": "DNA",
                        "seqLen": "100",
                        "strandedness": "double",
                        "circular": "0",
                        "operation": "makeRna",
                        "upstreamModification": "Unmodified",
                        "downstreamModification": "Unmodified",
                        "Node": {
                            "ID": "0",
                            "name": "root",
                            "type": "DNA",
                            "seqLen": "100",
                            "strandedness": "double",
                            "circular": "0",
                            "operation": "invalid",
                            "upstreamModification": "Unmodified",
                            "downstreamModification": "Unmodified",
                        },
                    },
                }
            }
        }
        tree = SgffHistoryTree(data)
        ancestors = tree.ancestors(0)
        assert len(ancestors) == 3
        assert ancestors[0].name == "root"
        assert ancestors[2].name == "grandchild"

    def test_to_dict(self):
        """Serialize tree back to dict"""
        data = {
            "HistoryTree": {
                "Node": {
                    "ID": "1",
                    "name": "test",
                    "type": "DNA",
                    "seqLen": "100",
                    "strandedness": "double",
                    "circular": "0",
                    "operation": "invalid",
                    "upstreamModification": "Unmodified",
                    "downstreamModification": "Unmodified",
                }
            }
        }
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
        data = {
            30: [
                {
                    8: [{"AdditionalSequenceProperties": {"UpstreamStickiness": "0"}}],
                    5: [{"Primers": {"nextValidID": "5"}}],
                    6: [{"Notes": {"Type": "Synthetic"}}],
                }
            ]
        }
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
        assert 8 in data
        assert 5 in data

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
        blocks = {
            16: [
                {"flags": 0, "blocks": {18: [{"bases": "ATCG"}]}},
                {"flags": 0, "blocks": {18: [{"bases": "GGGG"}]}},
            ]
        }
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
        assert HistoryOperation.INSERT_MULTI == "insertFragments"
        assert HistoryOperation.FLIP == "flip"
        assert HistoryOperation.MUTAGENESIS == "primerDirectedMutagenesis"
        assert HistoryOperation.CHANGE_METHYLATION == "changeMethylation"

    def test_unknown_operation_preserved(self):
        """Unknown operation values pass through instead of becoming INVALID"""
        op = HistoryOperation("customOperation")
        assert op == "customOperation"
        assert op.value == "customOperation"

    def test_categories(self):
        """All 26 known operations are defined"""
        # Count unique values
        values = {m.value for m in HistoryOperation}
        assert len(values) >= 25


class TestSgffPrimerList:
    def test_empty_blocks(self):
        """Empty blocks return empty list"""
        pl = SgffPrimerList({})
        assert len(pl) == 0

    def test_load_primers(self):
        """Load primers from block 5"""
        blocks = {5: [{"Primers": {"Primer": [{"name": "FWD", "sequence": "ATCG"}]}}]}
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
                {
                    "AlignableSequences": {
                        "Sequence": [{"name": "Ref", "sequence": "ATCG"}]
                    }
                }
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
        blocks = {
            16: [
                {"flags": 0, "blocks": {18: [{"bases": "ATCG"}]}},
                {"flags": 1, "blocks": {18: [{"bases": "GGGG"}]}},
            ]
        }
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
        blocks = {
            16: [
                {"flags": 0, "blocks": {18: [{"bases": "ATCG"}]}},
                {"flags": 0, "blocks": {18: [{"bases": "GGGG"}]}},
            ]
        }
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
        blocks = {
            16: [
                {"flags": 0, "blocks": {18: [{"bases": "A"}]}},
                {"flags": 0, "blocks": {18: [{"bases": "T"}]}},
                {"flags": 0, "blocks": {18: [{"bases": "C"}]}},
            ]
        }
        traces = SgffTraceList(blocks)
        bases = [t.bases for t in traces]
        assert bases == ["A", "T", "C"]

    def test_repr(self):
        """String representation"""
        blocks = {
            16: [
                {"flags": 0, "blocks": {18: [{"bases": "ATCG"}]}},
                {"flags": 0, "blocks": {18: [{"bases": "GGGG"}]}},
            ]
        }
        traces = SgffTraceList(blocks)
        assert "SgffTraceList" in repr(traces)
        assert "count=2" in repr(traces)


class TestSgffBindingSite:
    def test_from_dict_basic(self):
        """Parse a basic BindingSite dict"""
        data = {
            "location": "2252-2290",
            "boundStrand": "1",
            "annealedBases": "GTTGGGGTCTTTGCTCAG",
            "meltingTemperature": "71",
            "Component": {"hybridizedRange": "2252-2290", "bases": "GTTGGG"},
        }
        bs = SgffBindingSite.from_dict(data)
        assert bs.start == 2251  # 0-based
        assert bs.end == 2290
        assert bs.bound_strand == "-"
        assert bs.annealed_bases == "GTTGGGGTCTTTGCTCAG"
        assert bs.melting_temperature == 71.0
        assert bs.simplified is False
        assert "Component" in bs.extras

    def test_from_dict_simplified(self):
        """Simplified flag is parsed"""
        data = {
            "simplified": "1",
            "location": "100-200",
            "boundStrand": "0",
        }
        bs = SgffBindingSite.from_dict(data)
        assert bs.simplified is True
        assert bs.bound_strand == "+"

    def test_roundtrip(self):
        """to_dict reverses from_dict"""
        data = {
            "location": "2252-2290",
            "boundStrand": "1",
            "annealedBases": "GTTGGG",
            "meltingTemperature": "71",
            "Component": {"hybridizedRange": "2252-2290", "bases": "GTTGGG"},
        }
        bs = SgffBindingSite.from_dict(data)
        result = bs.to_dict()
        assert result["location"] == "2252-2290"
        assert result["boundStrand"] == "1"
        assert result["annealedBases"] == "GTTGGG"
        assert result["meltingTemperature"] == "71"
        assert result["Component"]["hybridizedRange"] == "2252-2290"

    def test_simplified_roundtrip(self):
        """Simplified flag survives roundtrip"""
        data = {"simplified": "1", "location": "100-200", "boundStrand": "0"}
        bs = SgffBindingSite.from_dict(data)
        result = bs.to_dict()
        assert result["simplified"] == "1"

    def test_forward_strand(self):
        """boundStrand 0 maps to +"""
        bs = SgffBindingSite.from_dict({"location": "1-10", "boundStrand": "0"})
        assert bs.bound_strand == "+"
        assert bs.to_dict()["boundStrand"] == "0"


class TestSgffPrimerExtras:
    def test_primer_extras_preserved(self):
        """Unknown attrs survive roundtrip; BindingSite is now parsed, not in extras"""
        data = {
            "name": "FWD",
            "sequence": "ATCG",
            "recentID": "5",
            "description": "Forward primer",
            "color": "#FF0000",
            "dateAdded": "2024-01-01",
            "BindingSite": {"location": "1-20", "boundStrand": "0"},
        }
        primer = SgffPrimer.from_dict(data)
        assert primer.name == "FWD"
        assert primer.extras["recentID"] == "5"
        assert primer.extras["color"] == "#FF0000"
        # BindingSite is now parsed into binding_sites, NOT in extras
        assert "BindingSite" not in primer.extras
        assert len(primer.binding_sites) == 1
        assert primer.binding_sites[0].start == 0
        assert primer.binding_sites[0].end == 20

        result = primer.to_dict()
        assert result["name"] == "FWD"
        assert result["recentID"] == "5"
        assert result["description"] == "Forward primer"
        assert result["color"] == "#FF0000"
        assert result["BindingSite"]["location"] == "1-20"

    def test_primer_clear_binding_sites(self):
        """clear_binding_sites clears the binding_sites list"""
        data = {
            "name": "FWD",
            "sequence": "ATCG",
            "BindingSite": {"location": "1-20", "boundStrand": "0"},
        }
        primer = SgffPrimer.from_dict(data)
        assert len(primer.binding_sites) == 1
        primer.clear_binding_sites()
        assert len(primer.binding_sites) == 0
        assert "BindingSite" not in primer.to_dict()

    def test_primer_computed_properties(self):
        """bind_position, bind_strand, melting_temperature from first non-simplified site"""
        data = {
            "name": "P1",
            "sequence": "ATCG",
            "BindingSite": [
                {
                    "location": "100-120",
                    "boundStrand": "1",
                    "meltingTemperature": "65",
                },
                {
                    "simplified": "1",
                    "location": "100-120",
                    "boundStrand": "1",
                    "meltingTemperature": "65",
                },
            ],
        }
        primer = SgffPrimer.from_dict(data)
        assert primer.bind_position == 99  # 0-based
        assert primer.bind_strand == "-"
        assert primer.melting_temperature == 65.0

    def test_primer_no_binding_sites_defaults(self):
        """Without binding sites, computed properties return defaults"""
        primer = SgffPrimer(name="P1", sequence="ATCG")
        assert primer.bind_position is None
        assert primer.bind_strand == "+"
        assert primer.melting_temperature is None

    def test_primer_multiple_sites_roundtrip(self):
        """Multiple BindingSite elements survive roundtrip"""
        data = {
            "name": "P1",
            "sequence": "ATCG",
            "BindingSite": [
                {"location": "100-120", "boundStrand": "0"},
                {"location": "500-520", "boundStrand": "1", "simplified": "1"},
            ],
        }
        primer = SgffPrimer.from_dict(data)
        assert len(primer.binding_sites) == 2

        result = primer.to_dict()
        assert isinstance(result["BindingSite"], list)
        assert len(result["BindingSite"]) == 2
        assert result["BindingSite"][0]["location"] == "100-120"
        assert result["BindingSite"][1]["location"] == "500-520"

    def test_primer_list_wrapper_extras(self):
        """Wrapper-level extras (HybridizationParams, nextValidID) preserved"""
        blocks = {
            5: [
                {
                    "Primers": {
                        "HybridizationParams": {"minContinuous": "10"},
                        "nextValidID": "3",
                        "Primer": [{"name": "P1", "sequence": "ATCG"}],
                    }
                }
            ]
        }
        pl = SgffPrimerList(blocks)
        assert len(pl) == 1
        assert pl._wrapper_extras["HybridizationParams"]["minContinuous"] == "10"
        assert pl._wrapper_extras["nextValidID"] == "3"

        # Sync and verify wrapper extras are preserved in block
        pl.add(SgffPrimer(name="P2", sequence="GGGG"))
        block_data = blocks[5][0]["Primers"]
        assert block_data["HybridizationParams"]["minContinuous"] == "10"
        assert block_data["nextValidID"] == "3"


class TestSgffAlignmentExtras:
    def test_alignment_extras_preserved(self):
        """Unknown attrs survive roundtrip via to_dict()"""
        data = {
            "name": "Ref",
            "sequence": "ATCG",
            "ID": "1",
            "sortOrder": "0",
            "trimmedRange": "5-100",
            "isTrace": "0",
        }
        alignment = SgffAlignment.from_dict(data)
        assert alignment.name == "Ref"
        assert alignment.extras["ID"] == "1"
        assert alignment.extras["sortOrder"] == "0"
        assert alignment.extras["trimmedRange"] == "5-100"

        result = alignment.to_dict()
        assert result["name"] == "Ref"
        assert result["ID"] == "1"
        assert result["trimmedRange"] == "5-100"

    def test_alignment_list_wrapper_extras(self):
        """Wrapper-level extras (trimStringency) preserved"""
        blocks = {
            17: [
                {
                    "AlignableSequences": {
                        "trimStringency": "0.05",
                        "Sequence": [{"name": "Ref", "sequence": "ATCG"}],
                    }
                }
            ]
        }
        al = SgffAlignmentList(blocks)
        assert len(al) == 1
        assert al._wrapper_extras["trimStringency"] == "0.05"

        # Sync and verify wrapper extras are preserved
        al.add(SgffAlignment(name="New", sequence="GGGG"))
        block_data = blocks[17][0]["AlignableSequences"]
        assert block_data["trimStringency"] == "0.05"


class TestSgffHistoryTreeNodeExtras:
    def test_tree_node_extras_preserved(self):
        """Unmodeled keys land in extras and survive roundtrip"""
        data = {
            "ID": "3",
            "name": "test.dna",
            "type": "DNA",
            "seqLen": "500",
            "strandedness": "double",
            "circular": "0",
            "operation": "makeDna",
            "upstreamModification": "Unmodified",
            "downstreamModification": "Unmodified",
            "hideHistory": "1",
            "customMapLabel": "My Label",
            "useCustomMapLabel": "1",
            "upstreamStickiness": "AATT",
            "downstreamStickiness": "TTAA",
            "RegeneratedSite": {"name": "BamHI"},
        }
        node = SgffHistoryTreeNode.from_dict(data)
        assert node.extras["hideHistory"] == "1"
        assert node.extras["customMapLabel"] == "My Label"
        assert node.extras["RegeneratedSite"]["name"] == "BamHI"

        result = node.to_dict()
        assert result["hideHistory"] == "1"
        assert result["customMapLabel"] == "My Label"
        assert result["upstreamStickiness"] == "AATT"
        assert result["RegeneratedSite"]["name"] == "BamHI"
        # Modeled fields still correct
        assert result["ID"] == "3"
        assert result["name"] == "test.dna"

    def test_tree_node_unmodified_omitted(self):
        """upstreamModification/downstreamModification omitted when Unmodified"""
        node = SgffHistoryTreeNode(
            id=1, name="t", type="DNA", seq_len=10,
            strandedness="double", circular=False,
            operation="makeDna",
            upstream_modification="Unmodified",
            downstream_modification="Unmodified",
        )
        result = node.to_dict()
        assert "upstreamModification" not in result
        assert "downstreamModification" not in result

    def test_tree_node_modification_emitted(self):
        """Non-default modifications are emitted"""
        node = SgffHistoryTreeNode(
            id=1, name="t", type="DNA", seq_len=10,
            strandedness="double", circular=False,
            operation="makeDna",
            upstream_modification="FivePrimePhosphorylated",
            downstream_modification="Unmodified",
        )
        result = node.to_dict()
        assert result["upstreamModification"] == "FivePrimePhosphorylated"
        assert "downstreamModification" not in result


class TestSgffInputSummaryExtras:
    def test_input_summary_extras_preserved(self):
        """Unmodeled keys (strainName, methylationChanges) survive roundtrip"""
        data = {
            "manipulation": "select",
            "val1": "10",
            "val2": "500",
            "strainName": "DH5alpha",
            "methylationChanges": "dam+",
            "dam2": "1",
        }
        summary = SgffInputSummary.from_dict(data)
        assert summary.extras["strainName"] == "DH5alpha"
        assert summary.extras["methylationChanges"] == "dam+"
        assert summary.extras["dam2"] == "1"

        result = summary.to_dict()
        assert result["manipulation"] == "select"
        assert result["strainName"] == "DH5alpha"
        assert result["methylationChanges"] == "dam+"
        assert result["dam2"] == "1"

    def test_input_summary_enzyme_keys_not_in_extras(self):
        """Enzyme name/siteCount keys are consumed, not in extras"""
        data = {
            "manipulation": "digest",
            "val1": "0",
            "val2": "100",
            "name1": "BamHI",
            "siteCount1": "2",
        }
        summary = SgffInputSummary.from_dict(data)
        assert len(summary.enzymes) == 1
        assert "name1" not in summary.extras
        assert "siteCount1" not in summary.extras


class TestSgffFeatureExtras:
    def test_feature_extras_preserved(self):
        """Unmodeled feature attrs survive roundtrip"""
        data = {
            "name": "GFP",
            "type": "CDS",
            "strand": "+",
            "segments": [{"range": "1-100", "color": "#00FF00"}],
            "qualifiers": {"note": "test"},
            "color": "#00FF00",
            "extras": {
                "recentID": "5",
                "readingFrame": "1",
                "translationMW": "26800",
                "hitsStopCodon": "1",
            },
        }
        feature = SgffFeature.from_dict(data)
        assert feature.extras["recentID"] == "5"
        assert feature.reading_frame == 1
        assert "readingFrame" not in feature.extras

        result = feature.to_dict()
        assert result["extras"]["recentID"] == "5"
        assert result["extras"]["translationMW"] == "26800"
        assert result["extras"]["readingFrame"] == "1"

    def test_reading_frame_none_by_default(self):
        """reading_frame is None when not present"""
        feature = SgffFeature.from_dict({"name": "X", "type": "gene"})
        assert feature.reading_frame is None

        result = feature.to_dict()
        assert "readingFrame" not in result["extras"]

    def test_reading_frame_negative_roundtrip(self):
        """Negative reading frame roundtrips correctly"""
        data = {
            "name": "CDS1",
            "type": "CDS",
            "strand": "-",
            "extras": {"readingFrame": "-1"},
        }
        feature = SgffFeature.from_dict(data)
        assert feature.reading_frame == -1

        result = feature.to_dict()
        assert result["extras"]["readingFrame"] == "-1"

    def test_feature_raw_qualifiers_preserved(self):
        """raw_qualifiers survive roundtrip"""
        raw_quals = [
            {"name": "note", "V": {"text": "A note"}},
            {"name": "product", "V": {"text": "GFP"}},
        ]
        data = {
            "name": "GFP",
            "type": "CDS",
            "segments": [],
            "qualifiers": {"note": "A note", "product": "GFP"},
            "raw_qualifiers": raw_quals,
        }
        feature = SgffFeature.from_dict(data)
        assert feature.raw_qualifiers is not None
        assert len(feature.raw_qualifiers) == 2

        result = feature.to_dict()
        assert result["raw_qualifiers"] == raw_quals


class TestSgffSegmentExtras:
    def test_segment_extras_preserved(self):
        """Unmodeled segment attrs survive roundtrip"""
        data = {
            "range": "1-100",
            "color": "#FF0000",
            "type": "standard",
            "name": "Seg1",
            "translated": "1",
        }
        segment = SgffSegment.from_dict(data)
        assert segment.type == "standard"
        assert segment.translated is True
        assert segment.extras["name"] == "Seg1"
        assert "type" not in segment.extras
        assert "translated" not in segment.extras

        result = segment.to_dict()
        assert result["type"] == "standard"
        assert result["translated"] == "1"
        assert result["name"] == "Seg1"
        assert result["range"] == "1-100"
        assert result["color"] == "#FF0000"

    def test_segment_defaults(self):
        """Default type is standard, translated is False"""
        segment = SgffSegment.from_dict({"range": "1-50"})
        assert segment.type == "standard"
        assert segment.translated is False

        result = segment.to_dict()
        assert result["type"] == "standard"
        assert "translated" not in result

    def test_segment_gap_type(self):
        """Gap segment type roundtrips"""
        segment = SgffSegment.from_dict({"range": "1-10", "type": "gap"})
        assert segment.type == "gap"
        assert segment.translated is False

        result = segment.to_dict()
        assert result["type"] == "gap"
        assert "translated" not in result

    def test_segment_translated_roundtrip(self):
        """Translated flag roundtrips as '1'"""
        segment = SgffSegment.from_dict({
            "range": "1-100",
            "type": "standard",
            "translated": "1",
        })
        assert segment.translated is True
        assert segment.type == "standard"

        result = segment.to_dict()
        assert result["translated"] == "1"
        assert result["type"] == "standard"


class TestSgffFeatureListExtras:
    def test_wrapper_extras_preserved(self):
        """Wrapper-level extras (nextValidID, recycledIDs) preserved"""
        blocks = {
            10: [
                {
                    "features": [
                        {"name": "A", "type": "gene", "segments": []},
                    ],
                    "wrapper_extras": {
                        "nextValidID": "10",
                        "recycledIDs": "3,5",
                    },
                }
            ]
        }
        fl = SgffFeatureList(blocks)
        assert len(fl) == 1
        assert fl._wrapper_extras["nextValidID"] == "10"
        assert fl._wrapper_extras["recycledIDs"] == "3,5"

        # Sync preserves wrapper extras
        fl.add(SgffFeature(name="B", type="CDS"))
        block_data = blocks[10][0]
        assert block_data["wrapper_extras"]["nextValidID"] == "10"


class TestSgffHistorySnapshot:
    def test_snapshot_creates_node(self):
        """Blocks with seq + features → correct snapshot node"""
        blocks = {
            0: [{"sequence": "ATCGATCG", "topology": "linear", "strandedness": "double"}],
            7: [{"HistoryTree": {"Node": {
                "ID": "1", "name": "test.dna", "type": "DNA",
                "seqLen": "8", "strandedness": "double",
                "circular": "0", "operation": "invalid",
            }}}],
            10: [{"features": [{"name": "GFP", "type": "CDS", "segments": []}]}],
        }
        history = SgffHistory(blocks)
        node = history.snapshot_current_state(blocks)

        assert node.index == 1
        assert node.sequence == "ATCGATCG"
        # Block 0 DNA is stored as compressed DNA (seq_type=1) in snapshots
        assert node.sequence_type == 1
        assert node.content is not None
        assert node.content.has_features

    def test_snapshot_compressed_dna(self):
        """Compressed DNA metadata preserved in snapshot"""
        blocks = {
            1: [{
                "sequence": "GATCGATC",
                "length": 8,
                "format_version": 30,
                "strandedness_flag": 1,
                "property_flags": 3,
                "header_seq_length": 8,
            }],
            7: [{"HistoryTree": {"Node": {
                "ID": "5", "name": "comp.dna", "type": "DNA",
                "seqLen": "8", "strandedness": "double",
                "circular": "0", "operation": "invalid",
            }}}],
        }
        history = SgffHistory(blocks)
        node = history.snapshot_current_state(blocks)

        assert node.sequence_type == 1
        assert node.format_version == 30
        assert node.property_flags == 3
        assert node.header_seq_length == 8
        assert node.index == 5

    def test_snapshot_empty_content(self):
        """Only sequence, no content blocks"""
        blocks = {
            0: [{"sequence": "ATCG", "topology": "linear", "strandedness": "single"}],
            7: [{"HistoryTree": {"Node": {
                "ID": "1", "name": "t", "type": "DNA",
                "seqLen": "4", "strandedness": "single",
                "circular": "0", "operation": "invalid",
            }}}],
        }
        history = SgffHistory(blocks)
        node = history.snapshot_current_state(blocks)

        assert node.sequence == "ATCG"
        assert node.content is None or not node.content.exists

    def test_next_id_with_tree(self):
        """next_id returns max(IDs) + 1"""
        blocks = {
            7: [{"HistoryTree": {"Node": {
                "ID": "3", "name": "root", "type": "DNA",
                "seqLen": "10", "strandedness": "double",
                "circular": "0", "operation": "invalid",
                "Node": {
                    "ID": "1", "name": "child", "type": "DNA",
                    "seqLen": "5", "strandedness": "double",
                    "circular": "0", "operation": "invalid",
                },
            }}}],
        }
        history = SgffHistory(blocks)
        assert history.next_id() == 4

    def test_next_id_no_tree(self):
        """next_id returns 1 when no tree exists"""
        history = SgffHistory({})
        assert history.next_id() == 1


class TestSgffHistoryRecordOperation:
    def _make_blocks(self, seq="ATCGATCG"):
        """Helper to create blocks with sequence, tree, and features."""
        return {
            0: [{"sequence": seq, "topology": "circular", "strandedness": "double"}],
            7: [{"HistoryTree": {"Node": {
                "ID": "1", "name": "test.dna", "type": "DNA",
                "seqLen": str(len(seq)), "strandedness": "double",
                "circular": "1", "operation": "invalid",
            }}}],
            10: [{"features": [{"name": "GFP", "type": "CDS", "segments": []}]}],
        }

    def test_record_operation_basic(self):
        """New root, old root demoted, block 11 created"""
        blocks = self._make_blocks()
        history = SgffHistory(blocks)

        assert len(history.tree) == 1
        assert len(history.nodes) == 0

        new_root = history.record_operation(
            blocks, "ATCGATCGGG", "insertFragment"
        )

        assert new_root is not None
        assert new_root.id == 2
        assert new_root.seq_len == 10
        assert new_root.operation == "insertFragment"
        assert new_root.circular is True
        assert len(new_root.children) == 1
        assert new_root.children[0].id == 1
        assert new_root.children[0].resurrectable is True
        assert len(history.tree) == 2
        assert len(history.nodes) == 1
        assert history.nodes[1].sequence == "ATCGATCG"

    def test_record_operation_preserves_content(self):
        """Snapshot captures features from content blocks"""
        blocks = self._make_blocks()
        history = SgffHistory(blocks)
        history.record_operation(blocks, "NEWSEQ", "replace")

        snapshot = history.nodes[1]
        assert snapshot.content is not None
        assert snapshot.content.has_features

    def test_record_operation_roundtrip(self):
        """write → read → verify full structure"""
        from sgffp.reader import SgffReader
        from sgffp.writer import SgffWriter
        from sgffp.internal import SgffObject, Cookie

        blocks = self._make_blocks("ATCGATCGATCG")
        sgff = SgffObject(cookie=Cookie(), blocks=blocks)

        sgff.history.record_operation(blocks, "ATCGATCGATCGAAA", "insertFragment")

        written = SgffWriter.to_bytes(sgff)
        restored = SgffReader.from_bytes(written)

        assert len(restored.history.tree) == 2
        assert len(restored.history.nodes) == 1
        assert restored.history.tree.root.seq_len == 15
        assert restored.history.tree.root.operation == "insertFragment"

    def test_record_operation_no_tree(self):
        """Graceful no-op when no history exists"""
        blocks = {0: [{"sequence": "ATCG"}]}
        history = SgffHistory(blocks)
        result = history.record_operation(blocks, "ATCGGG", "replace")
        assert result is None


class TestSgffObjectSetSequence:
    def test_set_sequence_records_history(self):
        """set_sequence creates a new history tree node"""
        from sgffp.internal import SgffObject, Cookie

        blocks = {
            0: [{"sequence": "ATCGATCG", "topology": "linear", "strandedness": "double"}],
            7: [{"HistoryTree": {"Node": {
                "ID": "1", "name": "test.dna", "type": "DNA",
                "seqLen": "8", "strandedness": "double",
                "circular": "0", "operation": "invalid",
            }}}],
        }
        sgff = SgffObject(cookie=Cookie(), blocks=blocks)
        old_tree_count = len(sgff.history.tree)

        sgff.set_sequence("ATCGATCGGG", operation="insertFragment")

        assert len(sgff.history.tree) == old_tree_count + 1
        assert sgff.sequence.value == "ATCGATCGGG"
        assert sgff.history.tree.root.seq_len == 10

    def test_set_sequence_no_history(self):
        """set_sequence works when no history exists"""
        from sgffp.internal import SgffObject

        sgff = SgffObject.new("ATCG")
        sgff.set_sequence("GGGG")
        assert sgff.sequence.value == "GGGG"


class TestHistoryOperationCategories:
    """Test all operation categories produce valid, roundtrippable history."""

    def _make_sgff(self, seq="ATCGATCGATCG", circular=True):
        """Create a minimal SgffObject with history."""
        from sgffp.internal import SgffObject, Cookie

        blocks = {
            0: [{"sequence": seq, "topology": "circular" if circular else "linear",
                 "strandedness": "double"}],
            7: [{"HistoryTree": {"Node": {
                "ID": "1", "name": "test.dna", "type": "DNA",
                "seqLen": str(len(seq)), "strandedness": "double",
                "circular": "1" if circular else "0",
                "operation": "invalid",
            }}}],
        }
        return SgffObject(cookie=Cookie(), blocks=blocks)

    def _roundtrip(self, sgff):
        """Write and re-read, return restored object."""
        from sgffp.reader import SgffReader
        from sgffp.writer import SgffWriter

        written = SgffWriter.to_bytes(sgff)
        return SgffReader.from_bytes(written)

    # --- Metadata operations (no sequence change) ---

    def test_change_methylation(self):
        """changeMethylation preserves sequence, grows tree by 1"""
        sgff = self._make_sgff()
        seq = sgff.sequence.value

        sgff.history.record_operation(
            sgff.blocks, seq, "changeMethylation",
            name="test.dna",
            InputSummary={"manipulation": "transformInto",
                          "strainName": "DH5alpha",
                          "methylationChanges": "Dam+,Dcm+"},
        )

        assert len(sgff.history.tree) == 2
        root = sgff.history.tree.root
        assert root.operation == "changeMethylation"
        assert root.seq_len == len(seq)

        restored = self._roundtrip(sgff)
        assert len(restored.history.tree) == 2
        assert restored.history.tree.root.operation == "changeMethylation"

    def test_change_topology(self):
        """changeTopology: linear → circular"""
        sgff = self._make_sgff(circular=False)
        seq = sgff.sequence.value

        sgff.history.record_operation(
            sgff.blocks, seq, "changeTopology",
            name="test.dna", circular=True,
        )

        root = sgff.history.tree.root
        assert root.operation == "changeTopology"
        assert root.circular is True
        assert root.children[0].circular is False

        restored = self._roundtrip(sgff)
        assert restored.history.tree.root.circular is True

    def test_change_strandedness(self):
        """changeStrandedness: double → single"""
        sgff = self._make_sgff()
        seq = sgff.sequence.value

        sgff.history.record_operation(
            sgff.blocks, seq, "changeStrandedness",
            name="test.dna", strandedness="single",
        )

        root = sgff.history.tree.root
        assert root.strandedness == "single"
        assert root.children[0].strandedness == "double"

    # --- Edit operations (sequence changes) ---

    def test_replace(self):
        """replace: point mutation"""
        sgff = self._make_sgff("ATCGATCGATCG")
        new_seq = "ATCGTTCGATCG"  # A→T at position 4

        sgff.set_sequence(new_seq, operation="replace")

        assert sgff.sequence.value == new_seq
        assert len(sgff.history.tree) == 2
        root = sgff.history.tree.root
        assert root.operation == "replace"
        assert root.seq_len == 12

        # Snapshot has old sequence
        assert sgff.history.nodes[1].sequence == "ATCGATCGATCG"

        restored = self._roundtrip(sgff)
        assert restored.sequence.value == new_seq
        assert len(restored.history.tree) == 2

    def test_insert_fragment(self):
        """insertFragment: insert 6bp BamHI site"""
        sgff = self._make_sgff("ATCGATCGATCG")
        new_seq = "ATCGGGATCCATCGATCG"  # Inserted GGATCC at pos 4

        sgff.set_sequence(new_seq, operation="insertFragment")

        root = sgff.history.tree.root
        assert root.operation == "insertFragment"
        assert root.seq_len == 18
        assert root.children[0].seq_len == 12

        restored = self._roundtrip(sgff)
        assert restored.sequence.value == new_seq

    def test_flip(self):
        """flip: reverse complement"""
        sgff = self._make_sgff("ATCGATCG")
        # Reverse complement of ATCGATCG = CGATCGAT
        new_seq = "CGATCGAT"

        sgff.history.record_operation(
            sgff.blocks, new_seq, "flip", name="test.dna",
        )
        sgff.sequence.value = new_seq

        root = sgff.history.tree.root
        assert root.operation == "flip"
        assert root.seq_len == 8

        restored = self._roundtrip(sgff)
        assert restored.sequence.value == new_seq

    def test_new_from_selection(self):
        """newFileFromSelection: extract subsequence"""
        sgff = self._make_sgff("ATCGATCGATCG", circular=False)
        new_seq = "GATCGA"  # positions 3-8

        sgff.history.record_operation(
            sgff.blocks, new_seq, "newFileFromSelection",
            name="test.dna",
            InputSummary={"manipulation": "select", "val1": "3", "val2": "8"},
        )
        sgff.sequence.value = new_seq

        root = sgff.history.tree.root
        assert root.operation == "newFileFromSelection"
        assert root.seq_len == 6

        restored = self._roundtrip(sgff)
        assert restored.sequence.value == new_seq

    def test_mutagenesis(self):
        """primerDirectedMutagenesis: SDM with primer"""
        sgff = self._make_sgff("ATCGATCGATCG")
        new_seq = "ATCGAAAGGATCG"  # TC→AAAG at pos 4-5 (12→13bp)

        sgff.history.record_operation(
            sgff.blocks, new_seq, "primerDirectedMutagenesis",
            name="test.dna",
            InputSummary={"manipulation": "replace", "val1": "4", "val2": "5"},
        )
        sgff.sequence.value = new_seq

        root = sgff.history.tree.root
        assert root.operation == "primerDirectedMutagenesis"
        assert root.seq_len == 13

        restored = self._roundtrip(sgff)
        assert restored.sequence.value == new_seq

    # --- Chaining operations ---

    def test_multiple_operations_chain(self):
        """Chain: invalid → insertFragment → changeMethylation → replace"""
        sgff = self._make_sgff("ATCGATCG")

        # 1. Insert fragment
        sgff.set_sequence("ATCGGGATCCATCG", operation="insertFragment")
        assert len(sgff.history.tree) == 2

        # 2. Change methylation (no sequence change)
        seq = sgff.sequence.value
        sgff.history.record_operation(
            sgff.blocks, seq, "changeMethylation", name="test.dna",
        )
        assert len(sgff.history.tree) == 3

        # 3. Point mutation
        sgff.set_sequence("ATCGGGATGCATCG", operation="replace")
        assert len(sgff.history.tree) == 4

        # Verify tree structure
        root = sgff.history.tree.root
        assert root.operation == "replace"
        child = root.children[0]
        assert child.operation == "changeMethylation"
        grandchild = child.children[0]
        assert grandchild.operation == "insertFragment"
        leaf = grandchild.children[0]
        assert leaf.operation == "invalid"

        # Verify all nodes have snapshots except root and leaf
        assert len(sgff.history.nodes) == 3  # nodes 1, 2, 3

        # Full roundtrip
        restored = self._roundtrip(sgff)
        assert len(restored.history.tree) == 4
        assert restored.sequence.value == "ATCGGGATGCATCG"
        assert restored.history.tree.root.operation == "replace"

    def test_multiple_operations_roundtrip_stability(self):
        """Double roundtrip: output is stable after 2 passes"""
        from sgffp.reader import SgffReader
        from sgffp.writer import SgffWriter

        sgff = self._make_sgff("ATCGATCG")
        sgff.set_sequence("ATCGATCGGG", operation="insertFragment")
        sgff.set_sequence("ATCGATCGGGAAA", operation="replace")

        bytes1 = SgffWriter.to_bytes(sgff)
        restored1 = SgffReader.from_bytes(bytes1)
        bytes2 = SgffWriter.to_bytes(restored1)

        assert bytes1 == bytes2

    # --- Real file operations ---

    def test_record_operation_on_real_file(self, pib2_dna):
        """Add operation to pIB2 and roundtrip"""
        from sgffp.reader import SgffReader
        from sgffp.writer import SgffWriter

        sgff = SgffReader.from_file(pib2_dna)
        original_count = len(sgff.history.tree)
        original_nodes = len(sgff.history.nodes)
        seq = sgff.sequence.value

        sgff.history.record_operation(
            sgff.blocks, seq, "changeMethylation",
            name=sgff.history.tree.root.name,
        )

        assert len(sgff.history.tree) == original_count + 1
        assert len(sgff.history.nodes) == original_nodes + 1

        written = SgffWriter.to_bytes(sgff)
        restored = SgffReader.from_bytes(written)

        assert len(restored.history.tree) == original_count + 1
        assert len(restored.history.nodes) == original_nodes + 1
        assert restored.history.tree.root.operation == "changeMethylation"

    def test_chain_operations_on_real_file(self, pib2_dna):
        """Chain two operations on pIB2"""
        from sgffp.reader import SgffReader
        from sgffp.writer import SgffWriter

        sgff = SgffReader.from_file(pib2_dna)
        original_count = len(sgff.history.tree)
        seq = sgff.sequence.value

        # 1. Methylation change
        sgff.history.record_operation(
            sgff.blocks, seq, "changeMethylation",
            name=sgff.history.tree.root.name,
        )

        # 2. Point mutation
        new_seq = seq[:100] + "T" + seq[101:]
        sgff.set_sequence(new_seq, operation="replace")

        assert len(sgff.history.tree) == original_count + 2

        written = SgffWriter.to_bytes(sgff)
        restored = SgffReader.from_bytes(written)
        assert len(restored.history.tree) == original_count + 2
        assert restored.sequence.value == new_seq
