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
    SgffPrimer,
    SgffPrimerList,
    SgffNotes,
    SgffProperties,
    SgffAlignment,
    SgffAlignmentList,
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
            "node_info": {"key": "value"},
        }
        node = SgffHistoryNode.from_dict(data)
        assert node.index == 5
        assert node.sequence == "ATCG"
        assert node.sequence_type == 0
        assert node.info == {"key": "value"}

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
